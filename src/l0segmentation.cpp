#include <Rcpp.h>
#include <math.h>
#include <limits>
#include <iostream>

using namespace Rcpp;

double inf = std::numeric_limits<double>::max();
double neginf = -inf;

struct SqrtErr{
    double knot; 
	double coef1,coef2;
	double constant;
};

double segMax(const SqrtErr * f){
	return -f->coef2/(2 * f->coef1);
}

double segEval(const SqrtErr * f, const double & b){
    return (f->coef1 * b + f->coef2) * b + f->constant;
}

void segSet(SqrtErr * f, const double & knot, const double & coef1, const double & coef2, const double & constant){
    f->knot = knot;
    f->coef1 = coef1;
    f->coef2 = coef2;
    f->constant = constant;
}

void addSeg(SqrtErr * f, int & len, const double & knot, const double & coef1, const double & coef2, const double & constant, const int & max){
    if (len+1 > max){
        stop("Function buffer is not big enough. Set max_seg_length to a bigger value (Right now %d).",max);
    }
    segSet((f+len),knot,coef1,coef2,constant);
    len = len + 1;
}

void segInit(SqrtErr * f, int & len, const double & knot, const double & y, const double & w, const int & max){
    addSeg(f,len,knot,-w,2*w*y,0,max);
}

void funcAdd(SqrtErr * d, int & d_len, SqrtErr * f, int & f_len, const double & y, const double & w){
    for(int i=0; i < f_len; i++){
        (d+i)->coef1 = (f+i)->coef1 - w;
        (d+i)->coef2 = (f+i)->coef2 + 2 * w * y;
        (d+i)->constant = (f+i)->constant;
        (d+i)->knot = (f+i)->knot;
    }
    d_len = f_len;
    f_len = 0;
}

void addRange(double * ranges, const double & val, int & ii, const int & max){
    if (ii+1 > max){
        stop("Range buffer is not big enough. Set average_range_length to a bigger value (Right now %d).",max);
    }
    ranges[ii] = val;
    ii += 1;
}

bool segSolve(const SqrtErr * f, const double & t, double & left, double & right){
    double a,b,c;
    if(f->coef1 != 0) a = f->coef1;
    else a = 1;
    b = f->coef2/a;
    c = (f->constant-t)/a;
    double delta = b*b - 4 * c;
    if (delta <= 0){
        return false;
    }
    delta = sqrt(delta);
    left = ( delta - b ) / ( 2 );
    right = ( -delta - b ) / ( 2 );
    if (left > right){
        double tmp = left;
        left = right;
        right = tmp;
    }
    return true;
}

void backtrace(const double * bstar, const double * ranges, const int * range_inds, const int & range_inds_len, double * sol){
    sol[range_inds_len] = bstar[range_inds_len];
    int s, e;
    for (int i = range_inds_len-1; i > -1; i--){
        s = range_inds[i]; 
        e = range_inds[i+1];
        sol[i] = sol[i+1];
        for (int j = s; j < e; j++){
            if (ranges[2*j] <= sol[i+1] && ranges[2*j+1] >= sol[i+1]){
                sol[i] = bstar[i];
                break;
            }
        }
    }
}

void maximize(const SqrtErr * f, const int & len, double & bprime, double & mprime){
    bprime = neginf;
    mprime = neginf;
    double bcurr;
    for (int i=0; i < len-1; i++){
        bcurr = segMax(f+i);
        if (bcurr > (f+i+1)->knot){
            bcurr = (f+i+1)->knot;
        }
        if (bcurr < (f+i)->knot){
            bcurr = (f+i)->knot;
        }
        if (mprime < segEval(f+i,bcurr)){
            mprime = segEval(f+i,bcurr);
            bprime = bcurr;
        }
    }
}

void flood(const double & threshold, const SqrtErr * in, const int & in_len, SqrtErr * out, int & out_len, double * ranges, int * range_inds, int & range_inds_len, const int & max_seg_length, const int & max_range_length){
    bool underwater = true, solution_exists;
    double tol = 1e-9;
    double left,right;
    out_len = 0;
    int ii = 2*range_inds[range_inds_len]; //index for the last element of the range array
    if (segEval(in,in->knot) < threshold){
        addSeg(out,out_len,in->knot,0,0,threshold,max_seg_length);
        addRange(ranges,neginf,ii,max_range_length);
    }
    else{
        addSeg(out,out_len,in->knot,in->coef1,in->coef2,threshold,max_seg_length);
    }
    for (int i=0; i < in_len-1; i++){
        solution_exists = segSolve((in+i), threshold, left, right);
        if (underwater && solution_exists && left < (in+i+1)->knot && left > (in+i)->knot - tol){
            addSeg(out,out_len,left,(in+i)->coef1,(in+i)->coef2,(in+i)->constant,max_seg_length);
            addRange(ranges,left,ii,max_range_length);
            underwater = false;
        }
        else if (!underwater){
            addSeg(out,out_len,(in+i)->knot,(in+i)->coef1,(in+i)->coef2,(in+i)->constant,max_seg_length);
        }
        if (!underwater && solution_exists && right < (in+i+1)->knot && right > (in+i)->knot){
            addSeg(out,out_len,right,0,0,threshold,max_seg_length);
            addRange(ranges,right,ii,max_range_length);
            underwater = true;
        }
    }
    if (ii % 2 != 0){
        addRange(ranges,inf,ii,max_range_length);
    }
    range_inds[range_inds_len+1] = ii / 2;
    range_inds_len += 1;
    addSeg(out,out_len,inf,0,0,neginf,max_seg_length);
}

void printRang(const double * ranges, const int * range_inds, const int & range_inds_len){
        for(int i=0; i<range_inds_len; i++){
                    std::cout << "Range " << i << ": ";
                            for(int j = 2*range_inds[i]; j < 2*range_inds[i+1]; j++){
                                            if (ranges[j] == inf) std::cout << "inf, ";
                                                        else if (ranges[j] == neginf) std::cout << "-inf, ";
                                                                    else std::cout << ranges[j] << ", ";
                                                                            }
                                    std::cout << std::endl;
                                        }
}

// [[Rcpp::export]]
NumericVector L0SqrtErrSeg(NumericVector y, NumericVector l2, Nullable<NumericVector> w = R_NilValue, int max_seg_length=3000, int average_range_length=6) {
    int N = y.size();

    if (l2.size() == 1){
        l2 = rep(l2[0],N-1);
    }
    else if (l2.size() != N-1){
        stop("Wrong number of lambda values given: It is %d, should have been %d.\n",l2.size(),N-1);
    }
    else if (is_true(any(l2 < 0))){
        stop("lambda values should be non-negative.\n");
    }

    NumericVector weights;
    if (w.isNotNull()){
        weights = w;
        if (weights.size() != N){
            stop("Wrong number of weights given: It is %d, should have been %d.\n",weights.size(),N);
        }
        if (is_true(any(weights <= 0))){
            stop("Weights should be positive.\n");
        }
    }
    else{
        weights = rep(1,N);
    }

    SqrtErr * d = new SqrtErr[max_seg_length];
    SqrtErr * f = new SqrtErr[max_seg_length];
    int f_len = 0, d_len = 0;

    int max_range_length = average_range_length*(N+1);
    double * ranges = new double[max_range_length];
    int * range_inds = new int[N];
    int range_inds_len = 0;
    range_inds[0] = 0;
    
    double mprime;
    double * bstar = new double[N];
    double * sol = new double[N];

    segInit(d, d_len, neginf, y[0], weights[0], max_seg_length);
    addSeg(d, d_len, inf, 0, 0, neginf, max_seg_length);

    for(int i = 1; i < N; i++){
        maximize(d, d_len, bstar[i-1], mprime);
        flood(mprime-l2[i-1], d, d_len, f, f_len, ranges, range_inds, range_inds_len, max_seg_length, max_range_length);
        funcAdd(d, d_len, f, f_len, y[i], weights[i]);
    }
    printRang(ranges,range_inds,range_inds_len);

    maximize(d, d_len, bstar[N-1], mprime);
    backtrace(bstar, ranges, range_inds, range_inds_len, sol);
    NumericVector z = NumericVector(sol,sol+N);

    //Clean up
    delete[] d;
    delete[] f;
    delete[] ranges;
    delete[] range_inds;
    delete[] bstar;
    delete[] sol;
    return z;
}

















