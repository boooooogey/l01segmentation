// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <math.h>
#include "util.h"
#include <boost/math/special_functions/lambert_w.hpp>

using namespace Rcpp;

struct PoisErr{
    double knot; 
	double coef1;
	double coef2;
	double constant;
    PoisErr& operator=(const PoisErr& other){
        if (this == &other)
            return *this;
        knot = other.knot;
        coef1 = other.coef1;
        coef2 = other.coef2;
        constant = other.constant;
        return *this;
    };
};

double segMax(const PoisErr * f){
    if (f->coef1 == 0){
        stop("Maximum achieved at infinity.");
    }
	return f->coef2/(-f->coef1);
}

double segEval(const PoisErr * f, const double & b){
    return f->coef1 * b + f->coef2 * log(b) + f->constant;
}

void segSet(PoisErr * f, const double & knot, const double & coef1, const double & coef2, const double & constant){
    f->knot = knot;
    f->coef1 = coef1;
    f->coef2 = coef2;
    f->constant = constant;
}

void segAdd(PoisErr * f, const double & coef1, const double & coef2, const double &  constant){
    f->coef1 += coef1;
    f->coef2 += coef2;
    f->constant += constant;
}

void segAdd(PoisErr * f1, const PoisErr * f2){
    f1->coef1 += f2->coef1;
    f1->coef2 += f2->coef2;
    f1->constant += f2->constant;
}

void setKnot(PoisErr * f, const double & knot){f->knot = knot;}
void setCoef1(PoisErr * f, const double & coef){f->coef1 = coef;}
void setCoef2(PoisErr * f, const double & coef){f->coef2 = coef;}
void setConstant(PoisErr * f, const double & constant){f->constant = constant;}

void addSeg(PoisErr * f, int & len, const double & knot, const double & coef1, const double & coef2, const double & constant, const int & max){
    if (len+1 > max){
        stop("Function buffer is not big enough. Set max_seg_length to a bigger value (Right now %d).",max);
    }
    segSet((f+len),knot,coef1,coef2,constant);
    len = len + 1;
}

void insertSeg(PoisErr * f, const int & i, const double & knot, const double & coef1, const double & coef2, const double & constant, const int & max){
    if (i+1 > max){
        stop("Function buffer is not big enough. Set max_seg_length to a bigger value (Right now %d).",max);
    }
    if (i < 0){
        stop("Function buffer is not big enough. Set max_seg_length to a bigger value (Right now %d).",max);
    }
    segSet((f+i),knot,coef1,coef2,constant);
}

void segInit(PoisErr * f, int & len, const double & knot, const double & y, const double & w, const int & max){
    addSeg(f,len,knot,-w,w*y,0,max);
}

void funcAdd(PoisErr * d, int & d_len, PoisErr * f, int & f_len, const double & y, const double & w){
    for(int i=0; i < f_len; i++){
        (d+i)->coef1 = (f+i)->coef1 - w;
        (d+i)->coef2 = (f+i)->coef2 + w * y;
        (d+i)->constant = (f+i)->constant;
        (d+i)->knot = (f+i)->knot;
    }
    d_len = f_len;
    f_len = 0;
}

bool segSolve(const PoisErr * f, const double & t, double & left, double & right){
    using boost::math::lambert_w0;
    using boost::math::lambert_wm1;
    double a = f->coef1, b=f->coef2, c = t - f->constant; 
    if (b == 0){
        left = right = c/a;
        return false;
    }
    double in = a * exp(c/b) / b;
    if ( in < -exp(-1) || in > 0){
        return false;
    }
    left = b * lambert_wm1(in) / a;
    right = b * lambert_w0(in) / a;
    if (left > right){
        double tmp = left;
        left = right;
        right = tmp;
    }
    return true;
}

void maximize(const PoisErr * f, const int & len, double & bprime, double & mprime){
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

void flood(const double & threshold, const PoisErr * in, const int & in_len, PoisErr * out, int & out_len, double * ranges, int * range_inds, int & range_inds_len, const int & max_seg_length, const int & max_range_length){
    bool underwater = true, solution_exists;
    double tol = 1e-9;
    double left,right;
    out_len = 0;
    int ii = 2*range_inds[range_inds_len]; //index for the last element of the range array
    if (segEval(in,in->knot) < threshold){
        addSeg(out,out_len,in->knot,0,0,threshold,max_seg_length);
        if(!addRange(ranges,0,ii,max_range_length)) stop("Range buffer is not big enough. Set average_range_length to a bigger value (Possibly to some value between 5-7).");
    }
    else{
        addSeg(out,out_len,in->knot,in->coef1,in->coef2,threshold,max_seg_length);
    }
    for (int i=0; i < in_len-1; i++){
        solution_exists = segSolve((in+i), threshold, left, right);
        if (underwater && solution_exists && left < (in+i+1)->knot && left > (in+i)->knot - tol){
            addSeg(out,out_len,left,(in+i)->coef1,(in+i)->coef2,(in+i)->constant,max_seg_length);
            if(!addRange(ranges,left,ii,max_range_length)) stop("Range buffer is not big enough. Set average_range_length to a bigger value (Possibly to some value between 5-7).");
            underwater = false;
        }
        else if (!underwater){
            addSeg(out,out_len,(in+i)->knot,(in+i)->coef1,(in+i)->coef2,(in+i)->constant,max_seg_length);
        }
        if (!underwater && solution_exists && right < (in+i+1)->knot && right > (in)->knot){
            addSeg(out,out_len,right,0,0,threshold,max_seg_length);
            if(!addRange(ranges,right,ii,max_range_length)) stop("Range buffer is not big enough. Set average_range_length to a bigger value (Possibly to some value between 5-7).");
            underwater = true;
        }
    }
    if (ii % 2 != 0){
        if(!addRange(ranges,inf,ii,max_range_length)) stop("Range buffer is not big enough. Set average_range_length to a bigger value (Possibly to some value between 5-7).");
    }
    range_inds[range_inds_len+1] = ii / 2;
    range_inds_len += 1;
    addSeg(out,out_len,inf,0,0,neginf,max_seg_length);
}

// [[Rcpp::export]]
List L0PoisErrSeg(NumericVector y, NumericVector l2, Nullable<NumericVector> w = R_NilValue, int max_seg_length=3000, int average_range_length=4) {
    int N = y.size();
    bool istherezeros = false;

    if (l2.size() == 1){
        l2 = rep(l2[0],N-1);
    }
    else if (l2.size() != N-1){
        stop("Wrong number of lambda values given: There are %d, should have been %d.\n",l2.size(),N-1);
    }
    else if (is_true(any(l2 < 0))){
        stop("lambda values should be non-negative.\n");
    }

    if (is_true(any(y < 0))){
        stop("Input must be nonnegative.\n");
    }
    if (is_true(any(y == 0))){
        y = y + 1e-9;
        istherezeros = true;
    }

    NumericVector weights;
    if (w.isNotNull()){
        weights = w;
        if (weights.size() != N){
            stop("Wrong number of weights given: There are %d, should have been %d.\n",weights.size(),N);
        }
        if (is_true(any(weights <= 0))){
            stop("Weights should be positive.\n");
        }
    }
    else{
        weights = rep(1,N);
    }

    PoisErr * d = new PoisErr[max_seg_length];
    PoisErr * f = new PoisErr[max_seg_length];
    int f_len = 0, d_len = 0;

    int max_range_length = (average_range_length+1)*N;
    double * ranges = new double[max_range_length];
    int * range_inds = new int[N];
    int range_inds_len = 0;
    range_inds[0] = 0;
    
    double mprime;
    double * bstar = new double[N];
    NumericVector z = NumericVector(N);

    segInit(d, d_len, 0, y[0], weights[0], max_seg_length);
    addSeg(d, d_len, inf, 0, 0, neginf, max_seg_length);

    for(int i = 1; i < N; i++){
        maximize(d, d_len, bstar[i-1], mprime);
        flood(mprime-l2[i-1], d, d_len, f, f_len, ranges, range_inds, range_inds_len, max_seg_length, max_range_length);
        funcAdd(d, d_len, f, f_len, y[i], weights[i]);
    }

    maximize(d, d_len, bstar[N-1], mprime);
    backtrace(bstar, ranges, range_inds, range_inds_len, z.begin());
    if (istherezeros){
        z = z - 1e-9;
    }

    //Clean up
    delete[] d;
    delete[] f;
    delete[] ranges;
    delete[] range_inds;
    delete[] bstar;
    return List::create(Named("x") = z, Named("maxval") = mprime);
}

// [[Rcpp::export]]
IntegerVector L0PoisBreakPoints(NumericVector y, NumericVector l2, Nullable<NumericVector> w = R_NilValue, int max_seg_length=3000, int average_range_length=4) {
    List out  = L0PoisErrSeg(y, l2, w, max_seg_length, average_range_length);
    NumericVector z = out["x"];
    IntegerVector ii(z.size());
    int N = 0;
    for(int i = 0; i < z.size()-1; i++){
        if( z[i] != z[i+1]){
            ii[N] = i+1;
            N++;
        }
    }
    return IntegerVector(ii.begin(),ii.begin()+N);
}
