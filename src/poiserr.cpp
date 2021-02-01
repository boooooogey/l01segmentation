// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <math.h>
#include <limits>
#include <iostream>
#include <boost/math/special_functions/lambert_w.hpp>

using namespace Rcpp;

double inf = std::numeric_limits<double>::max();
double neginf = -inf;

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

void addRange(double * ranges, const double & val, int & ii, const int & max){
    if (ii+1 > max){
        stop("Range buffer is not big enough. Set average_range_length to a bigger value (Possibly to some value between 5-7).");
    }
    ranges[ii] = val;
    ii += 1;
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

void printSeg(const PoisErr & f){
        if (f.knot == inf) std::cout << "x = inf, ";
            else if (f.knot == neginf) std::cout << "x = -inf, ";
                else std::cout << "x = " << f.knot << ", ";
                    if (f.constant == inf) std::cout << "c = inf, ";
                        else if (f.constant == neginf) std::cout << "c = -inf, ";
                            else std::cout << "c = " << f.constant << ", ";
                                if (f.coef1 == inf) std::cout << "coefficient1 = inf, ";
                                    else if (f.coef1 == neginf) std::cout << "coefficient1 = -inf, ";
                                        else std::cout << "coefficient1 = " << f.coef1 << ", ";
                                            if (f.coef2 == inf) std::cout << "coefficient2 = inf, ";
                                                else if (f.coef2 == neginf) std::cout << "coefficient2 = -inf, ";
                                                    else std::cout << "coefficient2 = " << f.coef2 << std::endl;
}

void printFunc(const PoisErr * f, const int & len){
        for(int i = 0; i < len-1; i++){
                    std::cout << i << ": ";
                            printSeg(f[i]);
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
        addRange(ranges,0,ii,max_range_length);
    }
    else{
        addSeg(out,out_len,in->knot,in->coef1,in->coef2,threshold,max_seg_length);
    }
    for (int i=0; i < in_len-1; i++){
        solution_exists = segSolve((in+i), threshold, left, right);
        //if (!solution_exists){
        //    std::cout << "step = " << i << std::endl;
        //    std::cout << "water level = " << threshold << std::endl;
        //    printFunc(in, in_len);
        //}
        if (underwater && solution_exists && left < (in+i+1)->knot && left > (in+i)->knot - tol){
            addSeg(out,out_len,left,(in+i)->coef1,(in+i)->coef2,(in+i)->constant,max_seg_length);
            addRange(ranges,left,ii,max_range_length);
            underwater = false;
        }
        else if (!underwater){
            addSeg(out,out_len,(in+i)->knot,(in+i)->coef1,(in+i)->coef2,(in+i)->constant,max_seg_length);
        }
        if (!underwater && solution_exists && right < (in+i+1)->knot && right > (in)->knot){
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
        //std::cout << i << ": " << std::endl;
        //std::cout << "d: " << std::endl;
        //printFunc(d,d_len);
        maximize(d, d_len, bstar[i-1], mprime);
        //std::cout << "maximum = (" << bstar[i-1] << ", " << mprime << ")" << std::endl;
        flood(mprime-l2[i-1], d, d_len, f, f_len, ranges, range_inds, range_inds_len, max_seg_length, max_range_length);
        //std::cout << "f: " << std::endl;
        //printFunc(f,f_len);
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

//IntegerVector L0PoisErrBreakPoints(NumericVector y, NumericVector l2, Nullable<NumericVector> w = R_NilValue, int max_seg_length=3000, int average_range_length=4) {
//    List out  = L0PoisErrSeg(y, l2, w, max_seg_length, average_range_length);
//    NumericVector z = out["x"];
//    IntegerVector ii(z.size());
//    int N = 0;
//    for(int i = 0; i < z.size()-1; i++){
//        if( z[i] != z[i+1]){
//            ii[N] = i+1;
//            N++;
//        }
//    }
//    return IntegerVector(ii.begin(),ii.begin()+N);
//}
//
//double clip(const double & x, const double & up, const double & low){
//    if (x > up) return up;
//    else if (x < low) return low;
//    else return x;
//}
//
//void backtrace(double * x, const double * ups, const double * lows, const int & len){
//    for(int i=len-2; i > -1; i--){
//        x[i] = clip(x[i+1],ups[i],lows[i]);
//    }
//}
//
//NumericVector L1PoisErrFil(NumericVector y, NumericVector l2, Nullable<NumericVector> weights = R_NilValue) {
//    int N = y.size();
//
//    if (l2.size() == 1){
//        l2 = rep(l2[0],N-1);
//    }
//    else if (l2.size() != N-1){
//        stop("Wrong number of lambda values given: There are %d, should have been %d.\n",l2.size(),N-1);
//    }
//    else if (is_true(any(l2 < 0))){
//        stop("lambda values should be non-negative.\n");
//    }
//
//    NumericVector w;
//    if (weights.isNotNull()){
//        w = weights;
//        if (w.size() != N){
//            stop("Wrong number of weights given: There are %d, should have been %d.\n",w.size(),N);
//        }
//        if (is_true(any(w <= 0))){
//            stop("Weights should be positive.\n");
//        }
//    }
//    else{
//        w = rep(1,N);
//    }
//
//    SqrtErr * f = new SqrtErr[N*2];
//    int max_seg_length = N*2;
//    SqrtErr first, last, low, high;
//    double * ups = new double[(N-1)*2];
//    double * lows = ups+N-1;
//    NumericVector z(N);
//    int l = N-1, u = N, up, lo;
//    lows[0] = -l2[0] + w[0] * y[0];
//    ups[0] = l2[0] + w[1] * y[0];
//    insertSeg(f, l, lows[0],  w[0], -w[0] * y[0]+l2[0], 0, max_seg_length);
//    insertSeg(f, u,  ups[0], -w[0],  w[0] * y[0]+l2[0], 0, max_seg_length);
//    segSet(&first, 0,  w[1], -l2[0]- w[1] * y[1], 0);
//    segSet( &last, 0, -w[1], -l2[0]+ w[1] * y[1], 0);
//    for(int i = 1; i < N-1; i++){
//        low = first;
//        lo = l;
//        while(segDerivative(&low, (f+lo)->knot, 1) < -l2[i]){
//            if (lo > u) break;
//            segAdd(&low, f+lo);
//            lo++;
//        }
//        lows[i] = linSolve(&low, -l2[i], 1);
//
//        high = last;
//        up = u;
//        while(segDerivative(&high, (f+up)->knot, -1) > l2[i]){
//            if (up < l) break;
//            segAdd(&high, f+up);
//            up--;
//        }
//        ups[i] = linSolve(&high, l2[i], -1);
//        
//        l = lo-1;
//        u = up+1;
//        insertSeg(f, l, lows[i],   low.coef1,  low.coef2+l2[i], 0, max_seg_length);
//        insertSeg(f, u,  ups[i],  high.coef1, high.coef2+l2[i], 0, max_seg_length);
//        segSet(&first, 0,  w[i+1], -l2[i]-w[i+1] * y[i+1], 0);
//        segSet( &last, 0, -w[i+1], -l2[i]+w[i+1] * y[i+1], 0);
//    }
//    low = first;
//    lo = l;
//    while(segDerivative(&low, (f+lo)->knot, 1) < 0){
//        if (lo > u) break;
//        segAdd(&low, f+lo);
//        lo++;
//    }
//    z[N-1] = linSolve(&low, 0, 1);
//    backtrace(z.begin(), ups, lows, z.size());
//    //clean up
//    delete[] f;
//    delete[] ups;
//    return z;
//}
//
