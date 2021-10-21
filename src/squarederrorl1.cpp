#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

struct SqrtErr{
    double knot; 
	double coef1,coef2;
	double constant;
    SqrtErr& operator=(const SqrtErr& other){
        if (this == &other)
            return *this;
        knot = other.knot;
        coef1 = other.coef1;
        coef2 = other.coef2;
        constant = other.constant;
        return *this;
    }
};

double segMax(const SqrtErr * f){
	return -f->coef2/(2 * f->coef1);
}

double segEval(const SqrtErr * f, const double & b){
    return (f->coef1 * b + f->coef2) * b + f->constant;
}

double segDerivative(const SqrtErr * f, const double & b, const int & sign){
    return sign * f->coef1 * b + sign * f->coef2;
}

void segSet(SqrtErr * f, const double & knot, const double & coef1, const double & coef2, const double & constant){
    f->knot = knot;
    f->coef1 = coef1;
    f->coef2 = coef2;
    f->constant = constant;
}

void segAdd(SqrtErr * f, const double & coef1, const double & coef2, const double &  constant){
    f->coef1 += coef1;
    f->coef2 += coef2;
    f->constant += constant;
}

void segAdd(SqrtErr * f1, const SqrtErr * f2){
    f1->coef1 += f2->coef1;
    f1->coef2 += f2->coef2;
    f1->constant += f2->constant;
}

void setKnot(SqrtErr * f, const double & knot){f->knot = knot;}
void setCoef1(SqrtErr * f, const double & coef1){f->coef1 = coef1;}
void setCoef2(SqrtErr * f, const double & coef2){f->coef2 = coef2;}
void setConstant(SqrtErr * f, const double & constant){f->constant = constant;}

void addSeg(SqrtErr * f, int & len, const double & knot, const double & coef1, const double & coef2, const double & constant, const int & max){
    if (len+1 > max){
        stop("Function buffer is not big enough. Set max_seg_length to a bigger value (Right now %d).",max);
    }
    segSet((f+len),knot,coef1,coef2,constant);
    len = len + 1;
}

void insertSeg(SqrtErr * f, const int & i, const double & knot, const double & coef1, const double & coef2, const double & constant, const int & max){
    if (i+1 > max){
        stop("Function buffer is not big enough. Set max_seg_length to a bigger value (Right now %d).",max);
    }
    if (i < 0){
        stop("Function buffer is not big enough. Set max_seg_length to a bigger value (Right now %d).",max);
    }
    segSet((f+i),knot,coef1,coef2,constant);
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

double linSolve(const SqrtErr *f, const double & t, const int & sign){
    // Solves t = coef1 * x + coef2
    // x = (t - coef2) /coef1
    return ( t - sign * f->coef2 ) / (sign * f->coef1);
}

double clip(const double & x, const double & up, const double & low){
    if (x > up) return up;
    else if (x < low) return low;
    else return x;
}

void backtrace(double * x, const double * ups, const double * lows, const int & len){
    for(int i=len-2; i > -1; i--){
        x[i] = clip(x[i+1],ups[i],lows[i]);
    }
}

// [[Rcpp::export]]
NumericVector L1GaussianApproximate(NumericVector y, NumericVector l2, Nullable<NumericVector> weights = R_NilValue) {
    int N = y.size();

    if (l2.size() == 1){
        l2 = rep(l2[0],N-1);
    }
    else if (l2.size() != N-1){
        stop("Wrong number of lambda values given: There are %d, should have been %d.\n",l2.size(),N-1);
    }
    else if (is_true(any(l2 < 0))){
        stop("lambda values should be non-negative.\n");
    }

    NumericVector w;
    if (weights.isNotNull()){
        w = weights;
        if (w.size() != N){
            stop("Wrong number of weights given: There are %d, should have been %d.\n",w.size(),N);
        }
        if (is_true(any(w <= 0))){
            stop("Weights should be positive.\n");
        }
    }
    else{
        w = rep(1,N);
    }

    SqrtErr * f = new SqrtErr[N*2];
    int max_seg_length = N*2;
    SqrtErr first, last, low, high;
    double * ups = new double[(N-1)*2];
    double * lows = ups+N-1;
    NumericVector z(N);
    int l = N-1, u = N, up, lo;
    lows[0] = -l2[0] + w[0] * y[0];
    ups[0] = l2[0] + w[1] * y[0];
    insertSeg(f, l, lows[0],  w[0], -w[0] * y[0]+l2[0], 0, max_seg_length);
    insertSeg(f, u,  ups[0], -w[0],  w[0] * y[0]+l2[0], 0, max_seg_length);
    segSet(&first, 0,  w[1], -l2[0]- w[1] * y[1], 0);
    segSet( &last, 0, -w[1], -l2[0]+ w[1] * y[1], 0);
    for(int i = 1; i < N-1; i++){
        low = first;
        lo = l;
        while(segDerivative(&low, (f+lo)->knot, 1) < -l2[i]){
            if (lo > u) break;
            segAdd(&low, f+lo);
            lo++;
        }
        lows[i] = linSolve(&low, -l2[i], 1);

        high = last;
        up = u;
        while(segDerivative(&high, (f+up)->knot, -1) > l2[i]){
            if (up < l) break;
            segAdd(&high, f+up);
            up--;
        }
        ups[i] = linSolve(&high, l2[i], -1);
        
        l = lo-1;
        u = up+1;
        insertSeg(f, l, lows[i],   low.coef1,  low.coef2+l2[i], 0, max_seg_length);
        insertSeg(f, u,  ups[i],  high.coef1, high.coef2+l2[i], 0, max_seg_length);
        segSet(&first, 0,  w[i+1], -l2[i]-w[i+1] * y[i+1], 0);
        segSet( &last, 0, -w[i+1], -l2[i]+w[i+1] * y[i+1], 0);
    }
    low = first;
    lo = l;
    while(segDerivative(&low, (f+lo)->knot, 1) < 0){
        if (lo > u) break;
        segAdd(&low, f+lo);
        lo++;
    }
    z[N-1] = linSolve(&low, 0, 1);
    backtrace(z.begin(), ups, lows, z.size());
    //clean up
    delete[] f;
    delete[] ups;
    return z;
}

// [[Rcpp::export]]
List L1GaussianApproximateCondensed(NumericVector y, NumericVector l2, Nullable<NumericVector> w = R_NilValue) {
    NumericVector z = L1GaussianApproximate(y, l2, w);
    IntegerVector ii(z.size());
    NumericVector v(z.size());
    int N = 0;
    v[0] = z[0];
    for(int i = 0; i < z.size()-1; i++){
        if( z[i] != z[i+1]){
            ii[N] = i+1;
            N++;
            v[N] = z[i+1];
        }
    }
    return List::create(Named("ii") = IntegerVector(ii.begin(),ii.begin()+N), Named("val") = NumericVector(v.begin(),v.begin()+N+1));
}
