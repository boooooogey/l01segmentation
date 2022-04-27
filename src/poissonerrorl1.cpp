#include <Rcpp.h>
#include <math.h>
#include <limits>

using namespace Rcpp;

template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void derivative_poisson( const double & a, const double & b, const double & x, double & deriv ){
    deriv = a * exp(x) + b;
}

void root_poisson( const double & a, const double & b, const double & l, double & root){
    double inf = std::numeric_limits<double>::infinity();
    double NaN = std::numeric_limits<double>::quiet_NaN();
    if (l - b < 0){
        root = -inf;
    }
    else if (l - b == 0){
        root = NaN;
    }
    else{
        root = log((l-b)/a);
    }
}

void backtrace_poisson( const double * tm, const double * tp, double * theta, const int & n ){
    for(int i = n-2; i >= 0 ; i--){
        if( theta[i+1] > tp[i] ) theta[i] = tp[i];
        else if( theta[i+1] < tm[i] ) theta[i] = tm[i];
        else theta[i] = theta[i+1];
    }
}

int backtraceCondensed_poisson( const double * tm, const double * tp, double * val, int * ii, const int & n ){
    int N = 1;
    for(int i = n-2; i >= 0 ; i--){
        if( val[N-1] > tp[i] ){ 
            ii[N-1] = i;
            val[N] = tp[i];
            N++;
        }
        else if( val[N-1] < tm[i] ){
            ii[N-1] = i;
            val[N] = tm[i];
            N++;
        }
    }
    return N;
}

void CoreLoop_poisson(const double * y, const double * w, const int & n, const double & lambda, double * tm, double * tp, double & sol){
    double inf = std::numeric_limits<double>::infinity();

    double * x = new double[2 * n];
    double * a = new double[2 * n];
    double * b = new double[2 * n];

    root_poisson( w[0], -y[0] * w[0], -lambda, tm[0] );
    root_poisson( w[0], -y[0] * w[0], lambda, tp[0] );

    int l = n - 1;
    int r = n;

    x[l] = tm[0];
    x[r] = tp[0];

    a[l] = w[0];
    b[l] = -w[0] * y[0] + lambda;

    a[r] = -w[0];
    b[r] = w[0] * y[0] + lambda;

    double afirst = w[1];
    double bfirst = -w[1] * y[1] - lambda;

    double alast = -w[1];
    double blast = w[1] * y[1] - lambda;

    double alo, blo, ahi, bhi;
    int lo, hi;
    double deriv;
    for(int k = 1; k < n-1; k++){
        alo = afirst;
        blo = bfirst;
        lo = l;
        for(int i = 0; i < r - l + 1; i++){
            derivative_poisson( alo, blo, x[lo], deriv );
            if( deriv > -lambda ) break;
            alo += a[lo];
            blo += b[lo];
            lo++;
        }
        root_poisson( alo, blo, -lambda, tm[k] );

        ahi = alast;
        bhi = blast;
        hi = r;
        for(int i = 0; i < r - l + 1; i++){
            if (x[hi] != x[hi]){
                derivative_poisson( -ahi, -bhi, -inf, deriv );
            }
            else{
                derivative_poisson( -ahi, -bhi, x[hi], deriv );
            }
            if( deriv <= lambda ) break;
            ahi += a[hi];
            bhi += b[hi];
            hi--;
        }
        root_poisson( -ahi, -bhi, lambda, tp[k] );

        l = lo - 1;
        r = hi + 1;

        x[l] = tm[k];
        a[l] = alo;
        b[l] = blo + lambda;

        x[r] = tp[k];
        a[r] = ahi;
        b[r] = bhi + lambda;

        afirst = w[k+1];
        bfirst = -w[k+1] * y[k+1] - lambda;

        alast = -w[k+1];
        blast = w[k+1] * y[k+1] - lambda;
    }
    alo = afirst;
    blo = bfirst;
    lo = l;
    for(int i = 0; i < r - l + 1; i++){
        if( x[lo] != -inf ){
            derivative_poisson( alo, blo, x[lo], deriv );
            if( deriv > 0 ) break;
        }
        alo += a[lo];
        blo += b[lo];
        lo++;
    }
    root_poisson( alo, blo, 0, sol );

    delete[] x;
    delete[] a;
    delete[] b;
}

//[[Rcpp::export]]
NumericVector L1PoissonApproximate(NumericVector y, NumericVector w, double lambda){
    int n = y.size();
    NumericVector x = NumericVector(n);
    NumericVector tm = NumericVector(n-1);
    NumericVector tp = NumericVector(n-1);
    CoreLoop_poisson(y.begin(), w.begin(), n, lambda, tm.begin(), tp.begin(), x[n-1]);
    backtrace_poisson( tm.begin(), tp.begin(), x.begin(), n );
    return x;
}

//[[Rcpp::export]]
List L1PoissonApproximateCondensed(NumericVector y, NumericVector w, double lambda){
    int n = y.size();
    NumericVector val = NumericVector(n);
    IntegerVector ii = IntegerVector(n);
    NumericVector tm = NumericVector(n-1);
    NumericVector tp = NumericVector(n-1);
    CoreLoop_poisson(y.begin(), w.begin(), n, lambda, tm.begin(), tp.begin(), val[0]);
    int k = backtraceCondensed_poisson( tm.begin(), tp.begin(), val.begin(), ii.begin(), n );
    return List::create( Named("val") = NumericVector( val.begin(), val.begin() + k ), Named("ii") = IntegerVector( ii.begin(), ii.begin() + k - 1 ) );
}
