#include <Rcpp.h>
#include <math.h>
#include <limits>

using namespace Rcpp;

template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void derivative( const double & a, const double & b, const double & c, const double & x, double & deriv ){
    double inf = std::numeric_limits<double>::infinity();
    if( x == 0 ) deriv = -1 * inf;
    else if( x == 1 ) deriv = inf;
    else if( a == 0 && b == 0 ) deriv = c;
    else if( a == 0 ){
        if ( x == 1 ) deriv = sgn(b) * inf;
        else deriv = b / ( 1 - x ) + c;
    }
    else if( b == 0 ){
        if ( x == 0 ) deriv = -1 * sgn(a) * inf;
        else deriv = -1 * a / x + c;
    }
    else deriv = -1 * a / x + b / ( 1 - x ) + c;
}

void root( const double & a, const double & b, const double & c, const double & l, double & root){
    if( a == 0 ){
        if( l == 0 ) root = 0;
        else root = ( l - b ) / l;
    }
    else if( b == 0 ){
        if( l == 0 ) root = 1;
        else root = -1 * a / l;
    }
    else if( l - c == 0 ) root = a / (double)( a + b );
    else{
        double aeq = l - c;
        double beq = b + a + - l + c;
        double ceq = a;
        double xplus = ( -1 * beq + sqrt( pow(beq, 2) + 4 * aeq * ceq ) ) / 2 / aeq;
        double xminus = ( -1 * beq - sqrt( pow(beq, 2) + 4 * aeq * ceq ) ) / 2 / aeq;
        if( xplus >= 0 && xplus <= 1 ) root = xplus;
        else root = xminus;
    }
}

bool in01( const double & x ){
    return x < 1 && x > 0;
}

void setKnots( double & tm, double & tp, const double & M ){
    if( !( in01(tm) || in01(tp) ) ){
        if( M != 0 ){
            tm = 1;
            tp = 1;
        }
        else{
            tm = 0;
            tp = 0;
        }
    }
    else if( !in01(tm) ) tm = 0;
    else if( !in01(tp) ) tp = 1;
}

void backtrace( const double * tm, const double * tp, double * theta, const int & n ){
    for(int i = n-2; i >= 0 ; i--){
        if( theta[i+1] > tp[i] ) theta[i] = tp[i];
        else if( theta[i+1] < tm[i] ) theta[i] = tm[i];
        else theta[i] = theta[i+1];
    }
}

int backtraceCondensed( const double * tm, const double * tp, double * val, int * ii, const int & n ){
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

void CoreLoop(const double * M, const double * C, const int & n, const double & lambda, double * tm, double * tp, double & sol){

    double * x = new double[2 * n];
    double * a = new double[2 * n];
    double * b = new double[2 * n];
    double * c = new double[2 * n];

    root( M[0], C[0] - M[0], 0, -lambda, tm[0] );
    root( M[0], C[0] - M[0], 0, lambda, tp[0] );

    int l = n - 1;
    int r = n;

    setKnots( tm[0], tp[0], M[0] );

    x[l] = tm[0];
    x[r] = tp[0];

    a[l] = M[0];
    b[l] = C[0] - M[0];
    c[l] = lambda;

    a[r] = -M[0];
    b[r] = -C[0] + M[0];
    c[r] = lambda;

    double afirst = M[1];
    double bfirst = C[1] - M[1];
    double cfirst = -lambda;

    double alast = -M[1];
    double blast = -C[1] + M[1];
    double clast = -lambda;

    double alo, blo, clo, ahi, bhi, chi;
    int lo, hi;
    double deriv;
    for(int k = 1; k < n-1; k++){
        alo = afirst;
        blo = bfirst;
        clo = cfirst;
        lo = l;
        for(int i = 0; i < r - l + 1; i++){
            derivative( alo, blo, clo, x[lo], deriv );
            if( deriv > -lambda ) break;
            alo += a[lo];
            blo += b[lo];
            clo += c[lo];
            lo++;
        }
        root( alo, blo, clo, -lambda, tm[k] );

        ahi = alast;
        bhi = blast;
        chi = clast;
        hi = r;
        for(int i = 0; i < r - l + 1; i++){
            derivative( -ahi, -bhi, -chi, x[hi], deriv );
            if( deriv < lambda ) break;
            ahi += a[hi];
            bhi += b[hi];
            chi += c[hi];
            hi--;
        }
        root( -ahi, -bhi, -chi, lambda, tp[k] );

        l = lo - 1;
        r = hi + 1;

        setKnots( tm[k], tp[k], M[k] );

        x[l] = tm[k];
        a[l] = alo;
        b[l] = blo;
        c[l] = clo + lambda;

        x[r] = tp[k];
        a[r] = ahi;
        b[r] = bhi;
        c[r] = chi + lambda;

        afirst = M[k+1];
        bfirst = C[k+1] - M[k+1];
        cfirst = -lambda;

        alast = -M[k+1];
        blast = -C[k+1] + M[k+1];
        clast = -lambda;
    }
    alo = afirst;
    blo = bfirst;
    clo = cfirst;
    lo = l;
    for(int i = 0; i < r - l + 1; i++){
        derivative( alo, blo, clo, x[lo], deriv );
        if( deriv > 0 ) break;
        alo += a[lo];
        blo += b[lo];
        clo += c[lo];
        lo++;
    }
    root( alo, blo, clo, 0, sol );

    delete[] x;
    delete[] a;
    delete[] b;
    delete[] c;
}

//[[Rcpp::export]]
NumericVector L1BinomialApproximate(NumericVector M, NumericVector C, double lambda){
    int n = M.size();
    NumericVector x = NumericVector(n);
    NumericVector tm = NumericVector(n-1);
    NumericVector tp = NumericVector(n-1);
    CoreLoop(M.begin(), C.begin(), n, lambda, tm.begin(), tp.begin(), x[n-1]);
    backtrace( tm.begin(), tp.begin(), x.begin(), n );
    return x;
}

//[[Rcpp::export]]
List L1BinomialApproximateCondensed(NumericVector M, NumericVector C, double lambda){
    int n = M.size();
    NumericVector val = NumericVector(n);
    IntegerVector ii = IntegerVector(n);
    NumericVector tm = NumericVector(n-1);
    NumericVector tp = NumericVector(n-1);
    CoreLoop(M.begin(), C.begin(), n, lambda, tm.begin(), tp.begin(), val[0]);
    int k = backtraceCondensed( tm.begin(), tp.begin(), val.begin(), ii.begin(), n );
    return List::create( Named("val") = NumericVector( val.begin(), val.begin() + k ), Named("ii") = IntegerVector( ii.begin(), ii.begin() + k - 1 ) );
}
