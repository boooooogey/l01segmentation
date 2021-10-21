// [[Rcpp::depends(BH)]]
#include "squarederror.hpp"
#include "poissonerror.hpp"
#include "exponentialerror.hpp"
#include "piecewisefunction.hpp"
#include "range.hpp"
#include "util.hpp"
#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector L0PoissonApproximate(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    NumericVector x = NumericVector(n);
    approximate<PoissonError>(n, y.begin(), l.begin(), w.begin(), x.begin());
    return x;
}

//[[Rcpp::export]]
NumericVector L0GaussianApproximate(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    NumericVector x = NumericVector(n);
    approximate<SquaredError>(n, y.begin(), l.begin(), w.begin(), x.begin());
    return x;
}

//[[Rcpp::export]]
NumericVector L0ExponentialApproximate(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    NumericVector x = NumericVector(n);
    approximate<ExponentialError>(n, y.begin(), l.begin(), w.begin(), x.begin());
    return x;
}

//[[Rcpp::export]]
List L0PoissonApproximateCondensed(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    int k;
    NumericVector val = NumericVector(n);
    IntegerVector start = IntegerVector(n);
    IntegerVector end = IntegerVector(n);
    approximate<PoissonError>(n, y.begin(), l.begin(), w.begin(), k, start.begin(), end.begin(), val.begin());
    return List::create(Named("values") = val[Rcpp::Range(0, k-1)] , Named("start") = start[Rcpp::Range(0, k-1)], Named("end") = end[Rcpp::Range(0, k-1)]);
}

//[[Rcpp::export]]
List L0GaussianApproximateCondensed(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    int k;
    NumericVector val = NumericVector(n);
    IntegerVector start = IntegerVector(n);
    IntegerVector end = IntegerVector(n);
    approximate<SquaredError>(n, y.begin(), l.begin(), w.begin(), k, start.begin(), end.begin(), val.begin());
    return List::create(Named("values") = val[Rcpp::Range(0, k-1)] , Named("start") = start[Rcpp::Range(0, k-1)], Named("end") = end[Rcpp::Range(0, k-1)]);
}

//[[Rcpp::export]]
List L0ExponentialApproximateCondensed(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    int k;
    NumericVector val = NumericVector(n);
    IntegerVector start = IntegerVector(n);
    IntegerVector end = IntegerVector(n);
    approximate<ExponentialError>(n, y.begin(), l.begin(), w.begin(), k, start.begin(), end.begin(), val.begin());
    return List::create(Named("values") = val[Rcpp::Range(0, k-1)] , Named("start") = start[Rcpp::Range(0, k-1)], Named("end") = end[Rcpp::Range(0, k-1)]);
}

//[[Rcpp::export]]
NumericVector L0PoissonApproximateN(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector x = NumericVector(n);
    approximate<PoissonError>(n, y.begin(), N, w.begin(), x.begin());
    return x;
}

//[[Rcpp::export]]
NumericVector L0GaussianApproximateN(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector x = NumericVector(n);
    approximate<SquaredError>(n, y.begin(), N, w.begin(), x.begin());
    return x;
}

//[[Rcpp::export]]
NumericVector L0ExponentialApproximateN(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector x = NumericVector(n);
    approximate<ExponentialError>(n, y.begin(), N, w.begin(), x.begin());
    return x;
}

//[[Rcpp::export]]
List L0PoissonApproximateNCondensed(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector val = NumericVector(N);
    IntegerVector start = IntegerVector(N);
    IntegerVector end = IntegerVector(N);
    approximate<PoissonError>(n, y.begin(), N, w.begin(), start.begin(), end.begin(), val.begin());
    return List::create(Named("values") = val , Named("start") = start, Named("end") = end);
}

//[[Rcpp::export]]
List L0GaussianApproximateNCondensed(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector val = NumericVector(N);
    IntegerVector start = IntegerVector(N);
    IntegerVector end = IntegerVector(N);
    approximate<SquaredError>(n, y.begin(), N, w.begin(), start.begin(), end.begin(), val.begin());
    return List::create(Named("values") = val , Named("start") = start, Named("end") = end);
}

//[[Rcpp::export]]
List L0ExponentialApproximateNCondensed(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector val = NumericVector(N);
    IntegerVector start = IntegerVector(N);
    IntegerVector end = IntegerVector(N);
    approximate<ExponentialError>(n, y.begin(), N, w.begin(), start.begin(), end.begin(), val.begin());
    return List::create(Named("values") = val , Named("start") = start, Named("end") = end);
}

//[[Rcpp::export]]
int L0PoissonBreakPoint(NumericVector y, NumericVector w){
    int n = y.size();
    return findbreakpoint<PoissonError>(n, y.begin(), w.begin());
}

//[[Rcpp::export]]
int L0GaussianBreakPoint(NumericVector y, NumericVector w){
    int n = y.size();
    return findbreakpoint<SquaredError>(n, y.begin(), w.begin());
}

//[[Rcpp::export]]
int L0ExponentialBreakPoint(NumericVector y, NumericVector w){
    int n = y.size();
    return findbreakpoint<ExponentialError>(n, y.begin(), w.begin());
}
