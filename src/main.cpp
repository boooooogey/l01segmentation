// [[Rcpp::depends(BH)]]
#include "binomialerror.hpp"
#include "exponentialerror.hpp"
#include "poissonerror.hpp"
#include "squarederror.hpp"
#include "piecewisefunction.hpp"
#include "range.hpp"
#include "util.hpp"
#include <Rcpp.h>
#include <set>
#include <vector>
//#include "gperftools/profiler.h"

using namespace Rcpp;

void add_to_lookup(std::set<int>& lookup, IntegerVector x) {
  int n = x.size();
  for (int i = 0; i < n; i++) {
    lookup.insert(x[i]);
  }
}

// [[Rcpp::export]]
std::set<int> combine_two_bp_sets_(IntegerVector x, IntegerVector y) {

  std::set<int> lookup;
  add_to_lookup(lookup, x);
  add_to_lookup(lookup, y);
  return lookup;
}

//// [[Rcpp::export]]
//SEXP start_profiler(SEXP str) {
//  ProfilerStart(as<const char*>(str));
//  return R_NilValue;
//}
//
//// [[Rcpp::export]]
//SEXP stop_profiler() {
//  ProfilerStop();
//  return R_NilValue;
//}

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
NumericVector L0BinomialApproximate(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size()/2;
    NumericVector x = NumericVector(n);
    approximate<BinomialError>(n, y.begin(), l.begin(), w.begin(), x.begin());
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
    return List::create(Named("value") = val[Rcpp::Range(0, k-1)] , Named("start") = start[Rcpp::Range(0, k-1)], Named("end") = end[Rcpp::Range(0, k-1)]);
}

//[[Rcpp::export]]
List L0GaussianApproximateCondensed(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    int k;
    NumericVector val = NumericVector(n);
    IntegerVector start = IntegerVector(n);
    IntegerVector end = IntegerVector(n);
    approximate<SquaredError>(n, y.begin(), l.begin(), w.begin(), k, start.begin(), end.begin(), val.begin());
    return List::create(Named("value") = val[Rcpp::Range(0, k-1)] , Named("start") = start[Rcpp::Range(0, k-1)], Named("end") = end[Rcpp::Range(0, k-1)]);
}

//[[Rcpp::export]]
List L0ExponentialApproximateCondensed(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size();
    int k;
    NumericVector val = NumericVector(n);
    IntegerVector start = IntegerVector(n);
    IntegerVector end = IntegerVector(n);
    approximate<ExponentialError>(n, y.begin(), l.begin(), w.begin(), k, start.begin(), end.begin(), val.begin());
    return List::create(Named("value") = val[Rcpp::Range(0, k-1)] , Named("start") = start[Rcpp::Range(0, k-1)], Named("end") = end[Rcpp::Range(0, k-1)]);
}

//[[Rcpp::export]]
List L0BinomialApproximateCondensed(NumericVector y, NumericVector l, NumericVector w){
    int n = y.size()/2;
    int k;
    NumericVector val = NumericVector(n);
    IntegerVector start = IntegerVector(n);
    IntegerVector end = IntegerVector(n);
    approximate<BinomialError>(n, y.begin(), l.begin(), w.begin(), k, start.begin(), end.begin(), val.begin());
    return List::create(Named("value") = val[Rcpp::Range(0, k-1)] , Named("start") = start[Rcpp::Range(0, k-1)], Named("end") = end[Rcpp::Range(0, k-1)]);
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
NumericVector L0BinomialApproximateN(NumericVector y, int N, NumericVector w){
    int n = y.size()/2;
    NumericVector x = NumericVector(n);
    approximate<BinomialError>(n, y.begin(), N, w.begin(), x.begin());
    return x;
}

//[[Rcpp::export]]
List L0PoissonApproximateNCondensed(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector val = NumericVector(N);
    IntegerVector start = IntegerVector(N);
    IntegerVector end = IntegerVector(N);
    approximate<PoissonError>(n, y.begin(), N, w.begin(), start.begin(), end.begin(), val.begin());
    return List::create(Named("value") = val , Named("start") = start, Named("end") = end);
}

//[[Rcpp::export]]
List L0GaussianApproximateNCondensed(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector val = NumericVector(N);
    IntegerVector start = IntegerVector(N);
    IntegerVector end = IntegerVector(N);
    approximate<SquaredError>(n, y.begin(), N, w.begin(), start.begin(), end.begin(), val.begin());
    return List::create(Named("value") = val , Named("start") = start, Named("end") = end);
}

//[[Rcpp::export]]
List L0ExponentialApproximateNCondensed(NumericVector y, int N, NumericVector w){
    int n = y.size();
    NumericVector val = NumericVector(N);
    IntegerVector start = IntegerVector(N);
    IntegerVector end = IntegerVector(N);
    approximate<ExponentialError>(n, y.begin(), N, w.begin(), start.begin(), end.begin(), val.begin());
    return List::create(Named("value") = val , Named("start") = start, Named("end") = end);
}

//[[Rcpp::export]]
List L0BinomialApproximateNCondensed(NumericVector y, int N, NumericVector w){
    int n = y.size()/2;
    NumericVector val = NumericVector(N);
    IntegerVector start = IntegerVector(N);
    IntegerVector end = IntegerVector(N);
    approximate<BinomialError>(n, y.begin(), N, w.begin(), start.begin(), end.begin(), val.begin());
    return List::create(Named("value") = val , Named("start") = start, Named("end") = end);
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

//[[Rcpp::export]]
int L0BinomialBreakPoint(NumericVector y, NumericVector w){
    int n = y.size()/2;
    return findbreakpoint<BinomialError>(n, y.begin(), w.begin());
}
