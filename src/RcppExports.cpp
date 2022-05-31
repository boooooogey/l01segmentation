// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// L1BinomialApproximate
NumericVector L1BinomialApproximate(NumericVector M, NumericVector C, NumericVector lambda);
RcppExport SEXP _l01segmentation_L1BinomialApproximate(SEXP MSEXP, SEXP CSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(L1BinomialApproximate(M, C, lambda));
    return rcpp_result_gen;
END_RCPP
}
// L1BinomialApproximateCondensed
List L1BinomialApproximateCondensed(NumericVector M, NumericVector C, NumericVector lambda);
RcppExport SEXP _l01segmentation_L1BinomialApproximateCondensed(SEXP MSEXP, SEXP CSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(L1BinomialApproximateCondensed(M, C, lambda));
    return rcpp_result_gen;
END_RCPP
}
// blockcoordinatedescent
List blockcoordinatedescent(const arma::mat& Yhat, const double& lambda, const arma::vec& w, const double mintimer, const double tol);
RcppExport SEXP _l01segmentation_blockcoordinatedescent(SEXP YhatSEXP, SEXP lambdaSEXP, SEXP wSEXP, SEXP mintimerSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Yhat(YhatSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type mintimer(mintimerSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(blockcoordinatedescent(Yhat, lambda, w, mintimer, tol));
    return rcpp_result_gen;
END_RCPP
}
// L0PoissonApproximate
NumericVector L0PoissonApproximate(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0PoissonApproximate(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0PoissonApproximate(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0GaussianApproximate
NumericVector L0GaussianApproximate(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0GaussianApproximate(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0GaussianApproximate(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0ExponentialApproximate
NumericVector L0ExponentialApproximate(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0ExponentialApproximate(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0ExponentialApproximate(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0BinomialApproximate
NumericVector L0BinomialApproximate(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0BinomialApproximate(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0BinomialApproximate(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0PoissonApproximateCondensed
List L0PoissonApproximateCondensed(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0PoissonApproximateCondensed(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0PoissonApproximateCondensed(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0GaussianApproximateCondensed
List L0GaussianApproximateCondensed(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0GaussianApproximateCondensed(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0GaussianApproximateCondensed(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0ExponentialApproximateCondensed
List L0ExponentialApproximateCondensed(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0ExponentialApproximateCondensed(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0ExponentialApproximateCondensed(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0BinomialApproximateCondensed
List L0BinomialApproximateCondensed(NumericVector y, NumericVector l, NumericVector w);
RcppExport SEXP _l01segmentation_L0BinomialApproximateCondensed(SEXP ySEXP, SEXP lSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0BinomialApproximateCondensed(y, l, w));
    return rcpp_result_gen;
END_RCPP
}
// L0PoissonApproximateN
NumericVector L0PoissonApproximateN(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0PoissonApproximateN(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0PoissonApproximateN(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0GaussianApproximateN
NumericVector L0GaussianApproximateN(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0GaussianApproximateN(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0GaussianApproximateN(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0ExponentialApproximateN
NumericVector L0ExponentialApproximateN(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0ExponentialApproximateN(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0ExponentialApproximateN(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0BinomialApproximateN
NumericVector L0BinomialApproximateN(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0BinomialApproximateN(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0BinomialApproximateN(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0PoissonApproximateNCondensed
List L0PoissonApproximateNCondensed(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0PoissonApproximateNCondensed(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0PoissonApproximateNCondensed(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0GaussianApproximateNCondensed
List L0GaussianApproximateNCondensed(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0GaussianApproximateNCondensed(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0GaussianApproximateNCondensed(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0ExponentialApproximateNCondensed
List L0ExponentialApproximateNCondensed(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0ExponentialApproximateNCondensed(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0ExponentialApproximateNCondensed(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0BinomialApproximateNCondensed
List L0BinomialApproximateNCondensed(NumericVector y, int N, NumericVector w);
RcppExport SEXP _l01segmentation_L0BinomialApproximateNCondensed(SEXP ySEXP, SEXP NSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0BinomialApproximateNCondensed(y, N, w));
    return rcpp_result_gen;
END_RCPP
}
// L0PoissonBreakPoint
int L0PoissonBreakPoint(NumericVector y, NumericVector w);
RcppExport SEXP _l01segmentation_L0PoissonBreakPoint(SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0PoissonBreakPoint(y, w));
    return rcpp_result_gen;
END_RCPP
}
// L0GaussianBreakPoint
int L0GaussianBreakPoint(NumericVector y, NumericVector w);
RcppExport SEXP _l01segmentation_L0GaussianBreakPoint(SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0GaussianBreakPoint(y, w));
    return rcpp_result_gen;
END_RCPP
}
// L0ExponentialBreakPoint
int L0ExponentialBreakPoint(NumericVector y, NumericVector w);
RcppExport SEXP _l01segmentation_L0ExponentialBreakPoint(SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0ExponentialBreakPoint(y, w));
    return rcpp_result_gen;
END_RCPP
}
// L0BinomialBreakPoint
int L0BinomialBreakPoint(NumericVector y, NumericVector w);
RcppExport SEXP _l01segmentation_L0BinomialBreakPoint(SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L0BinomialBreakPoint(y, w));
    return rcpp_result_gen;
END_RCPP
}
// L1PoissonApproximate
NumericVector L1PoissonApproximate(NumericVector y, NumericVector w, double lambda);
RcppExport SEXP _l01segmentation_L1PoissonApproximate(SEXP ySEXP, SEXP wSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(L1PoissonApproximate(y, w, lambda));
    return rcpp_result_gen;
END_RCPP
}
// L1PoissonApproximateCondensed
List L1PoissonApproximateCondensed(NumericVector y, NumericVector w, double lambda);
RcppExport SEXP _l01segmentation_L1PoissonApproximateCondensed(SEXP ySEXP, SEXP wSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(L1PoissonApproximateCondensed(y, w, lambda));
    return rcpp_result_gen;
END_RCPP
}
// binVector
List binVector(NumericVector y, int bin_size);
RcppExport SEXP _l01segmentation_binVector(SEXP ySEXP, SEXP bin_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type bin_size(bin_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(binVector(y, bin_size));
    return rcpp_result_gen;
END_RCPP
}
// L1GaussianApproximate
NumericVector L1GaussianApproximate(NumericVector y, NumericVector l2, Nullable<NumericVector> weights);
RcppExport SEXP _l01segmentation_L1GaussianApproximate(SEXP ySEXP, SEXP l2SEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(L1GaussianApproximate(y, l2, weights));
    return rcpp_result_gen;
END_RCPP
}
// L1GaussianApproximateCondensed
List L1GaussianApproximateCondensed(NumericVector y, NumericVector l2, Nullable<NumericVector> w);
RcppExport SEXP _l01segmentation_L1GaussianApproximateCondensed(SEXP ySEXP, SEXP l2SEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L1GaussianApproximateCondensed(y, l2, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_l01segmentation_L1BinomialApproximate", (DL_FUNC) &_l01segmentation_L1BinomialApproximate, 3},
    {"_l01segmentation_L1BinomialApproximateCondensed", (DL_FUNC) &_l01segmentation_L1BinomialApproximateCondensed, 3},
    {"_l01segmentation_blockcoordinatedescent", (DL_FUNC) &_l01segmentation_blockcoordinatedescent, 5},
    {"_l01segmentation_L0PoissonApproximate", (DL_FUNC) &_l01segmentation_L0PoissonApproximate, 3},
    {"_l01segmentation_L0GaussianApproximate", (DL_FUNC) &_l01segmentation_L0GaussianApproximate, 3},
    {"_l01segmentation_L0ExponentialApproximate", (DL_FUNC) &_l01segmentation_L0ExponentialApproximate, 3},
    {"_l01segmentation_L0BinomialApproximate", (DL_FUNC) &_l01segmentation_L0BinomialApproximate, 3},
    {"_l01segmentation_L0PoissonApproximateCondensed", (DL_FUNC) &_l01segmentation_L0PoissonApproximateCondensed, 3},
    {"_l01segmentation_L0GaussianApproximateCondensed", (DL_FUNC) &_l01segmentation_L0GaussianApproximateCondensed, 3},
    {"_l01segmentation_L0ExponentialApproximateCondensed", (DL_FUNC) &_l01segmentation_L0ExponentialApproximateCondensed, 3},
    {"_l01segmentation_L0BinomialApproximateCondensed", (DL_FUNC) &_l01segmentation_L0BinomialApproximateCondensed, 3},
    {"_l01segmentation_L0PoissonApproximateN", (DL_FUNC) &_l01segmentation_L0PoissonApproximateN, 3},
    {"_l01segmentation_L0GaussianApproximateN", (DL_FUNC) &_l01segmentation_L0GaussianApproximateN, 3},
    {"_l01segmentation_L0ExponentialApproximateN", (DL_FUNC) &_l01segmentation_L0ExponentialApproximateN, 3},
    {"_l01segmentation_L0BinomialApproximateN", (DL_FUNC) &_l01segmentation_L0BinomialApproximateN, 3},
    {"_l01segmentation_L0PoissonApproximateNCondensed", (DL_FUNC) &_l01segmentation_L0PoissonApproximateNCondensed, 3},
    {"_l01segmentation_L0GaussianApproximateNCondensed", (DL_FUNC) &_l01segmentation_L0GaussianApproximateNCondensed, 3},
    {"_l01segmentation_L0ExponentialApproximateNCondensed", (DL_FUNC) &_l01segmentation_L0ExponentialApproximateNCondensed, 3},
    {"_l01segmentation_L0BinomialApproximateNCondensed", (DL_FUNC) &_l01segmentation_L0BinomialApproximateNCondensed, 3},
    {"_l01segmentation_L0PoissonBreakPoint", (DL_FUNC) &_l01segmentation_L0PoissonBreakPoint, 2},
    {"_l01segmentation_L0GaussianBreakPoint", (DL_FUNC) &_l01segmentation_L0GaussianBreakPoint, 2},
    {"_l01segmentation_L0ExponentialBreakPoint", (DL_FUNC) &_l01segmentation_L0ExponentialBreakPoint, 2},
    {"_l01segmentation_L0BinomialBreakPoint", (DL_FUNC) &_l01segmentation_L0BinomialBreakPoint, 2},
    {"_l01segmentation_L1PoissonApproximate", (DL_FUNC) &_l01segmentation_L1PoissonApproximate, 3},
    {"_l01segmentation_L1PoissonApproximateCondensed", (DL_FUNC) &_l01segmentation_L1PoissonApproximateCondensed, 3},
    {"_l01segmentation_binVector", (DL_FUNC) &_l01segmentation_binVector, 2},
    {"_l01segmentation_L1GaussianApproximate", (DL_FUNC) &_l01segmentation_L1GaussianApproximate, 3},
    {"_l01segmentation_L1GaussianApproximateCondensed", (DL_FUNC) &_l01segmentation_L1GaussianApproximateCondensed, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_l01segmentation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
