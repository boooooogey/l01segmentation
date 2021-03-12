// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
// L0PoisErrSeg
List L0PoisErrSeg(NumericVector y, NumericVector l2, Nullable<NumericVector> w, int max_seg_length, int average_range_length);
RcppExport SEXP _l01segmentation_L0PoisErrSeg(SEXP ySEXP, SEXP l2SEXP, SEXP wSEXP, SEXP max_seg_lengthSEXP, SEXP average_range_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type max_seg_length(max_seg_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type average_range_length(average_range_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(L0PoisErrSeg(y, l2, w, max_seg_length, average_range_length));
    return rcpp_result_gen;
END_RCPP
}
// L0PoisBreakPoints
List L0PoisBreakPoints(NumericVector y, NumericVector l2, Nullable<NumericVector> w, int max_seg_length, int average_range_length);
RcppExport SEXP _l01segmentation_L0PoisBreakPoints(SEXP ySEXP, SEXP l2SEXP, SEXP wSEXP, SEXP max_seg_lengthSEXP, SEXP average_range_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type max_seg_length(max_seg_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type average_range_length(average_range_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(L0PoisBreakPoints(y, l2, w, max_seg_length, average_range_length));
    return rcpp_result_gen;
END_RCPP
}
// L0SqrtErrSeg
List L0SqrtErrSeg(NumericVector y, NumericVector l2, Nullable<NumericVector> w, int max_seg_length, int average_range_length);
RcppExport SEXP _l01segmentation_L0SqrtErrSeg(SEXP ySEXP, SEXP l2SEXP, SEXP wSEXP, SEXP max_seg_lengthSEXP, SEXP average_range_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type max_seg_length(max_seg_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type average_range_length(average_range_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(L0SqrtErrSeg(y, l2, w, max_seg_length, average_range_length));
    return rcpp_result_gen;
END_RCPP
}
// L0SqrtErrBreakPoints
List L0SqrtErrBreakPoints(NumericVector y, NumericVector l2, Nullable<NumericVector> w, int max_seg_length, int average_range_length);
RcppExport SEXP _l01segmentation_L0SqrtErrBreakPoints(SEXP ySEXP, SEXP l2SEXP, SEXP wSEXP, SEXP max_seg_lengthSEXP, SEXP average_range_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type max_seg_length(max_seg_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type average_range_length(average_range_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(L0SqrtErrBreakPoints(y, l2, w, max_seg_length, average_range_length));
    return rcpp_result_gen;
END_RCPP
}
// L1SqrtErrFil
NumericVector L1SqrtErrFil(NumericVector y, NumericVector l2, Nullable<NumericVector> weights);
RcppExport SEXP _l01segmentation_L1SqrtErrFil(SEXP ySEXP, SEXP l2SEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(L1SqrtErrFil(y, l2, weights));
    return rcpp_result_gen;
END_RCPP
}
// L1SqrtErrBreakPoints
List L1SqrtErrBreakPoints(NumericVector y, NumericVector l2, Nullable<NumericVector> w);
RcppExport SEXP _l01segmentation_L1SqrtErrBreakPoints(SEXP ySEXP, SEXP l2SEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(L1SqrtErrBreakPoints(y, l2, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_l01segmentation_blockcoordinatedescent", (DL_FUNC) &_l01segmentation_blockcoordinatedescent, 5},
    {"_l01segmentation_L0PoisErrSeg", (DL_FUNC) &_l01segmentation_L0PoisErrSeg, 5},
    {"_l01segmentation_L0PoisBreakPoints", (DL_FUNC) &_l01segmentation_L0PoisBreakPoints, 5},
    {"_l01segmentation_L0SqrtErrSeg", (DL_FUNC) &_l01segmentation_L0SqrtErrSeg, 5},
    {"_l01segmentation_L0SqrtErrBreakPoints", (DL_FUNC) &_l01segmentation_L0SqrtErrBreakPoints, 5},
    {"_l01segmentation_L1SqrtErrFil", (DL_FUNC) &_l01segmentation_L1SqrtErrFil, 3},
    {"_l01segmentation_L1SqrtErrBreakPoints", (DL_FUNC) &_l01segmentation_L1SqrtErrBreakPoints, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_l01segmentation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
