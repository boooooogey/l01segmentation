# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

gmm <- function(data, k, km_iter = 10L, em_iter = 5L, var_floor = 1e-10L, verbose = FALSE) {
    .Call(`_l01segmentation_gmm`, data, k, km_iter, em_iter, var_floor, verbose)
}

blockcoordinatedescent <- function(Yhat, lambda, w, mintimer = 5, tol = 1e-6) {
    .Call(`_l01segmentation_blockcoordinatedescent`, Yhat, lambda, w, mintimer, tol)
}

L0PoisErrSeg <- function(y, l2, w = NULL, max_seg_length = 3000L, average_range_length = 4L) {
    .Call(`_l01segmentation_L0PoisErrSeg`, y, l2, w, max_seg_length, average_range_length)
}

L0PoisBreakPoints <- function(y, l2, w = NULL, max_seg_length = 3000L, average_range_length = 4L) {
    .Call(`_l01segmentation_L0PoisBreakPoints`, y, l2, w, max_seg_length, average_range_length)
}

L0SqrtErrSeg <- function(y, l2, w = NULL, max_seg_length = 3000L, average_range_length = 4L) {
    .Call(`_l01segmentation_L0SqrtErrSeg`, y, l2, w, max_seg_length, average_range_length)
}

L0SqrtErrBreakPoints <- function(y, l2, w = NULL, max_seg_length = 3000L, average_range_length = 4L) {
    .Call(`_l01segmentation_L0SqrtErrBreakPoints`, y, l2, w, max_seg_length, average_range_length)
}

L1SqrtErrFil <- function(y, l2, weights = NULL) {
    .Call(`_l01segmentation_L1SqrtErrFil`, y, l2, weights)
}

L1SqrtErrBreakPoints <- function(y, l2, w = NULL) {
    .Call(`_l01segmentation_L1SqrtErrBreakPoints`, y, l2, w)
}

