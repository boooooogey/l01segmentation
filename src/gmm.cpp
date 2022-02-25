#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

List gmm(arma::mat const & data, int const & k, int const & km_iter = 10, int const & em_iter = 5, double const & var_floor = 1e-10, bool verbose = false) {
    arma::gmm_full model;
    int n = data.n_cols;
    bool status = model.learn(data, k, arma::eucl_dist, arma::random_subset, km_iter, em_iter, var_floor, verbose);
    if(status == false) warning("learning failed");
    arma::mat z(k,n);
    for(int i = 0; i < k; i++){
        z.row(i) = model.log_p(data, i);
    }
    return List::create(Named("centers") = model.means, Named("covariance") = model.fcovs, Named("z") = z, Named("weights") = model.hefts, Named("loglikelihood")=model.sum_log_p(data));
}
