#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

// Beta-binomial IRLS kernels for BARCS.
//
// These two functions were previously compiled inside the CB2 package and
// reached through its shared library.  They are BARCS's own inner loop, so
// they live here; CB2 remains a dependency for the original two-group
// workflow and for quant(), which needs its AdaptiveHash reference index.

//' Weighted least-squares cross-products for beta-binomial IRLS
//'
//' Computes the symmetric information matrix \eqn{X^\top W X} and score
//' target \eqn{X^\top W z} in compiled code. This is an internal kernel used
//' by `bbreg()`.
//'
//' @param x Numeric design matrix.
//' @param weight Positive working-weight vector.
//' @param response Numeric working-response vector.
//' @return A list with `information` and `score_target`.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List bb_wls_system_cpp(const arma::mat &x,
                             const arma::vec &weight,
                             const arma::vec &response) {
  if (x.n_rows != weight.n_elem || x.n_rows != response.n_elem) {
    Rcpp::stop("The design, weight, and response dimensions do not match.");
  }
  if (!weight.is_finite() || arma::any(weight <= 0) ||
      !response.is_finite() || !x.is_finite()) {
    Rcpp::stop("The weighted least-squares inputs must be finite and weights positive.");
  }
  const arma::mat weighted_x = x.each_col() % weight;
  const arma::mat information = x.t() * weighted_x;
  const arma::vec score_target = x.t() * (weight % response);
  return Rcpp::List::create(
    Rcpp::_["information"] = information,
    Rcpp::_["score_target"] = score_target
  );
}

//' Compiled weighted least-squares solve for beta-binomial IRLS
//'
//' Forms \eqn{X^\top W X} and \eqn{X^\top W z}, solves the symmetric system,
//' and optionally returns its inverse in one compiled call.
//'
//' @param x Numeric design matrix.
//' @param weight Positive working-weight vector.
//' @param response Numeric working-response vector.
//' @param covariance Whether to return the inverse information matrix.
//' @return A list with `coefficient` and, when requested, `covariance`.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List bb_wls_solve_cpp(const arma::mat &x,
                            const arma::vec &weight,
                            const arma::vec &response,
                            const bool covariance = false) {
  if (x.n_rows != weight.n_elem || x.n_rows != response.n_elem) {
    Rcpp::stop("The design, weight, and response dimensions do not match.");
  }
  const arma::mat weighted_x = x.each_col() % weight;
  const arma::mat information = x.t() * weighted_x;
  const arma::vec score_target = x.t() * (weight % response);
  arma::vec coefficient;
  const bool solved = arma::solve(
    coefficient, information, score_target,
    arma::solve_opts::likely_sympd
  );
  if (!solved || !coefficient.is_finite()) {
    Rcpp::stop("The weighted least-squares system is singular.");
  }
  if (!covariance) {
    return Rcpp::List::create(Rcpp::_["coefficient"] = coefficient);
  }
  arma::mat inverse;
  if (!arma::inv_sympd(inverse, information)) {
    Rcpp::stop("The weighted information matrix is not positive definite.");
  }
  return Rcpp::List::create(
    Rcpp::_["coefficient"] = coefficient,
    Rcpp::_["covariance"] = inverse
  );
}
