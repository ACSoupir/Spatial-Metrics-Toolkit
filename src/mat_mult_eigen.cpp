#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd mat_mult_eigen(const Eigen::MatrixXd& A) {
  return A * A;
}