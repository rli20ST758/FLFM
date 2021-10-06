// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]


#include <RcppArmadillo.h>
#include <RcppEigen.h>


// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
  arma::mat C = A * B;
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}



// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP eigenMapMatMult3(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
  Eigen::MatrixXd D = A * B * C;
  return Rcpp::wrap(D);
}
