
 #include <RcppEigen.h>

 
 // [[Rcpp::depends(RcppEigen)]]
 
 using Eigen::Map;                       // 'maps' rather than copies
 using Eigen::MatrixXd;                  // variable size matrix, double precision
 using Eigen::VectorXd;                  // variable size vector, double precision
 using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
 using namespace Rcpp;
 using namespace Eigen;
 
 // [[Rcpp::export]]
 List getEigenDecomp(Map<MatrixXd> M) {
   SelfAdjointEigenSolver<MatrixXd> es(M);
   return List::create(Named("values") = es.eigenvalues(),
                        Named("vectors") = es.eigenvectors());
 }
 
 
