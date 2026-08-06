#include "GlmNetCpp.h"
#include "GlmNetCvCpp.h"
#include "MyClass.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <tuple>

// IMPORTANT:  Need to set R environment setting in order to build
// with C++14 features: Sys.setenv( "PKG_CXXFLAGS"="-std=c++14" )

// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// Remark:  using namespace Rcpp;  // Not required for RcppEigen; poor practice as
// well to expose the entire namespace; should be avoided!

// Enable C++14 via this plugin
using std::tuple;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;      // Still need to scope in function below
using Eigen::VectorXd;      // Still need to scope in function below
using Eigen::Map;           // Still need to scope in function below


// These 1st functions -- scalarMultiplication(.) and addReals(.) call functions
// in Dan's toy C++ examples, as a test of the Rcpp interface with Eigen.
// They are to be removed before submission of the code to Google.
// [[Rcpp::export]]
Eigen::MatrixXd scalarMultiplication(double c, const Eigen::Map<Eigen::MatrixXd>& M) {
  MyClass mc;
  return mc.scMlt(c, M);    // Same as return c * M in R

}

// [[Rcpp::export]]
double addReals(double x, double y) {
  MyClass mc;
  return mc.add(x, y);    // Same as return

}

// Actual Package functions:
// [[Rcpp::export]]
Eigen::MatrixXd fitGlmFixed(const Eigen::MatrixXd& predictor_matrix,
                            const Eigen::VectorXd& response_vector,
                            double alpha = 1, int num_lambda = 100, int glm_type = 1,
                            int max_iter = 100,
                            double abs_tol = 1.0e-4,
                            double rel_tol = 1.0e-2,
                            bool normalize_grad = false,
                            int k_fold = 5)
{
  GlmNetCvCpp gnc(predictor_matrix, response_vector);

  Eigen::MatrixXd X = gnc.FitGlmFixed();
  return X;
}

// [[Rcpp::export]]  
Eigen::MatrixXd fitGlmCv(const Eigen::MatrixXd& predictor_matrix,
                            const Eigen::VectorXd& response_vector,
                            double alpha = 1, int num_lambda = 100, int glm_type = 1,
                            int max_iter = 100,
                            double abs_tol = 1.0e-4,
                            double rel_tol = 1.0e-2,
                            bool normalize_grad = false,
                            int k_fold = 5,
                            bool has_intercept = true,
                            int k_fold_iter = 5)
{
    GlmNetCvCpp gnc(predictor_matrix, response_vector,
                    alpha,
                    num_lambda,
                    glm_type,
                    max_iter,
                    abs_tol,
                    rel_tol,
                    normalize_grad,
                    k_fold,
                    has_intercept,
                    k_fold_iter);

    Eigen::MatrixXd X = gnc.FitGlmCv();
    return X;
}

// [[Rcpp::export]]
double ExpNegativeLogLikelihood_cpp(const Eigen::VectorXd& x,
    const Eigen::MatrixXd& predictor_matrix,
                         const Eigen::VectorXd& response_vector,
                         double alpha = 1, int num_lambda = 100, int glm_type = 1,
                         int max_iter = 100,
                         double abs_tol = 1.0e-4,
                         double rel_tol = 1.0e-2,
                         bool normalize_grad = false,
                         int k_fold = 5)
{
  GlmNetCpp gnc(predictor_matrix, response_vector);

  double X = gnc.ExpNegativeLogLikelihood(x);
  return X;
}

// [[Rcpp::export]]
Eigen::VectorXd GradExpNegativeLogLikelihood_cpp(const Eigen::VectorXd& x,
                                    const Eigen::MatrixXd& predictor_matrix,
                                    const Eigen::VectorXd& response_vector,
                                    double alpha = 1, int num_lambda = 100, int glm_type = 1,
                                    int max_iter = 100,
                                    double abs_tol = 1.0e-4,
                                    double rel_tol = 1.0e-2,
                                    bool normalize_grad = false,
                                    int k_fold = 5)
{
  GlmNetCpp gnc(predictor_matrix, response_vector);

  Eigen::VectorXd X = gnc.GradExpNegativeLogLikelihood(x);
  return X;
}

// [[Rcpp::export]]
Eigen::VectorXd ProxGradDescent_cpp(const Eigen::MatrixXd& predictor_matrix,
                            const Eigen::VectorXd& response_vector,
                            double lambda = 0,
                            double alpha = 1, int glm_type = 1,
                            int max_iter = 100,
                            double abs_tol = 1.0e-4,
                            double rel_tol = 1.0e-2,
                            bool normalize_grad = false,
                            int k_fold = 5)
{
  GlmNetCpp gnc(predictor_matrix, response_vector);

  Eigen::VectorXd X = gnc.ProxGradDescent(lambda);
  return X;
}





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# clAdd(52, 54)
set.seed(42)
X <- matrix(rnorm(4*4), 4, 4)

scal <- 2.0
scal * X
scalarMultiplication(scal, X)

p <- 2.0
q <- 3.0
p + q
addReals(p, q)
  */
