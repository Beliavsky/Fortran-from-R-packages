#include "MyClass.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;

double MyClass::add(double x, double y)
{
  double z = x + y;
  return z;
}

MatrixXd MyClass::scMlt(double c, const MatrixXd& mtx)
{
  return c * mtx;
}
