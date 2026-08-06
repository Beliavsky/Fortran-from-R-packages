#ifndef MY_CLASS_H
#define MY_CLASS_H

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
class MyClass
{
public:

  double add(double x, double y);

  Eigen::MatrixXd scMlt(double c, const Eigen::MatrixXd& mtx);
};


#endif
