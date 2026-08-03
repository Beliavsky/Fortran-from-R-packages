#ifndef SUBNP_STATE_H
#define SUBNP_STATE_H

#include <RcppArmadillo.h>

struct subnp_state {
    arma::vec pars;
    arma::vec lagrange_mults;
    arma::vec scaled_eval;
    arma::mat augmented_hessian;
    arma::vec lower;
    arma::vec upper;
    arma::vec ineq_lower;
    arma::vec ineq_upper;
    Rcpp::IntegerVector problem_indicators;
    double lambda;
    arma::vec scaling_factors;
    double penalty_param;
    double tol;
    double ftol;
    double min_iter;
    int trace;
    int nfeval;
    Rcpp::Function solnp_fun, solnp_gradfun, solnp_eqfun, solnp_ineqfun, solnp_eqjac, solnp_ineqjac;
    subnp_state(Rcpp::List& state); // Declaration only
};

Rcpp::List csubnp_cpp(subnp_state& s);

#endif // SUBNP_STATE_H
