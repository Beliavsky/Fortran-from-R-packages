#ifndef SUBNP_HELPERS_H
#define SUBNP_HELPERS_H

#include <RcppArmadillo.h>
#include <utility>
#include <limits>
#include <string>


inline void print_progress(int iter, double obj, double constr_norm, double rel_obj_change,
                           double step_norm, double penalty_param
) {
    Rcpp::Rcout << "Iter: " << std::setw(4) << iter
                << "  Obj: " << std::setprecision(6) << std::setw(8) << obj
                << "  ||Constr||: " << std::scientific << std::setprecision(3) << constr_norm
                << "  RelObj: " << std::scientific << std::setprecision(3) << rel_obj_change
                << "  Step: " << std::scientific << std::setprecision(2) << step_norm
                << "  Penalty: " << std::scientific << std::setprecision(2) << penalty_param
                << std::endl;
}

// Norm
inline double vnorm(const arma::vec& v) {
    return std::sqrt(dot(v, v));
}

// Solver Warnings
inline void solver_warnings(const std::string& message_indicator) {
    std::string ans;
    if (message_indicator == "M1") {
        ans = "\nsolnp: Redundant constraints were found.";
    } else if (message_indicator == "M2") {
        ans = "\nLinearized problem has no feasible solution. The problem may not be feasible.";
    } else if (message_indicator == "M3") {
        ans = "\nMinor optimization routine did not converge in the specified number of minor iterations.";
    }
    if (!ans.empty()) {
        Rcpp::warning(ans);
    }
}

// Augmented Lagrangian
inline double augmented_lagrangian(const arma::vec& scaled_value,
                                   const arma::vec& lagrange_mults,
                                   double rho, arma::uword n_constraints)
{
    arma::vec ceval = scaled_value.subvec(1, n_constraints); // elements 2:(n_constraints+1) in R
    return scaled_value(0) - arma::dot(lagrange_mults, ceval) + rho * arma::dot(ceval, ceval);
}

// Weighted system with QR
inline std::pair<arma::vec, bool> qr_solve_weighted_system(const arma::mat& augmented_jacobian,
                                                           const arma::vec& dx, const arma::rowvec& cx)
{
    arma::mat Dx = arma::diagmat(dx);
    arma::mat At = (augmented_jacobian * Dx).t();
    arma::vec rhs = dx % cx.t();

    arma::vec y;
    bool success = false;
    try {
        arma::mat Q, R;
        arma::qr(Q, R, At);
        arma::vec Qt_rhs = Q.t() * rhs;
        y = arma::solve(R, Qt_rhs, arma::solve_opts::fast + arma::solve_opts::no_approx);
        success = true;
    } catch (const std::runtime_error&) {
        //y.reset();
        y.zeros();
        success = false;
    }
    return std::make_pair(y, success);
}

// Transposed system with QR
inline std::pair<arma::vec, bool>
    qr_solve_transposed_system(const arma::mat& cz,
                               const arma::mat& augmented_jacobian,
                               const arma::vec& delta_gradient)
    {
        arma::vec y;
        bool success = false;
        try {
            arma::mat At = cz.t() * augmented_jacobian.t();
            arma::mat Q, R;
            arma::qr(Q, R, At);
            arma::vec Qt_yg = Q.t() * delta_gradient;
            y = arma::solve(R, Qt_yg, arma::solve_opts::fast + arma::solve_opts::no_approx);
            success = true;
        } catch (const std::runtime_error&) {
            //y.reset();
            y.zeros();
            success = false;
        }
        return std::make_pair(y, success);
    }

// Weighted system with solve
inline std::pair<arma::vec, bool> solve_weighted_system(const arma::mat& augmented_jacobian,
                                                        const arma::vec& dx, const arma::rowvec& cx)
{
    arma::mat Dx = arma::diagmat(dx);
    arma::mat At = (augmented_jacobian * Dx).t();
    arma::vec rhs = dx % cx.t();

    arma::vec y;
    bool success = false;
    try {
        y = arma::solve(At, rhs, arma::solve_opts::fast + arma::solve_opts::no_approx);
        success = true;
    } catch (const std::runtime_error&) {
        //y.reset();
        y.zeros();
        success = false;
    }
    return std::make_pair(y, success);
}

// Cholesky
inline std::pair<arma::mat, bool> cholesky(const arma::mat& A) {
    arma::mat cz;
    bool chol_success = false;
    try {
        cz = arma::chol(A);
        chol_success = true;
    } catch (const std::runtime_error&) {
        cz.reset();
        cz.zeros();
        chol_success = false;
    }
    return std::make_pair(cz, chol_success);
}


inline Rcpp::List compute_kkt_diagnostics(
        const arma::vec& parameters,
        const arma::vec& lagrange_mults,
        int n_eq,
        int n_ineq,
        const Rcpp::Function& gradient_fun,
        const Rcpp::Function& eq_j,
        const Rcpp::Function& ineq_j,
        const Rcpp::Function& eq_f,
        const Rcpp::Function& ineq_f,
        const arma::vec& ineq_lower,
        const arma::vec& ineq_upper,
        double tol,
        int error_code
) {
    double NAval = NA_REAL;
    if (error_code != 0) {
        // Compute only primal constraint violations
        double eq_violation = 0.0, ineq_violation = 0.0;
        if (n_eq > 0)
            eq_violation = arma::norm(Rcpp::as<arma::vec>(eq_f(parameters)), "inf");
        if (n_ineq > 0) {
            arma::vec ineq_vals = Rcpp::as<arma::vec>(ineq_f(parameters));
            arma::vec lower_viol = arma::clamp(ineq_lower - ineq_vals, 0.0, arma::datum::inf);
            arma::vec upper_viol = arma::clamp(ineq_vals - ineq_upper, 0.0, arma::datum::inf);
            arma::vec all_viol = arma::max(lower_viol, upper_viol);
            ineq_violation = all_viol.max();
        }
        return Rcpp::List::create(
            Rcpp::_["kkt_stationarity"] = NAval,
            Rcpp::_["eq_violation"] = eq_violation,
            Rcpp::_["ineq_violation"] = ineq_violation,
            Rcpp::_["dual_feas_violation"] = NAval,
            Rcpp::_["compl_slackness"] = NAval
        );
    }
    arma::vec grad = Rcpp::as<arma::vec>(gradient_fun(parameters));

    arma::vec lagr_mult_eq, lagr_mult_ineq;
    if (n_eq > 0) lagr_mult_eq = lagrange_mults.subvec(0, n_eq - 1);
    if (n_ineq > 0) lagr_mult_ineq = lagrange_mults.subvec(n_eq, n_eq + n_ineq - 1);

    arma::mat A_eq, A_ineq;
    if (n_eq > 0)   A_eq   = Rcpp::as<arma::mat>(eq_j(parameters));
    if (n_ineq > 0) A_ineq = Rcpp::as<arma::mat>(ineq_j(parameters));

    arma::vec lagr_grad = grad;
    if (n_eq > 0)   lagr_grad -= A_eq.t() * lagr_mult_eq;
    if (n_ineq > 0) lagr_grad -= A_ineq.t() * lagr_mult_ineq;

    double kkt_stationarity = arma::norm(lagr_grad, "inf");

    double eq_violation = 0.0, ineq_violation = 0.0;
    arma::vec all_viol;
    if (n_ineq > 0) {
        arma::vec ineq_vals = Rcpp::as<arma::vec>(ineq_f(parameters));
        arma::vec lower_viol = arma::clamp(ineq_lower - ineq_vals, 0.0, arma::datum::inf);
        arma::vec upper_viol = arma::clamp(ineq_vals - ineq_upper, 0.0, arma::datum::inf);
        all_viol = arma::max(lower_viol, upper_viol);
        ineq_violation = all_viol.max(); // or arma::norm(all_viol, "inf")
    }

    if (n_eq > 0)
        eq_violation = arma::norm(Rcpp::as<arma::vec>(eq_f(parameters)), "inf");

    double dual_violation = 0.0;
    if (n_ineq > 0)
        dual_violation = arma::norm(arma::min(lagr_mult_ineq, 0.0), "inf");

    double comp_slack = 0.0;
    if (n_ineq > 0 && all_viol.n_elem == lagr_mult_ineq.n_elem)
        comp_slack = arma::norm(lagr_mult_ineq % all_viol, "inf");

    return Rcpp::List::create(
        Rcpp::_["kkt_stationarity"] = kkt_stationarity,
        Rcpp::_["eq_violation"] = eq_violation,
        Rcpp::_["ineq_violation"] = ineq_violation,
        Rcpp::_["dual_feas_violation"] = dual_violation,
        Rcpp::_["compl_slackness"] = comp_slack
    );
}


inline double compute_stationarity(
        const arma::vec& parameters,
        const arma::vec& lagrange_mults,
        int n_eq,
        int n_ineq,
        const Rcpp::Function& gradient_fun,
        const Rcpp::Function& eq_j,
        const Rcpp::Function& ineq_j
) {
    arma::vec grad = Rcpp::as<arma::vec>(gradient_fun(parameters));
    arma::vec lagr_mult_eq, lagr_mult_ineq;
    if (n_eq > 0) lagr_mult_eq = lagrange_mults.subvec(0, n_eq - 1);
    if (n_ineq > 0) lagr_mult_ineq = lagrange_mults.subvec(n_eq, n_eq + n_ineq - 1);
    arma::mat A_eq, A_ineq;

    if (n_eq > 0)   A_eq   = Rcpp::as<arma::mat>(eq_j(parameters));
    if (n_ineq > 0) A_ineq = Rcpp::as<arma::mat>(ineq_j(parameters));
    arma::vec lagr_grad = grad;
    if (n_eq > 0)   lagr_grad -= A_eq.t() * lagr_mult_eq;
    if (n_ineq > 0) lagr_grad -= A_ineq.t() * lagr_mult_ineq;

    double kkt_stationarity = arma::norm(lagr_grad, "inf");
    return kkt_stationarity;
}

#endif // SUBNP_HELPERS_H
