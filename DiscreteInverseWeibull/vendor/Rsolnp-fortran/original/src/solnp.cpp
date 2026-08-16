// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "subnp_state.h"
#include "subnp_helpers.h"

subnp_state::subnp_state(Rcpp::List& state)
    : pars(state["augmented_parameters"]),
      lagrange_mults(state["lagrange_mults"]),
      scaled_eval(state["scaled_eval"]),
      augmented_hessian(state["augmented_hessian"]),
      lower(state["lower_tmp"]),
      upper(state["upper_tmp"]),
      problem_indicators(state["problem_indicators"]),
      lambda(Rcpp::as<double>(state["lambda"])),
      penalty_param(state["penalty_param"]),
      tol(state["tol"]),
      ftol(state["ftol"]),
      min_iter(state["min_iter"]),
      trace(state["trace"]),
      solnp_fun(state["solnp_fun"]),
      solnp_gradfun(state["solnp_gradfun"]),
      solnp_eqfun(state["solnp_eqfun"]),
      solnp_ineqfun(state["solnp_ineqfun"]),
      solnp_eqjac(state["solnp_eqjac"]),
      solnp_ineqjac(state["solnp_ineqjac"])
{
    // Defensive assignment for possibly missing or NULL ineq_lower and ineq_upper
    if (state.containsElementNamed("ineq_lower") && !Rf_isNull(state["ineq_lower"])) {
        ineq_lower = Rcpp::as<arma::vec>(state["ineq_lower"]);
    } else {
        ineq_lower = arma::vec();
    }

    if (state.containsElementNamed("ineq_upper") && !Rf_isNull(state["ineq_upper"])) {
        ineq_upper = Rcpp::as<arma::vec>(state["ineq_upper"]);
    } else {
        ineq_upper = arma::vec();
    }

    // Defensive assignment for possibly missing or NULL scaling_factors
    if (state.containsElementNamed("scaling_factors") && !Rf_isNull(state["scaling_factors"])) {
        scaling_factors = Rcpp::as<arma::vec>(state["scaling_factors"]);
    } else {
        scaling_factors = arma::vec();
    }
    nfeval = 0.0;
}

// [[Rcpp::export(.csolnp)]]
Rcpp::List csolnp(Rcpp::List state) {
    int major_iteration_count = 0;
    int max_major_iterations = Rcpp::as<int>(state["max_major_iterations"]);
    double penalty_param = Rcpp::as<double>(state["penalty_param"]);
    double tol = Rcpp::as<double>(state["tol"]);
    //double ftol = Rcpp::as<double>(state["ftol"]);
    int n_eq = Rcpp::as<int>(state["n_eq"]);
    int n_ineq = Rcpp::as<int>(state["n_ineq"]);
    int num_parameters = Rcpp::as<int>(state["num_parameters"]);
    int total_constraints = Rcpp::as<int>(state["total_constraints"]);
    Rcpp::LogicalVector problem_indicators = state["problem_indicators"];
    arma::vec scaled_eval = Rcpp::as<arma::vec>(state["scaled_eval"]);
    arma::vec augmented_parameters = Rcpp::as<arma::vec>(state["augmented_parameters"]);
    arma::vec lagrange_mults = Rcpp::as<arma::vec>(state["lagrange_mults"]);
    arma::mat augmented_hessian = Rcpp::as<arma::mat>(state["augmented_hessian"]);
    double lambda = Rcpp::as<double>(state["lambda"]);
    arma::vec lower_tmp = Rcpp::as<arma::vec>(state["lower_tmp"]);
    arma::vec upper_tmp = Rcpp::as<arma::vec>(state["upper_tmp"]);
    arma::vec ineq_lower = (state.containsElementNamed("ineq_lower") && !Rf_isNull(state["ineq_lower"])) ? Rcpp::as<arma::vec>(state["ineq_lower"]) : arma::vec();
    arma::vec ineq_upper = (state.containsElementNamed("ineq_upper") && !Rf_isNull(state["ineq_upper"])) ? Rcpp::as<arma::vec>(state["ineq_upper"]) : arma::vec();
    arma::mat all_bounds = Rcpp::as<arma::mat>(state["all_bounds"]);
    Rcpp::Function objective_fun = state["solnp_fun"];
    Rcpp::Function gradient_fun = state["solnp_gradfun"];
    Rcpp::Function eq_f = state["solnp_eqfun"];
    Rcpp::Function ineq_f = state["solnp_ineqfun"];
    Rcpp::Function eq_j = state["solnp_eqjac"];
    Rcpp::Function ineq_j = state["solnp_ineqjac"];
    int trace = state["trace"];
    double previous_objective_value = Rcpp::as<double>(state["previous_objective_value"]);
    arma::vec status_vector = Rcpp::as<arma::vec>(state["status_vector"]);
    arma::vec historical_objective_values = Rcpp::as<arma::vec>(state["historical_objective_values"]);
    subnp_state subnp(state);
    subnp.lower = lower_tmp;
    subnp.upper = upper_tmp;
    subnp.ineq_lower = ineq_lower;
    subnp.ineq_upper = ineq_upper;
    subnp.penalty_param = penalty_param;
    subnp.solnp_fun = objective_fun;
    subnp.solnp_gradfun = gradient_fun;
    subnp.solnp_eqfun = eq_f;
    subnp.solnp_ineqfun = ineq_f;
    subnp.solnp_eqjac = eq_j;
    subnp.solnp_ineqjac = ineq_j;
    subnp.nfeval = 0;
    int n_fun_eval = 0;
    int error_code = 0;
    int convergence = 0;
    double gtol = 1e-6;
    const double min_penalty = std::sqrt(tol);
    double grad_norm = 0.0;
    arma::vec previous_parameters = augmented_parameters.subvec(n_ineq, n_ineq + num_parameters - 1);
    double best_feas_obj = std::numeric_limits<double>::infinity();
    //arma::vec best_feas_params;
    //double best_feas_constr = std::numeric_limits<double>::infinity();
    double step_norm = 0.0;
    double current_objective_value = std::numeric_limits<double>::infinity();
    arma::vec tmp_obj = arma::zeros(1);
    arma::vec best_lagrange_mults;
    while (major_iteration_count < max_major_iterations) {
        // Build/update control list per iteration
        subnp.penalty_param = penalty_param;
        major_iteration_count++;
        // 1. Objective/equality scaling
        arma::vec objective_and_eq_scale;
        if (problem_indicators[6] > 0) {
            double eq_max = arma::abs(scaled_eval.subvec(1, n_eq)).max(); // indices 1...n_eq (inclusive)
            objective_and_eq_scale = arma::vec(n_eq + 1, arma::fill::ones) * eq_max;
            objective_and_eq_scale(0) = scaled_eval(0);
        } else {
            objective_and_eq_scale = arma::vec(1);
            objective_and_eq_scale(0) = 1.0;
        }

        // 2. Scaling factors
        arma::vec scaling_factors;
        if (problem_indicators[10] == 0) {
            scaling_factors = arma::join_vert(objective_and_eq_scale, augmented_parameters);
        } else {
            scaling_factors = arma::join_vert(objective_and_eq_scale, arma::ones(augmented_parameters.n_elem));
        }
        // Clamp scaling factors
        for (auto& val : scaling_factors) {
            val = std::min(std::max(std::abs(val), tol), 1.0 / tol);
        }

        // 3. Update subnp_state with fields that change per iteration
        subnp.pars = augmented_parameters;
        subnp.lagrange_mults = lagrange_mults;
        subnp.scaled_eval = scaled_eval;
        subnp.augmented_hessian = augmented_hessian;
        subnp.lambda = lambda;
        subnp.scaling_factors = scaling_factors;

        // 4. Call inner solver
        Rcpp::List subnp_results = csubnp_cpp(subnp);
        error_code = Rcpp::as<int>(subnp_results["solnp_error"]);
        n_fun_eval += Rcpp::as<int>(subnp_results["nfeval"]);
        // 5. Check for error from subproblem
        if (error_code == 1) {
            max_major_iterations = major_iteration_count;
        }

        // 6. Update parameter state
        augmented_parameters = Rcpp::as<arma::vec>(subnp_results["p"]);

        lagrange_mults = Rcpp::as<arma::vec>(subnp_results["y"]);
        augmented_hessian = Rcpp::as<arma::mat>(subnp_results["augmented_hessian"]);
        lambda = Rcpp::as<double>(subnp_results["lambda"]);

        // 7. Extract current parameters from augmented vector
        arma::vec current_parameters = augmented_parameters.subvec(n_ineq, n_ineq + num_parameters - 1);

        step_norm = arma::norm(current_parameters - previous_parameters, 2);
        previous_parameters = current_parameters;

        // 8. Evaluate objective
        current_objective_value = Rcpp::as<double>(objective_fun(current_parameters));
        n_fun_eval += 1;
        if (trace > 0) {
            print_progress(major_iteration_count, current_objective_value, status_vector(2),
                           status_vector(0), step_norm, penalty_param);
        }
        // this requires to also add grad_norm to be robust
        // if (status_vector(2) < ftol) {
        //   if (current_objective_value < best_feas_obj) {
        //     Rcpp::Rcout<<"current_objective_value"<<current_objective_value<<std::endl;
        //     best_feas_obj = current_objective_value;
        //     best_feas_params = augmented_parameters;
        //     best_feas_constr = status_vector(2);
        //     best_lagrange_mults = lagrange_mults;
        //   }
        // }

        // 10. Evaluate constraints
        arma::vec combined(1);
        combined(0) = current_objective_value;
        if (n_eq > 0) {
          arma::vec current_eq_values = Rcpp::as<arma::vec>(eq_f(current_parameters));
          combined = arma::join_vert(combined, current_eq_values);
        }
        if (n_ineq > 0) {
          arma::vec current_ineq_values = Rcpp::as<arma::vec>(ineq_f(current_parameters));
          combined = arma::join_vert(combined, current_ineq_values);
        }
        scaled_eval = combined;

        // 11. Objective relative change
        status_vector(0) = (previous_objective_value - scaled_eval(0)) /
            std::max(std::abs(scaled_eval(0)), 1.0);
        previous_objective_value = scaled_eval(0);

        // 12. Constraint violation and penalty parameter update
        if (total_constraints > 0) {
            arma::vec current_constraint_violations = scaled_eval.subvec(1, total_constraints);

            if (problem_indicators[3] > 0) {
                arma::vec temp_ineq_slack_lb = current_constraint_violations.subvec(n_eq, total_constraints - 1) - all_bounds.col(0).rows(0, n_ineq - 1);
                arma::vec temp_ineq_slack_ub = all_bounds.col(1).rows(0, n_ineq - 1) - current_constraint_violations.subvec(n_eq, total_constraints - 1);
                arma::vec all_slack = arma::join_vert(temp_ineq_slack_lb, temp_ineq_slack_ub);
                if (all_slack.min() > 0) {
                    augmented_parameters.subvec(0, n_ineq - 1) = current_constraint_violations.subvec(n_eq, total_constraints - 1);
                }
                current_constraint_violations.subvec(n_eq, total_constraints - 1) -= augmented_parameters.subvec(0, n_ineq - 1);
            }
            // Norm of constraint violations
            status_vector(2) = vnorm(current_constraint_violations);
            if (error_code == 0) {
              grad_norm = compute_stationarity(current_parameters, lagrange_mults, n_eq, n_ineq, gradient_fun, eq_j, ineq_j);
              if (status_vector(2) < tol && grad_norm < gtol) {
                penalty_param = 0;
                lambda = std::min(lambda, tol);
              } else {
                if (status_vector(2) < 5 * status_vector(1)) {
                  penalty_param = std::max(penalty_param / 5.0, min_penalty);
                }
                if (status_vector(2) > 10 * status_vector(1)) {
                  penalty_param = 5.0 * std::max(penalty_param, min_penalty);
                }
              }
            } else {
              if (status_vector(2) < tol) {
                penalty_param = 0;
                lambda = std::min(lambda, tol);
              } else {
                if (status_vector(2) < 5 * status_vector(1)) {
                  penalty_param = std::max(penalty_param / 5.0, min_penalty);
                }
                if (status_vector(2) > 10 * status_vector(1)) {
                  penalty_param = 5.0 * std::max(penalty_param, min_penalty);
                }
              }
            }
            if (penalty_param == 0 && status_vector(2) > tol) {
              penalty_param = std::max(1.0, min_penalty);
              if (trace > 0) {
                Rcpp::Rcout << "Penalty restored to " << penalty_param << " due to loss of feasibility." << std::endl;
              }
            }
            // Reset multipliers/Hessian if needed
            if (std::max(tol + status_vector(0), status_vector(1) - status_vector(2)) <= 0) {
                lagrange_mults.zeros();
                augmented_hessian = arma::diagmat(augmented_hessian.diag());
            }
            status_vector(1) = status_vector(2);
        }

        // 13. Convergence check
        if (vnorm(arma::vec({status_vector(0), status_vector(1)})) <= tol) {
          max_major_iterations = major_iteration_count;
          convergence = 0;
        } else {
          convergence = 1;
        }
        if (error_code == 1) {
          convergence = 2;
        }

        // 14. Track objective
        tmp_obj(0) = current_objective_value;
        historical_objective_values = arma::join_vert(historical_objective_values, tmp_obj);
    }
    best_feas_obj = current_objective_value;

    // if (best_feas_obj < std::numeric_limits<double>::infinity() && std::abs(current_objective_value - best_feas_obj) > 1e-10) {
    //   augmented_parameters = best_feas_params;
    //   lagrange_mults = best_lagrange_mults;
    //   status_vector(2) = best_feas_constr;
    //   Rcpp::warning("Final iterate is infeasible. Returning best feasible solution found.");
    // } else {
    //   best_feas_obj = current_objective_value;
    // }
    arma::vec optimal_parameters = augmented_parameters.subvec(n_ineq, n_ineq + num_parameters - 1);
    Rcpp::List kkt_diagnostics = compute_kkt_diagnostics(optimal_parameters, lagrange_mults, n_eq, n_ineq,
                                                         gradient_fun, eq_j, ineq_j, eq_f, ineq_f, ineq_lower, ineq_upper, tol, error_code);
    // --- Return results ---
    return Rcpp::List::create(
        Rcpp::_["parameters"] = augmented_parameters,
        Rcpp::_["lagrange_mults"] = lagrange_mults,
        Rcpp::_["augmented_hessian"] = augmented_hessian,
        Rcpp::_["lambda"] = lambda,
        Rcpp::_["status_vector"] = status_vector,
        Rcpp::_["major_iteration_count"] = major_iteration_count,
        Rcpp::_["n_fun_evaluations"] = n_fun_eval,
        Rcpp::_["error_code"] = error_code,
        Rcpp::_["convergence"] = convergence,
        Rcpp::_["historical_objective_values"] = historical_objective_values,
        Rcpp::_["best_objective"] = best_feas_obj,
        Rcpp::_["kkt_diagnostics"] = kkt_diagnostics
    );
}

