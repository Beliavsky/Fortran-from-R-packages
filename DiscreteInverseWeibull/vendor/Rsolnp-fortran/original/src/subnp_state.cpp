#include "subnp_state.h"
#include "subnp_helpers.h"
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

Rcpp::List csubnp_cpp(subnp_state& s)
{
    Rcpp::IntegerVector setup = s.problem_indicators;
    arma::uword n_eq = setup[7];     // Be careful: R is 1-based, C++ is 0-based!
    arma::uword n_ineq = setup[4];
    arma::uword n_pars = setup[0];
    arma::uword n_constraints = n_eq + n_ineq;
    arma::uword n_pic = n_pars + n_ineq;
    arma::uword mm = 0;              // used later for bounds
    double rho = s.penalty_param;
    int maxit = s.min_iter;
    double tol = s.tol;
    int trace = s.trace;
    bool success;
    int minit = 0;
    int nfeval = s.nfeval;
    double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
    arma::vec working_params = s.pars;
    arma::vec line_search_steps = arma::zeros(3);
    arma::vec step_objectives = arma::zeros(3);
    arma::vec lagrange_mults = s.lagrange_mults;
    arma::vec scaled_value = s.scaled_eval;
    arma::vec obm = arma::zeros(1 + n_eq + n_ineq);
    arma::mat augmented_hessian = s.augmented_hessian;
    double lambda = s.lambda;

    arma::vec r, u, dx, y;
    arma::vec JTy, JTr;
    arma::mat J;
    // initialize some variables
    double fun_value = 0.0;
    double j = 0.0;
    double reduce = std::numeric_limits<double>::infinity();
    arma::vec grad_f;
    arma::vec grad_s;
    arma::vec scaled_grad_f;
    arma::vec scaled_grad_x;
    arma::vec x;
    arma::vec eq_value;
    arma::vec ineq_value;
    arma::mat Aineq;
    arma::mat Aeq;
    arma::mat Alin_x;
    arma::mat param_trials;
    arma::mat augmented_jacobian;
    arma::vec bfgs_scaling_factors = arma::zeros(2);
    int status_flag = 1;
    double step_interval_width = 1.0;
    double max_step_obj, min_step_obj;
    bool condition_1, condition_2, condition_3;
    arma::vec lower = s.lower;
    arma::vec upper = s.upper;
    arma::vec ineq_lower, ineq_upper;
    arma::vec scaling_factors = s.scaling_factors;

    if (n_ineq > 0) {
        ineq_lower = s.ineq_lower;
        ineq_upper = s.ineq_upper;
    }
    scaled_value /= scaling_factors.subvec(0, n_constraints);
    working_params /= scaling_factors.subvec(n_eq + 1, n_constraints + n_pars);

    arma::mat param_bounds(n_pic, 2, arma::fill::zeros);
    if (n_ineq > 0) {
        param_bounds.submat(0, 0, n_ineq-1, 0) = ineq_lower;
        param_bounds.submat(0, 1, n_ineq-1, 1) = ineq_upper;
    }
    if (n_pars > 0) {
        param_bounds.submat(n_ineq, 0, n_pic-1, 0) = lower;
        param_bounds.submat(n_ineq, 1, n_pic-1, 1) = upper;
    }

    if (setup(10) > 0) {
        mm = (setup(9) == 0 ? n_ineq : n_pic);
        arma::vec scale_vec = scaling_factors.subvec(n_eq + 1, n_eq + mm);
        param_bounds.cols(0, 1).each_col() /= scale_vec;
    }


    auto solnp_error_return = [&](const arma::vec& current_params, arma::vec& y_local, arma::mat& hess_local,
                                  double lambda_local, const arma::vec& scaling_factors_local,
                                  arma::uword n_eq_local, arma::uword n_constraints_local, arma::uword n_pars_local,
                                  bool constraints_exist, const std::string& error_message, int nfeval, int error_code = 1) {
        arma::vec scaling_tmp = scaling_factors_local.subvec(n_eq_local + 1, n_constraints_local + n_pars_local);
        arma::vec p_out = current_params % scaling_tmp;
        if (constraints_exist) {
            y_local = y_local.zeros();
        }
        hess_local = scaling_factors_local(0) * hess_local / (scaling_tmp * scaling_tmp.t());
        return Rcpp::List::create(
            Rcpp::_["p"] = p_out,
            Rcpp::_["y"] = y_local,
            Rcpp::_["augmented_hessian"] = hess_local,
            Rcpp::_["lambda"] = lambda_local,
            Rcpp::_["solnp_error"] = error_code,
            Rcpp::_["message"] = error_message,
            Rcpp::_["nfeval"] = nfeval
        );
    };

    if (n_constraints > 0) {
        lagrange_mults = scaling_factors.subvec(1, n_constraints) % lagrange_mults / scaling_factors[0];
    }
    augmented_hessian = augmented_hessian % (scaling_factors.subvec(n_eq + 1, n_constraints + n_pars) * scaling_factors.subvec(n_eq + 1, n_constraints + n_pars).t()) / scaling_factors(0);

    j = scaled_value(0);

    if (setup[6] > 0 && setup[3] > 0) {
        arma::mat Aeq_block = arma::zeros(n_eq, n_ineq + n_pars);
        arma::mat Aineq_block = join_rows(-arma::eye(n_ineq, n_ineq), arma::zeros(n_ineq, n_pars));
        augmented_jacobian = join_cols(Aeq_block, Aineq_block);
    } else if (setup[3] > 0) {
        augmented_jacobian = join_rows(-arma::eye(n_ineq, n_ineq), arma::zeros(n_ineq, n_pars));
    } else if (setup[6] > 0) {
        augmented_jacobian = arma::zeros(n_eq, n_pars);
    }

    arma::vec grad_augmented_obj = arma::zeros(n_pic);
    arma::vec p = working_params.subvec(0,n_pic - 1);
    arma::vec constraint, constraint_offset;

    if (n_constraints > 0) {
        constraint = scaled_value.subvec(1, n_constraints);
        x = working_params.subvec(n_ineq, n_pic-1) % scaling_factors.subvec(n_constraints+1, n_constraints+n_pars);
        if (setup[6] > 0) {
            Aeq = Rcpp::as<mat>(s.solnp_eqjac(x));
            arma::vec r = 1.0 / scaling_factors.subvec(1, n_eq); // [n_eq]
            arma::vec c = scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars); // [n_pars]
            Aeq = Aeq % (r * c.t());
            augmented_jacobian.submat(0, n_ineq, n_eq-1, n_pic-1) = Aeq;
        }
        if (setup[3] > 0) {
            Aineq = Rcpp::as<mat>(s.solnp_ineqjac(x));
            arma::vec r = 1.0 / scaling_factors.subvec(n_eq + 1, n_constraints);              // length: n_ineq
            arma::vec c = scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars);  // length: n_pars
            Aineq = Aineq % (r * c.t());
            augmented_jacobian.submat(n_eq, n_ineq, n_constraints-1, n_pic-1) = Aineq;
            augmented_jacobian.submat(n_eq, 0, n_constraints-1, n_ineq-1) = -arma::eye(n_ineq, n_ineq);
        }

        arma::vec g_user = Rcpp::as<vec>(s.solnp_gradfun(x));
        g_user /= (scaling_factors.subvec(n_constraints+1, n_constraints+n_pars) * scaling_factors(0));
        grad_augmented_obj.subvec(n_ineq, n_pic-1) = g_user;
        if (setup[3] > 0) {
            constraint.subvec(n_eq, n_eq + n_ineq - 1) -= working_params.subvec(0, n_ineq-1);
        }
        constraint_offset = augmented_jacobian * working_params - constraint;

        status_flag = -1;
        line_search_steps(0) = tol - arma::abs(constraint).max();

        if (line_search_steps(0) <= 0) {
            status_flag = 1;
            if (setup[10] == 0) {
                working_params -= augmented_jacobian.t() * arma::solve(augmented_jacobian * augmented_jacobian.t(), constraint, arma::solve_opts::fast + arma::solve_opts::no_approx);
                line_search_steps(1) = 1;
            }
        }
        if (line_search_steps(0) <= 0) {
            if (working_params.n_elem <= n_pic) {
                // Grow the vector to (n_pic + 1), filling new elements with 0
                arma::vec tmp = arma::zeros<arma::vec>(n_pic + 1);
                tmp.subvec(0, working_params.n_elem - 1) = working_params;
                working_params = tmp;
            }
            working_params(n_pic) = 1.0;
            augmented_jacobian = arma::join_rows(augmented_jacobian, -constraint);
            arma::rowvec cx = arma::join_rows(arma::rowvec(n_pic, arma::fill::zeros), arma::rowvec({1.0}));
            dx = arma::ones<arma::vec>(n_pic + 1);
            step_interval_width = 1.0;
            minit = 0;
            while(step_interval_width >= tol) {
                minit += 1;
                arma::mat gap = arma::join_rows(working_params.subvec(0, mm-1) - param_bounds.col(0), param_bounds.col(1) - working_params.subvec(0, mm-1));
                for (arma::uword k = 0; k < gap.n_rows; ++k) {
                    gap.row(k) = arma::sort(gap.row(k));
                }
                dx.subvec(0, mm - 1) = gap.col(0);
                dx(n_pic)  = working_params(n_pic);
                if (setup[9] == 0) {
                    dx.subvec(mm, n_pic - 1) = std::max(dx.subvec(0, mm-1).max(), 100.0) * arma::ones<arma::vec>(n_pic - mm);
                }
                auto result =  qr_solve_weighted_system(augmented_jacobian, dx, cx);
                success = result.second;
                if (!success) {
                    return solnp_error_return(working_params, y, augmented_hessian, lambda, scaling_factors, n_eq, n_constraints, n_pars, n_constraints > 0, "QR decomposition failed", nfeval, 1);
                }
                y = result.first;
                arma::vec v = dx % (dx % (cx.t() - augmented_jacobian.t() * y));
                if (v(n_pic) > 0) {
                    double z = working_params(n_pic) / v(n_pic);
                    for (arma::uword i = 0; i < mm; ++i) {
                        if (v(i) < 0) {
                            z = std::min(z, -(param_bounds(i, 1) - working_params(i)) / v(i)); // R: param_bounds[i, 2]
                        } else if (v(i) > 0) {
                            z = std::min(z, (working_params(i) - param_bounds(i, 0)) / v(i)); // R: param_bounds[i, 1]
                        }
                    }
                    if (z >= working_params(n_pic) / v(n_pic)) {
                        working_params -= z * v;
                    } else {
                        working_params -= 0.9 * z * v;
                    }
                    step_interval_width = working_params(n_pic);
                    if (minit >= 10) {
                        step_interval_width = 0;
                    }
                } else {
                    step_interval_width = 0;
                    minit = 10;
                }
            }

            if (minit >= 10) {
                if (trace > 1) solver_warnings("M2");
            }
            augmented_jacobian = augmented_jacobian.cols(0, n_pic - 1);
            constraint_offset = augmented_jacobian * working_params.subvec(0, n_pic - 1);
        }
    }

    p = working_params.subvec(0, n_pic - 1);
    y = arma::zeros<arma::vec>(1);
    if (status_flag > 0) {
        x = p.subvec(n_ineq, n_pic - 1) % scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars);
        fun_value = as<double>(s.solnp_fun(x));
        nfeval += 1;
        eq_value = (n_eq > 0) ? as<arma::vec>(s.solnp_eqfun(x)) : arma::zeros(0);
        ineq_value = (n_ineq > 0) ? as<arma::vec>(s.solnp_ineqfun(x)) : arma::zeros(0);
        scaled_value(0) = fun_value;
        if (n_eq > 0) scaled_value.subvec(1, n_eq) = eq_value;
        if (n_ineq > 0) scaled_value.subvec(n_eq + 1, n_constraints) = ineq_value;
        scaled_value = scaled_value / scaling_factors.subvec(0, n_constraints);
    }

    j = scaled_value(0);

    if (setup[3] > 0) {
        scaled_value.subvec(n_eq + 1,n_constraints) -=  p.subvec(0, n_ineq - 1);
    }

    if (n_constraints > 0) {
        scaled_value.subvec(1, n_constraints) = scaled_value.subvec(1, n_constraints) - augmented_jacobian * p + constraint_offset;
        arma::vec ceval = scaled_value.subvec(1, n_constraints);
        j = augmented_lagrangian(scaled_value, lagrange_mults, rho, n_constraints);
    }


    minit = 0;
    arma::vec delta_position = arma::zeros(n_pic);
    arma::vec delta_gradient = arma::zeros(n_pic);
    arma::vec scaled_value_1, scaled_value_2, scaled_value_3;

    while (minit < maxit) {
        minit++;
        if (status_flag > 0) {
            x = p.subvec(n_ineq,n_pic - 1) % scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars);
            fun_value = as<double>(s.solnp_fun(x));
            nfeval += 1;
            eq_value = (n_eq > 0) ? as<arma::vec>(s.solnp_eqfun(x)) : arma::zeros(0);
            ineq_value = (n_ineq > 0) ? as<arma::vec>(s.solnp_ineqfun(x)) : arma::zeros(0);
            grad_f = as<arma::vec>(s.solnp_gradfun(x));
            scaled_grad_f = grad_f % scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars);
            obm(0) = fun_value;
            if (n_eq > 0) obm.subvec(1, n_eq) = eq_value;
            if (n_ineq > 0) obm.subvec(n_eq + 1, n_constraints) = ineq_value;
            obm = obm / scaling_factors.subvec(0, n_constraints);
            scaled_grad_x = scaled_grad_f / scaling_factors(0);
            if (n_eq > 0) {
                Aeq = Rcpp::as<mat>(s.solnp_eqjac(x));
                arma::vec r = 1.0 / scaling_factors.subvec(1, n_eq); // [n_eq]
                arma::vec c = scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars); // [n_pars]
                Aeq = Aeq % (r * c.t());
            }
            if (n_ineq > 0) {
                Aineq = Rcpp::as<mat>(s.solnp_ineqjac(x));
                arma::vec r = 1.0 / scaling_factors.subvec(n_eq + 1, n_constraints);              // length: n_ineq
                arma::vec c = scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars);  // length: n_pars
                Aineq = Aineq % (r * c.t());
            }
            if (n_constraints > 0) {
                Alin_x = augmented_jacobian.cols(n_ineq, n_ineq + n_pars - 1);
                if (n_eq > 0 && n_ineq > 0) {
                    J = arma::join_vert(Aeq, Aineq) - Alin_x;
                } else if (n_eq > 0) {
                    J = Aeq - Alin_x;
                } else if (n_ineq > 0) {
                    J = Aineq - Alin_x;
                }
                r = obm.subvec(1, n_constraints);
                if (n_ineq > 0) {
                    r.subvec(n_eq, n_constraints - 1) -= p.subvec(0, n_ineq  - 1);
                    r -= augmented_jacobian * p;
                    r += constraint_offset;
                }
                JTy = J.t() * lagrange_mults;
                JTr = J.t() * r;
                scaled_grad_x = scaled_grad_x - JTy + 2 * rho * JTr;
                if (n_ineq > 0) {
                    grad_s = lagrange_mults.subvec(n_eq, n_constraints - 1) - 2 * rho * r.subvec(n_eq, n_constraints - 1);
                }
            }
            if (n_ineq > 0) grad_augmented_obj.subvec(0, n_ineq - 1) = grad_s;
            grad_augmented_obj.subvec(n_ineq, n_ineq + n_pars - 1) = scaled_grad_x;

        }
        if (setup[3] > 0) {
            grad_augmented_obj.subvec(0, n_ineq - 1).fill(0.0);
        }

        if (minit > 1) {
            delta_gradient = grad_augmented_obj - delta_gradient;
            delta_position = p - delta_position;
            bfgs_scaling_factors(0) = arma::as_scalar(delta_position.t() * augmented_hessian * delta_position);
            bfgs_scaling_factors(1) = arma::as_scalar(delta_position.t() * delta_gradient);
            if ((bfgs_scaling_factors(0) * bfgs_scaling_factors(1)) > 0) {
                delta_position = augmented_hessian * delta_position;
                augmented_hessian = augmented_hessian - (delta_position * delta_position.t()) / bfgs_scaling_factors(0)
                    + (delta_gradient * delta_gradient.t()) / bfgs_scaling_factors(1);
            }

        }
        dx = 0.01 * arma::ones<arma::vec>(n_pic);
        if (setup[10] > 0) {
            // Compute gap: [p[1:mm] - param_bounds[,1], param_bounds[,2] - p[1:mm]]
            arma::mat gap = arma::join_rows(
                p.subvec(0, mm-1) - param_bounds.col(0),
                param_bounds.col(1) - p.subvec(0, mm-1)
            );
            // Row-wise sort (for each row, sort the two elements)
            for (arma::uword i = 0; i < gap.n_rows; ++i)
                gap.row(i) = arma::sort(gap.row(i));

            // Take the smallest element of each row, add sqrt(eps)
            arma::vec gap_vec = gap.col(0) + epsilon * arma::ones<arma::vec>(mm);

            // Update first mm elements of dx
            dx.subvec(0, mm-1) = arma::ones<arma::vec>(mm) / gap_vec;

            if (setup[9] == 0) {
                double min_dx = dx.subvec(0, mm-1).min();
                double val = std::min(min_dx, 0.01); // matching R's min(c(dx[1:mm, 1], 0.01))
                dx.subvec(mm, n_pic-1) = val * arma::ones<arma::vec>(n_pic - mm);
            }
        }

        step_interval_width = -1;
        lambda /= 10.0;
        // ToDo: add these as control options
        int step_count = 0.0;
        const double step_eps = 1e-12;
        double previous_step_interval_width = step_interval_width;
        while (step_interval_width <= 0) {
            arma::mat diag_dx2 = arma::diagmat(arma::square(dx));
            arma::mat hess_reg = augmented_hessian + lambda * diag_dx2;
            arma::mat cz_temp;
            arma::mat cz;
            auto result = cholesky(hess_reg);
            bool chol_success = result.second;
            if (!chol_success) {
                return solnp_error_return(p, y, augmented_hessian, lambda, scaling_factors, n_eq, n_constraints, n_pars, n_constraints > 0, "Cholesky decomposition failed", nfeval, 1);
            }

            cz_temp = result.first;
            bool inv_success = false;
            try {
                cz = arma::inv(cz_temp, arma::inv_opts::allow_approx);
                inv_success = true;
            } catch (std::runtime_error& e) {
                cz_temp.reset();
                inv_success = false;
            }
            if (!inv_success) {
                return solnp_error_return(p, y, augmented_hessian, lambda, scaling_factors, n_eq, n_constraints, n_pars, n_constraints > 0, "Cholesky Matrix inversion failed", nfeval, 1);
            }
            delta_gradient = cz.t() * grad_augmented_obj;
            // || augmented_jacobian.n_rows == 0

            if (n_constraints == 0) {
                u = -cz * delta_gradient;
            } else {
                auto result = qr_solve_transposed_system(cz, augmented_jacobian, delta_gradient);
                success = result.second;
                if (!success) {
                    return solnp_error_return(p, y, augmented_hessian, lambda, scaling_factors, n_eq, n_constraints, n_pars, n_constraints > 0, "QR decomposition failed", nfeval, 1);
                }
                y = result.first;
                u = -cz * (delta_gradient - (cz.t() * augmented_jacobian.t()) * y);
            }

            working_params = u.subvec(0, n_pic-1) + p;

            if (setup[10] == 0) {  // R's setup[12], C++ is zero-based
                step_interval_width = 1;
            } else {
                arma::vec diff1 = working_params.subvec(0, mm-1) - param_bounds.col(0);
                arma::vec diff2 = param_bounds.col(1) - working_params.subvec(0, mm-1);
                step_interval_width = std::min(diff1.min(), diff2.min());
                lambda *= 3;
            }
            if (std::abs(step_interval_width - previous_step_interval_width) < step_eps) {
                step_count += 1;
            } else {
                previous_step_interval_width = step_interval_width;
            }
            if (step_count > 25) {
                return solnp_error_return(p, y, augmented_hessian, lambda, scaling_factors, n_eq, n_constraints, n_pars, n_constraints > 0, "Infeasible Line Search", nfeval, 1);
            }
        }

        line_search_steps(0) = 0;
        scaled_value_1 = scaled_value;
        scaled_value_2 = scaled_value_1;
        scaled_value_3 = 0.0 * scaled_value;
        step_objectives(0) = j;
        step_objectives(1) = j;
        param_trials = arma::join_rows(p, p);
        line_search_steps(2) = 1.0;
        param_trials = arma::join_rows(param_trials, working_params);
        x = param_trials.col(2).rows(n_ineq, n_pic - 1) % scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars);
        fun_value = as<double>(s.solnp_fun(x));
        nfeval += 1;
        eq_value = (n_eq > 0) ? as<arma::vec>(s.solnp_eqfun(x)) : arma::zeros(0);
        ineq_value = (n_ineq > 0) ? as<arma::vec>(s.solnp_ineqfun(x)) : arma::zeros(0);
        scaled_value_3(0) = fun_value;
        if (n_eq > 0) scaled_value_3.subvec(1, n_eq) = eq_value;
        if (n_ineq > 0) scaled_value_3.subvec(n_eq + 1, n_constraints) = ineq_value;
        scaled_value_3 = scaled_value_3 / scaling_factors.subvec(0, n_constraints);


        step_objectives(2) = scaled_value_3(0);

        if (setup[3] > 0) {
            scaled_value_3.subvec(n_eq + 1,n_constraints) -= param_trials.col(2).rows(0,n_ineq - 1);
        }
        if (n_constraints > 0) {
            scaled_value_3.subvec(1, n_constraints) = scaled_value_3.subvec(1, n_constraints) - augmented_jacobian * param_trials.col(2) + constraint_offset;
            step_objectives(2) = augmented_lagrangian(scaled_value_3, lagrange_mults, rho, n_constraints);
        }

        step_interval_width = 1;

        while (step_interval_width > tol) {
            line_search_steps(1) = (line_search_steps(0) + line_search_steps(2)) / 2.0;
            param_trials.col(1) = (1.0 - line_search_steps(1)) * p + line_search_steps(1) * working_params;
            x = param_trials.col(1).rows(n_ineq, n_pic - 1) % scaling_factors.subvec(n_constraints + 1, n_constraints + n_pars);
            fun_value = as<double>(s.solnp_fun(x));
            nfeval += 1;
            eq_value = (n_eq > 0) ? as<arma::vec>(s.solnp_eqfun(x)) : arma::zeros(0);
            ineq_value = (n_ineq > 0) ? as<arma::vec>(s.solnp_ineqfun(x)) : arma::zeros(0);
            scaled_value_2(0) = fun_value;
            if (n_eq > 0) scaled_value_2.subvec(1, n_eq) = eq_value;
            if (n_ineq > 0) scaled_value_2.subvec(n_eq + 1, n_constraints) = ineq_value;
            scaled_value_2 = scaled_value_2 / scaling_factors.subvec(0, n_constraints);
            step_objectives(1) = scaled_value_2(0);
            if (setup[3] > 0) {
                scaled_value_2.subvec(n_eq + 1,n_constraints) -= param_trials.col(1).rows(0,n_ineq - 1);
            }
            if (n_constraints > 0) {
                scaled_value_2.subvec(1, n_constraints) = scaled_value_2.subvec(1, n_constraints) - augmented_jacobian * param_trials.col(1) + constraint_offset;
                step_objectives(1) = augmented_lagrangian(scaled_value_2, lagrange_mults, rho, n_constraints);
            }
            max_step_obj = step_objectives.max();
            if (max_step_obj < j) {
                min_step_obj = step_objectives.min();
                step_interval_width = tol * (max_step_obj - min_step_obj) / (j - max_step_obj);
            }
            condition_1 = step_objectives(1) >= step_objectives(0);
            condition_2 = step_objectives(0) <= step_objectives(2) && step_objectives(1) < step_objectives(0);
            condition_3 = step_objectives(1) <  step_objectives(0) && step_objectives(0) > step_objectives(2);

            if (condition_1) {
                step_objectives(2) = step_objectives(1);
                scaled_value_3 = scaled_value_2;
                line_search_steps(2) = line_search_steps(1);
                param_trials.col(2) = param_trials.col(1);
            }

            if (condition_2) {
                step_objectives(2) = step_objectives(1);
                scaled_value_3 = scaled_value_2;
                line_search_steps(2) = line_search_steps(1);
                param_trials.col(2) = param_trials.col(1);
            }

            if (condition_3) {
                step_objectives(0) = step_objectives(1);
                scaled_value_1 = scaled_value_2;
                line_search_steps(0) = line_search_steps(1);
                param_trials.col(0) = param_trials.col(1);
            }

            if (step_interval_width >= tol) {
                step_interval_width = line_search_steps(2) - line_search_steps(0);
            }
        }
        delta_position = p;
        delta_gradient = grad_augmented_obj;
        status_flag = 1;

        min_step_obj = step_objectives.min();
        if (j <= min_step_obj) {
            maxit = minit;
        }
        reduce = (j - min_step_obj) / (1.0 + std::abs(j));
        if (reduce < tol) {
            maxit = minit; // break early
        }
        condition_1 = step_objectives(0) < step_objectives(1);
        condition_2 = step_objectives(2) < step_objectives(1) && step_objectives(0) >= step_objectives(1);
        condition_3 = step_objectives(0) >=  step_objectives(1) && step_objectives(2) >= step_objectives(1);
        if (condition_1) {
            j = step_objectives(0);
            p = param_trials.col(0);
            scaled_value = scaled_value_1;
        }
        if (condition_2) {
            j = step_objectives(2);
            p = param_trials.col(2);
            scaled_value = scaled_value_3;
        }
        if (condition_3) {
            j = step_objectives(1);
            p = param_trials.col(1);
            scaled_value = scaled_value_2;
        }

    }
    p = p % scaling_factors.subvec(n_eq + 1, n_constraints + n_pars);
    if (n_constraints > 0) {
        y = scaling_factors(0) * y / scaling_factors.subvec(1, n_constraints);
    }

    arma::mat scaling_outer = scaling_factors.subvec(n_eq + 1, n_constraints + n_pars) * scaling_factors.subvec(n_eq + 1, n_constraints + n_pars).t(); // outer product
    augmented_hessian = scaling_factors(0) * augmented_hessian / scaling_outer;

    if (reduce > tol) {
        if (trace > 1) solver_warnings("M3");
    }
    return Rcpp::List::create(
        Rcpp::_["p"] = p,
        Rcpp::_["y"] = y,
        Rcpp::_["augmented_hessian"] = augmented_hessian,
        Rcpp::_["lambda"] = lambda,
        Rcpp::_["solnp_error"] = 0,
        Rcpp::_["message"] = "Success",
        Rcpp::_["nfeval"] = nfeval);
}
