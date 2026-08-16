.rsolnp <- function(pars, fn, gr = NULL, eq_fn = NULL, eq_b = NULL, eq_jac = NULL,
         ineq_fn = NULL, ineq_lower = NULL, ineq_upper = NULL,
         ineq_jac = NULL, lower = NULL, upper = NULL, control = list(), ...)
{
    # Start timer for performance measurement
    start_time <- Sys.time()
    # Store original parameter names
    parameter_names <- names(pars)

    if (is.null(lower)) {
        lower <- -1000 * abs(pars)
        warning("\nlower values are NULL. Setting to -1000 x abs(pars).")
    }

    if (is.null(upper)) {
        upper <- 1000 * abs(pars)
        warning("\nupper values are NULL. Setting to 1000 x abs(pars).")
    }

    nudge_factor <- 1e-5

    # Clamp parameters at or below lower bound
    idx_lower <- which(pars <= lower)
    if (length(idx_lower) > 0) {
        tmp_p <- lower[idx_lower]
        # For strict zero bounds, nudge slightly up
        tmp_p[tmp_p == 0] <- 1e-12
        # Nudge up by a tiny fraction
        tmp_p <- tmp_p + abs(tmp_p) * nudge_factor
        pars[idx_lower] <- tmp_p
        warning("\nsome parameters at lower bounds. Adjusting by a factor of 1e-5.")
    }

    # Clamp parameters at or above upper bound
    idx_upper <- which(pars >= upper)
    if (length(idx_upper) > 0) {
        tmp_p <- upper[idx_upper]
        # For strict zero upper bounds (rare), nudge slightly down
        tmp_p[tmp_p == 0] <- -1e-12
        # Nudge down by a tiny fraction
        tmp_p <- tmp_p - abs(tmp_p) * nudge_factor
        pars[idx_upper] <- tmp_p
        warning("\nsome parameters at upper bounds. Adjusting by a factor of 1e-5.")
    }

    # Create a dedicated environment for csolnp's internal state
    .solnp_environment <- environment()
    assign("parameter_names", parameter_names, envir = .solnp_environment)
    # Initialize function call and error counters
    assign(".solnp_nfn", 0, envir = .solnp_environment)
    assign(".solnp_errors", 0, envir = .solnp_environment)
    # Indicator vector for problem characteristics
    # [1] Number of parameters (np)
    # [4] Has inequality constraints
    # [5] Number of inequality constraints
    # [7] Has equality constraints
    # [8] Number of equality constraints
    # [10] Has parameter upper/lower bounds
    # [11] Has any bounds (parameter or inequality)
    problem_indicators <- solnp_problem_setup(pars, fn, gr, eq_fn, eq_b, eq_jac, ineq_fn, ineq_lower, ineq_upper, ineq_jac, lower, upper, ...)
    num_parameters <- problem_indicators[1]
    objective_list <- solnp_objective_wrapper(fn = fn, gr = gr, n = num_parameters, ...)
    inequality_list <- solnp_ineq_wrapper(ineq_fn, ineq_jac, n_ineq = problem_indicators[5], n = problem_indicators[1], ...)
    equality_list <- solnp_eq_wrapper(eq_fn, eq_b, eq_jac, n_eq = problem_indicators[8], n = problem_indicators[1], ...)
    lower_tmp <- lower
    upper_tmp <- upper
    pars_tmp <- pars
    objective_fun <- objective_list$objective_fun
    gradient_fun <- objective_list$gradient_fun

    # Objective function checks and initial evaluation
    ineq_f <- inequality_list$ineq_f
    ineq_j <- inequality_list$ineq_j
    eq_f <- equality_list$eq_f
    eq_j <- equality_list$eq_j
    aug_fn_gr <- augmented_function_gradient(gr = objective_list$gradient_fun, n_ineq = problem_indicators[5], n = problem_indicators[1])
    aug_ineq_jac <- aug_eq_jac <- NULL
    if (!is.null(ineq_fn)) {
        aug_ineq_jac <- augmented_inequality_jacobian(ineq_jac = inequality_list$ineq_j, n_ineq = problem_indicators[5], n = problem_indicators[1])
    } else {
        aug_ineq_jac <- function(pars) {return(NULL)}
    }
    if (!is.null(eq_fn)) {
        aug_eq_jac <- augmented_equality_jacobian(eq_jac = equality_list$eq_j, n_eq = problem_indicators[8], n_ineq = problem_indicators[5], n = problem_indicators[1])
    } else {
        aug_eq_jac <- function(pars) {return(NULL)}
    }
    # Inequality constraint checks and initial evaluation
    initial_ineq_values <- NULL
    initial_ineq_guess <- NULL

    if (!is.null(ineq_f)) {
        initial_ineq_values <- ineq_f(pars_tmp)
        initial_ineq_guess <- (ineq_lower + ineq_upper) / 2
    }
    # Equality constraint checks and initial evaluation
    initial_eq_values <- NULL
    if (!is.null(eq_f)) {
        initial_eq_values <- eq_f(pars_tmp)
    }

    initial_objective_value <- objective_fun(pars_tmp)
    # Combine inequality and parameter bounds into a single matrix
    all_bounds <- rbind(cbind(ineq_lower, ineq_upper),
                        cbind(lower_tmp, upper_tmp))

    # Parse and apply control parameters
    control_params <- .solnp_ctrl(control)
    penalty_param <- control_params$rho
    max_major_iterations <- control_params$max_iter
    max_minor_iterations <- control_params$min_iter
    tolerance <- control_params$tol
    trace_progress <- control_params$trace
    # Calculate total number of constraints
    n_ineq <- problem_indicators[5]
    n_eq <- problem_indicators[8]
    total_constraints <- n_ineq + n_eq

    # Initialize objective function values and internal status vector
    current_objective_value <- previous_objective_value <- initial_objective_value
    status_vector <- 0 * .ones(3, 1) # [obj_change, constraint_norm_prev, constraint_norm_curr]


    # Initialize Lagrange multipliers and constraints
    if (total_constraints > 0) {
        lagrange_mults <- 0 * .ones(total_constraints, 1)
        # Concatenate equality and inequality constraint values
        constraint_values <- c(initial_eq_values, initial_ineq_values)
        # Adjust inequality constraint values based on initial guess
        if (problem_indicators[4]) {
            # Calculate slack for inequality constraints relative to bounds
            temp_ineq_slack_lb <- constraint_values[(n_eq + 1):total_constraints] - ineq_lower
            temp_ineq_slack_ub <- ineq_upper - constraint_values[(n_eq + 1):total_constraints]
            # Check if all initial inequalities are within bounds (positive slack)
            min_slack <- apply(cbind(temp_ineq_slack_lb, temp_ineq_slack_ub), 1, min)
            if (all(min_slack > 0)) {
                initial_ineq_guess <- constraint_values[(n_eq + 1):total_constraints]
            }
            # Adjust inequality constraint values by subtracting the initial guess
            constraint_values[(n_eq + 1):total_constraints] <-
                constraint_values[(n_eq + 1):total_constraints] - initial_ineq_guess
        }
        # Store the initial norm of the constraint violations
        status_vector[2] <- .vnorm(constraint_values)

        # If initial constraints are "small" (within 10*tolerance), set penalty to 0
        if (max(status_vector[2] - 10 * tolerance, n_ineq, na.rm = TRUE) <= 0) {
            penalty_param <- 0
        }
    } else {
        lagrange_mults <- 0
    }

    # Starting augmented parameter vector (initial_ineq_guess for inequalities, then original parameters)
    augmented_parameters <- c(initial_ineq_guess, pars_tmp)

    # Initialize Hessian approximation (identity matrix for augmented parameters)
    augmented_hessian <- diag(num_parameters + n_ineq)

    # Initialize mu (related to optimality conditions, often small)
    lambda <- num_parameters # Or some other initial value

    # Store initial objective and constraint values
    scaled_eval <- c(initial_objective_value, initial_eq_values, initial_ineq_values)

    # Initialize major iteration counter
    major_iteration_count <- 0

    # Store historical objective values
    historical_objective_values <- c(initial_objective_value)

    penalty_param <- control_params$rho
    tol <- control_params$tol
    min_iter <- control_params$min_iter
    trace <- control_params$trace
    # Main optimization loop (Major Iterations)
    while (major_iteration_count < max_major_iterations) {
        ctrl <- list(rho = penalty_param, min_iter = min_iter, tol = tol, trace = trace)
        major_iteration_count <- major_iteration_count + 1
        assign(".solnp_iter", major_iteration_count, envir = .solnp_environment) # Update global iteration count
        # Control parameters for the inner optimization (.subnp)
        #
        # Determine scaling factors for objective, constraints, and parameters
        objective_and_eq_scale <- NULL
        if (problem_indicators[7]) { # If equality constraints exist
            # Scale objective by its value, and equality constraints by max absolute value
            objective_and_eq_scale <- c(scaled_eval[1], rep(1, n_eq) * max(abs(scaled_eval[2:(n_eq + 1)])))
        } else {
            objective_and_eq_scale <- 1 # No scaling for equality constraints
        }

        if (!problem_indicators[11]) { # If no parameter or inequality bounds
            scaling_factors <- c(objective_and_eq_scale, augmented_parameters)
        } else {
            # If bounds exist, scale parameters by 1 (or a default value)
            scaling_factors <- c(objective_and_eq_scale, rep(1, length(augmented_parameters)))
        }

        # Ensure scaling factors are not too small or too large
        scaling_factors <- apply(matrix(scaling_factors, ncol = 1), 1,
                                 FUN = function(x) min(max(abs(x), tolerance), 1/tolerance))
        L <- list()
        L$objective_fun <- objective_fun
        L$equality_fun <- eq_f
        L$inequality_fun <- ineq_f
        L$augmented_gradient <- aug_fn_gr
        L$augmented_inequality_jac <- aug_ineq_jac
        L$augmented_equality_jac <- aug_eq_jac
        L$gradient_fun <- gradient_fun
        L$ineq_jac <- ineq_j
        L$eq_jac <- eq_j
        L$lower <- lower_tmp
        L$upper <- upper_tmp
        L$inequality_lower <- ineq_lower
        L$inequality_upper <- ineq_upper
        L$problem_indicators <- problem_indicators
        .env <- new.env()
        subnp_results <- .rsubnp(pars = augmented_parameters, lagrange_mults = lagrange_mults, scaled_eval = scaled_eval,
                                augmented_hessian = augmented_hessian, lambda = lambda, scaling_factors = scaling_factors,
                                ctrl = ctrl, fun_list = L, .env = .solnp_environment)
        if (get(".solnp_errors", envir = .solnp_environment) == 1) {
            max_major_iterations <- major_iteration_count # Terminate major iterations
        }

        # Update parameters, Lagrange multipliers, Hessian, and mu from subnp results
        augmented_parameters <- subnp_results$p
        lagrange_mults <- subnp_results$y
        augmented_hessian <- subnp_results$augmented_hessian
        lambda <- subnp_results$lambda # mu is returned as lambda from subnp

        # Extract original parameters from augmented parameters
        current_parameters <- augmented_parameters[(n_ineq + 1):(n_ineq + num_parameters)]

        # Re-evaluate objective function with updated parameters
        current_objective_value <- objective_fun(current_parameters)

        # Report progress if tracing is enabled
        if (trace_progress) {
            .report(major_iteration_count, current_objective_value, current_parameters)
        }

        # Re-evaluate constraint functions with updated parameters
        current_eq_values <- eq_f(current_parameters)
        current_ineq_values <- ineq_f(current_parameters)

        # Update combined objective and constraint values
        scaled_eval <- c(current_objective_value, current_eq_values, current_ineq_values)

        # Calculate relative change in objective function
        status_vector[1] <- (previous_objective_value - scaled_eval[1]) / max(abs(scaled_eval[1]), 1)
        previous_objective_value <- scaled_eval[1]

        # Update constraint violation norms and penalty parameter
        if (total_constraints > 0) {
            current_constraint_violations <- scaled_eval[2:(total_constraints + 1)]

            if (problem_indicators[4]) {
                # Recalculate slack for inequality constraints relative to bounds
                temp_ineq_slack_lb <- current_constraint_violations[(n_eq + 1):total_constraints] - all_bounds[1:n_ineq, 1]
                temp_ineq_slack_ub <- all_bounds[1:n_ineq, 2] - current_constraint_violations[(n_eq + 1):total_constraints]
                # If inequalities are within bounds, update the initial guess for inequalities in augmented parameters
                if (min(c(temp_ineq_slack_lb, temp_ineq_slack_ub)) > 0) {
                    augmented_parameters[1:n_ineq] <- current_constraint_violations[(n_eq + 1):total_constraints]
                }
                # Re-adjust inequality constraint values by subtracting the current guess for inequalities
                current_constraint_violations[(n_eq + 1):total_constraints] <-
                    current_constraint_violations[(n_eq + 1):total_constraints] - augmented_parameters[1:n_ineq]
            }
            # Calculate current norm of constraint violations
            status_vector[3] <- .vnorm(current_constraint_violations)

            # Adjust penalty parameter based on constraint violation norms
            if (status_vector[3] < 10 * tolerance) {
                penalty_param <- 0
                lambda <- min(lambda, tolerance)
            }

            if (status_vector[3] < 5 * status_vector[2]) {
                penalty_param <- penalty_param / 5
            }

            if (status_vector[3] > 10 * status_vector[2]) {
                penalty_param <- 5 * max(penalty_param, sqrt(tolerance))
            }

            # If changes in objective and constraints are small, reset Lagrange multipliers and Hessian
            if (max(c(tolerance + status_vector[1], status_vector[2] - status_vector[3])) <= 0) {
                lagrange_mults <- 0
                augmented_hessian <- diag(diag(augmented_hessian)) # Reset Hessian to diagonal
            }
            status_vector[2] <- status_vector[3] # Update previous constraint norm
        }

        # Check for convergence based on combined norms of objective change and constraint violations
        if (.vnorm(c(status_vector[1], status_vector[2])) <= tolerance) {
            max_major_iterations <- major_iteration_count # Terminate major iterations
        }

        # Store current objective value for historical tracking
        historical_objective_values <- c(historical_objective_values, current_objective_value)
    }

    # Final adjustments and results preparation
    if (problem_indicators[4]) { # If inequality constraints were present
        initial_ineq_guess <- augmented_parameters[1:n_ineq]
    }

    # Extract the optimized parameters (excluding initial inequality guess)
    optimized_parameters <- augmented_parameters[(n_ineq + 1):(n_ineq + num_parameters)]
    # Determine convergence status
    convergence_status <- 0
    if (get(".solnp_errors", envir = .solnp_environment) == 1) {
        convergence_status <- 2 # Solution not reliable due to Hessian inversion problems
        if (trace_progress) {
            cat("\nsolnp: Solution not reliable....Problem Inverting Hessian.\n")
        }
    } else {
        if (.vnorm(c(status_vector[1], status_vector[2])) <= tolerance) {
            convergence_status <- 0 # Converged successfully
            if (trace_progress) {
                cat(paste("\nsolnp: Completed in ", major_iteration_count, " iterations\n", sep = ""))
            }
        } else {
            convergence_status <- 1 # Exited after maximum iterations, tolerance not achieved
            if (trace_progress) {
                cat(paste("\nsolnp: Exiting after maximum number of iterations. Tolerance not achieved\n", sep = ""))
            }
        }
    }
    # Calculate elapsed time and retrieve final function evaluation count
    final_function_evaluations <- get(".solnp_nfn", envir = .solnp_environment)
    elapsed_time <- Sys.time() - start_time

    # Prepare results list
    results <- list(
        pars = optimized_parameters,
        convergence = convergence_status,
        objective = tail(as.numeric(historical_objective_values), 1),
        objective_history = as.numeric(historical_objective_values),
        lagrange = lagrange_mults,
        hessian = augmented_hessian,
        ineq_initial_values = initial_ineq_guess,
        n_eval = final_function_evaluations,
        outer_iterations = major_iteration_count,
        elapsed_time = elapsed_time,
        vscale = scaling_factors
    )
    return(results)
}

.rsubnp <- function(pars, lagrange_mults, scaled_eval, augmented_hessian, lambda, scaling_factors, ctrl, fun_list, .env, ...) {
    .solnp_fun <- fun_list$objective_fun      # J(P)
    .solnp_eqfun <- fun_list$equality_fun       # EC(P)
    .solnp_ineqfun <- fun_list$inequality_fun   # IC(P)
    .solnp_gradfun <- fun_list$gradient_fun
    .solnp_eqjac   <- fun_list$eq_jac
    .solnp_ineqjac <- fun_list$ineq_jac
    lower <- fun_list$lower
    upper <- fun_list$upper
    ineq_lower <- fun_list$inequality_lower
    ineq_upper <- fun_list$inequality_upper
    setup <- fun_list$problem_indicators

    # pars [nineq + np]
    rho <- ctrl$rho
    maxit <- ctrl$min_iter
    tol <- ctrl$tol
    trace <- ctrl$trace
    # [1] length of pars
    # [2] has function gradient?
    # [3] has hessian?
    # [4] has ineq?
    # [5] ineq length
    # [6] has jacobian (inequality)
    # [7] has eq?
    # [8] eq length
    # [9] has jacobian (equality)
    # [10] has upper / lower bounds
    # [11] has either lower/upper bounds or ineq


    n_eq <- setup[8]
    n_ineq <- setup[ 5 ]
    n_pars <- setup[ 1 ]
    status_flag <- 1
    line_search_steps <- c(0,0,0)
    n_constraints <- n_eq + n_ineq
    n_pic <- n_pars + n_ineq
    working_params <- pars
    # param_bounds [ 2 x (nineq + np) ]
    param_bounds <- rbind(cbind(ineq_lower, ineq_upper), cbind(lower, upper))
    step_objectives <- numeric()
    param_trials <- matrix()
    bfgs_scaling_factors <- numeric()

    if (n_constraints == 0) {
        augmented_jacobian <- matrix(0, 0, n_pars)  # avoid fake constraint
    }
    # scale the cost, the equality constraints, the inequality constraints,
    # the parameters (inequality parameters AND actual parameters),
    # and the parameter bounds if there are any
    # Also make sure the parameters are no larger than (1-tol) times their bounds
    # scaled_eval [ 1 neq nineq]

    scaled_eval <- scaled_eval / scaling_factors[1:(n_constraints + 1)]
    # working_params [np]
    working_params <- working_params / scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]
    if (setup[11] > 0) {
        if (setup[10] == 0) {
            mm <- n_ineq
        } else {
            mm <- n_pic
        }
        param_bounds <- param_bounds / cbind(scaling_factors[(n_eq + 2):(n_eq + mm + 1)], scaling_factors[(n_eq + 2):(n_eq + mm + 1)])
    }

    # scale the lagrange multipliers and the Hessian
    if (n_constraints > 0) {
        # lagrange_mults [total constraints = nineq + neq]
        # scale here is [tc] and dot multiplied by lagrange_mults
        lagrange_mults <- scaling_factors[2:(n_constraints + 1)] * lagrange_mults / scaling_factors[1]
    }
    # lagrange_mults = [zeros 3x1]

    # h is [ (np+nineq) x (np+nineq) ]
    #columnvector %*% row vector (size h) then dotproduct h then dotdivide scale[1]
    augmented_hessian <- augmented_hessian * (scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)] %*% t(scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)])) / scaling_factors[1]
    # h[ 8x8 eye]
    j <- scaled_eval[1]
    if (setup[7] > 0 && setup[4] > 0) {
        # Both equality and inequality constraints
        Aeq_block   <- cbind(matrix(0, nrow = n_eq, ncol = n_ineq), matrix(0, nrow = n_eq, ncol = n_pars))
        Aineq_block <- cbind(-diag(n_ineq), matrix(0, nrow = n_ineq, ncol = n_pars))
        augmented_jacobian <- rbind(Aeq_block, Aineq_block)
    } else if (setup[4] > 0) {
        # Only inequality constraints
        augmented_jacobian <- cbind(-diag(n_ineq), matrix(0, nrow = n_ineq, ncol = n_pars))
    } else if (setup[7] > 0) {
        # Only equality constraints
        augmented_jacobian <- cbind(matrix(0, nrow = n_eq, ncol = n_ineq), matrix(0, nrow = n_eq, ncol = n_pars))
    } else {
        # No constraints — explicitly handle this
        augmented_jacobian <- matrix(0, 0, n_ineq + n_pars)
    }
    # gradient
    grad_augmented_obj <- 0 * .ones(n_pic, 1)
    p <- working_params[1:n_pic]
    if (n_constraints > 0) {
        # [ nc ]
        constraint <- scaled_eval[2:(n_constraints + 1)]
        # constraint [5 0 11 3x1]
        # gradient routine
        x <- working_params[(n_ineq + 1):n_pic] * scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)]
        if (setup[7] > 0) {
            Aeq <- .solnp_eqjac(x)
            Aeq <- Aeq * outer(1 / scaling_factors[2:(n_eq + 1)], scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)])
            augmented_jacobian[1:n_eq, (n_ineq + 1):n_pic] <- Aeq
        }
        if (setup[4] > 0) {
            Aineq <- .solnp_ineqjac(x)
            Aineq <- Aineq * outer(1 / scaling_factors[(n_eq + 2):(n_constraints + 1)], scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)])
            augmented_jacobian[(n_eq + 1):n_constraints, (n_ineq + 1):n_pic] <- Aineq
            augmented_jacobian[(n_eq + 1):n_constraints, 1:n_ineq] <- -diag(n_ineq)
        }
        g_user <- .solnp_gradfun(x)
        g_user <- g_user / scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)] / scaling_factors[1]
        grad_augmented_obj[(n_ineq + 1):n_pic, 1] <- g_user


        if (setup[4] > 0) {
            constraint[(n_eq + 1):(n_eq + n_ineq)] <- constraint[(n_eq + 1):(n_eq + n_ineq)] - working_params[1:n_ineq]
        }
        # solver messages
        if (.solvecond(augmented_jacobian) > 1 / .eps ) {
            if (trace) .solver_warnings("M1")
        }

        # augmented_jacobian(matrix) x columnvector - columnvector
        # constraint_offset [nc,1]
        constraint_offset <- augmented_jacobian %*% working_params - constraint
        status_flag <- -1
        line_search_steps[1] = tol - max(abs(constraint))

        if (line_search_steps[1] <= 0) {
            status_flag <- 1
            if (setup[11] == 0) {
                # augmented_jacobian %*% t(augmented_jacobian) gives [nc x nc]
                # t(augmented_jacobian) %*% above gives [(np+nc) x 1]
                working_params <- working_params - t(augmented_jacobian) %*% solve(augmented_jacobian %*% t(augmented_jacobian), constraint)
                line_search_steps[1] <- 1
            }
        }

        if (line_search_steps[1] <= 0) {
            # this expands working_params to [nc+np+1]
            working_params[n_pic + 1] <- 1
            augmented_jacobian <- cbind(augmented_jacobian, -constraint)
            # cx is rowvector
            cx <- cbind(.zeros(1, n_pic), 1)
            dx <- .ones(n_pic + 1, 1)
            step_interval_width <- 1
            minit <- 0
            while (step_interval_width >= tol) {
                minit <- minit + 1
                # gap [(nc + np) x 2]
                gap <- cbind(working_params[1:mm] - param_bounds[, 1], param_bounds[, 2] - working_params[1:mm])
                # this sorts every row
                gap <- t(apply(gap, 1, FUN = function(x) sort(x)))
                dx[1:mm] <- gap[, 1]
                # expand dx by 1
                dx[n_pic + 1] <- working_params[n_pic + 1]
                if (setup[10] == 0) {
                    dx[(mm + 1):n_pic] <- max(c(dx[1:mm], 100)) * .ones(n_pic - mm, 1)
                }
                # t( augmented_jacobian %*% diag( as.numeric(dx) ) ) gives [(np+nc + 1 (or more) x nc]
                # dx * t(cx) dot product of columnvectors
                # qr.solve returns [nc x 1]
                y <- try(qr.solve(t(augmented_jacobian %*% diag(as.numeric(dx), length(dx), length(dx))), dx * t(cx)), silent = TRUE)
                if (inherits(y, "try-error")) {
                    p <- working_params * scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]  # unscale the parameter vector
                    if (n_constraints > 0) {
                        y <- 0 # unscale the lagrange multipliers
                    }
                    augmented_hessian <- scaling_factors[1] * augmented_hessian / (scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)] %*% t(scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1) ]))
                    ans = list(p = p, y = y, augmented_hessian = augmented_hessian, lambda = lambda)
                    assign(".solnp_errors", 1, envir = .env)
                    return(ans)
                }
                    v <- dx * (dx * (t(cx) - t(augmented_jacobian) %*% y))
                if (v[n_pic + 1] > 0 ) {
                    z <- working_params[n_pic + 1] / v[n_pic + 1]
                    for (i in 1:mm) {
                        if (v[i] < 0) {
                            z <- min(z, -(param_bounds[i, 2] - working_params[i]) / v[i])
                        } else if (v[i] > 0) {
                            z <- min(z, (working_params[i] - param_bounds[i , 1]) / v[i])
                        }
                    }
                    if (z >= working_params[n_pic + 1] / v[n_pic + 1]) {
                        working_params <- working_params - z * v
                    } else {
                        working_params <- working_params - 0.9 * z * v
                    }
                    step_interval_width <- working_params[n_pic + 1]
                    if (minit >= 10) {
                        step_interval_width <- 0
                    }
                } else {
                    step_interval_width <- 0
                    minit <- 10
                }
            }
            if (minit >= 10) {
                if (trace) .solver_warnings("M2")
            }
            augmented_jacobian <- matrix(augmented_jacobian[, 1:n_pic], ncol = n_pic)
            constraint_offset <- augmented_jacobian %*% working_params[1:n_pic]
        }
    }

    p <- working_params[1:n_pic]
    y <- 0

    if (status_flag > 0) {
        tmp_value <- p[(n_ineq + 1):n_pic] * scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)]
        fun_value <- .solnp_fun(tmp_value)
        eq_value <- .solnp_eqfun(tmp_value)
        ineq_value <- .solnp_ineqfun(tmp_value)
        ctmp <- get(".solnp_nfn", envir =  .env)
        assign(".solnp_nfn", ctmp + 1, envir = .env)
        scaled_eval <- c(fun_value, eq_value, ineq_value) / scaling_factors[1:(n_constraints + 1)]
    }
    j <- scaled_eval[1]

    if (setup[4] > 0) {
        scaled_eval[(n_eq + 2):(n_constraints + 1)] <- scaled_eval[(n_eq + 2):(n_constraints + 1)] - p[1:n_ineq]
    }

    if (n_constraints > 0) {
        scaled_eval[2:(n_constraints + 1)] <- scaled_eval[2:(n_constraints + 1)] - augmented_jacobian %*% p + constraint_offset
        j <- scaled_eval[1] - t(lagrange_mults) %*% matrix(scaled_eval[2:(n_constraints + 1)], ncol = 1) + rho * .vnorm(scaled_eval[2:(n_constraints + 1)])^2
    }

    minit <- 0
    delta_position <- rep(0, n_pic)   # <- NEW
    delta_gradient <- rep(0, n_pic)   # <- NEW
    reduce <- Inf
    while (minit < maxit) {
        minit <- minit + 1
        if (status_flag > 0) {
            # Unscale actual parameters
            Alin_x <- augmented_jacobian[ , (n_ineq + 1):(n_ineq + n_pars), drop = FALSE]
            x <- p[(n_ineq + 1):n_pic] * scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)]
            # Evaluate f(x) gradient
            fun_value <- .solnp_fun(x)
            ineq_value <- .solnp_ineqfun(x)
            eq_value <- .solnp_eqfun(x)
            grad_f <- .solnp_gradfun(x)
            grad_f.p <- grad_f * scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)]
            # Evaluate Jacobians (they may return NULL)
            obm <- c(fun_value, eq_value, ineq_value) / scaling_factors[1:(n_constraints + 1)]

            Aeq.j   <- matrix(0, 0, n_pars)
            Aineq.j <- matrix(0, 0, n_pars)
            gradx <- grad_f.p / scaling_factors[1]   # only ∇f term
            grads <- numeric(0)

            # Scale Jacobians
            if (n_eq > 0) {
                Aeq.raw <- .solnp_eqjac(x)
                if (!is.null(Aeq.raw)) {
                    #Aeq <- Aeq * outer(1 / scaling_factors[2:(n_eq + 1)], scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)])
                    Aeq.p <- sweep(Aeq.raw, 2, scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)], '*')
                    Aeq.j <- sweep(Aeq.p, 1, scaling_factors[2:(n_eq + 1)], '/')
                }
            }

            if (n_ineq) {
                Aineq.raw <- .solnp_ineqjac(x)
                if (!is.null(Aineq.raw)) {
                    Aineq.p <- sweep(Aineq.raw, 2, scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)], '*')
                    Aineq.j <- sweep(Aineq.p, 1, scaling_factors[(n_eq + 2):(n_constraints + 1)], '/')
                }
            }

            J <- rbind(Aeq.j, Aineq.j) - Alin_x

            if (n_constraints > 0) {
                ## residual vector r  (length n_constraints)
                r <- obm[2:(n_constraints + 1)]

                if (n_ineq > 0) {
                    r[(n_eq + 1):n_constraints] <- r[(n_eq + 1):n_constraints] - p[1:n_ineq]
                }
                r <- r - augmented_jacobian %*% p + constraint_offset

                ## multiplier & penalty pieces (only if J has rows)
                if (nrow(J) > 0) {
                    JTy <- drop(crossprod(J, lagrange_mults))   # length n_pars
                    JTr <- drop(crossprod(J, r))                # length n_pars
                } else {
                    JTy <- 0
                    JTr <- 0
                }
                gradx <- gradx - JTy + 2 * rho * JTr
                ## slack variables
                if (n_ineq > 0) {
                    idx   <- (n_eq + 1):n_constraints
                    rg    <- r[idx]                   # g̃ − s
                    grads <- lagrange_mults[idx] - 2 * rho * rg
                }
            }
            if (n_ineq > 0)  grad_augmented_obj[1:n_ineq] <- grads
            grad_augmented_obj[(n_ineq + 1):(n_ineq + n_pars)] <- gradx
        }
        if (setup[4] > 0) {
            grad_augmented_obj[1:n_ineq] <- 0
        }

        if (minit > 1) {
            delta_gradient <- grad_augmented_obj - delta_gradient
            delta_position <- p - delta_position
            bfgs_scaling_factors[1] <- t(delta_position) %*% augmented_hessian %*% delta_position
            bfgs_scaling_factors[2] <- t(delta_position) %*% delta_gradient
            if ((bfgs_scaling_factors[1] * bfgs_scaling_factors[2]) > 0) {
                delta_position <- augmented_hessian %*% delta_position
                augmented_hessian <- augmented_hessian - (delta_position %*% t(delta_position)) / bfgs_scaling_factors[1] + (delta_gradient %*% t(delta_gradient)) / bfgs_scaling_factors[2]
            }
        }

        dx <- 0.01 * .ones(n_pic, 1)
        if (setup[11] > 0) {
            gap <- cbind(p[1:mm] - param_bounds[, 1], param_bounds[, 2] - p[1:mm])
            gap <- t(apply(gap, 1, FUN = function(x) sort(x)))
            gap <- gap[, 1] + sqrt(.eps) * .ones(mm, 1)
            dx[1:mm, 1] <- .ones(mm, 1) / gap
            if (setup[10] == 0) {
                dx[(mm + 1):n_pic, 1] <- min(c(dx[1:mm, 1], 0.01)) * .ones(n_pic - mm, 1)
            }
        }
        step_interval_width <- -1
        lambda <- lambda / 10

        while (step_interval_width <= 0) {
            cz <- try(chol(augmented_hessian + lambda * diag(as.numeric(dx * dx), length(dx), length(dx))), silent = TRUE)
            if (inherits(cz, "try-error")) {
                p <- p * scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]  # unscale the parameter vector
                if (n_constraints > 0) {
                    y <- 0 # unscale the lagrange multipliers
                }
                augmented_hessian <- scaling_factors[1] * augmented_hessian / (scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)] %*% t(scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]))
                ans <- list(p = p, y = y, augmented_hessian = augmented_hessian, lambda = lambda)
                assign(".solnp_errors", 1, envir = .env)
                return(ans)
            }
            cz <- try(solve(cz), silent = TRUE)
            if (inherits(cz, "try-error")) {
                p <- p * scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]  # unscale the parameter vector
                if (n_constraints > 0) {
                    y <- 0 # unscale the lagrange multipliers
                }
                augmented_hessian <- scaling_factors[1] * augmented_hessian / (scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)] %*% t(scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]))
                ans <- list(p = p, y = y, augmented_hessian = augmented_hessian, lambda = lambda)
                assign(".solnp_errors", 1, envir = .env)
                return(ans)
            }
            delta_gradient <- t(cz) %*% grad_augmented_obj

            if (n_constraints == 0 || nrow(augmented_jacobian) == 0) {
                u <- -cz %*% delta_gradient
            } else {
                y <- try(qr.solve(t(cz) %*% t(augmented_jacobian), delta_gradient), silent = TRUE)
                if (inherits(y, "try-error")) {
                    p <- p * scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]  # unscale the parameter vector
                    if (n_constraints > 0) {
                        # y = scaling_factors[ 1 ] * y / scaling_factors[ 2:(nc + 1) ] # unscale the lagrange multipliers
                        y <- 0
                    }
                    augmented_hessian <- scaling_factors[1] * augmented_hessian / (scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)] %*% t(scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]))
                    ans <- list(p = p, y = y, augmented_hessian = augmented_hessian, lambda = lambda)
                    assign(".solnp_errors", 1, envir = .env)
                    return(ans)
                }
                u <- -cz %*% (delta_gradient - (t(cz) %*% t(augmented_jacobian)) %*% y)
            }
            working_params <- u[1:n_pic] + p

            if (setup[11] == 0) {
                step_interval_width <- 1
            } else {
                step_interval_width <- min(c(working_params[1:mm] - param_bounds[, 1], param_bounds[, 2] - working_params[1:mm]))
                lambda <- 3 * lambda
            }
        }

        line_search_steps[1] <- 0
        scaled_eval_1 <- scaled_eval
        scaled_eval_2 <- scaled_eval_1
        step_objectives[1] <- j
        step_objectives[2] <- j
        param_trials <- cbind(p, p)
        line_search_steps[3] <- 1.0
        param_trials <- cbind(param_trials, working_params)
        tmpv <- param_trials[(n_ineq + 1):n_pic, 3 ] * scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)]
        ### CONTINUE HERE
        fun_value <- .solnp_fun(tmpv)
        eq_value <- .solnp_eqfun(tmpv)
        ineq_value <- .solnp_ineqfun(tmpv)
        ctmp <- get(".solnp_nfn", envir =  .env)
        assign(".solnp_nfn", ctmp + 1, envir = .env)
        scaled_eval_3 <- c(fun_value, eq_value, ineq_value) / scaling_factors[1:(n_constraints + 1)]
        step_objectives[3] <- scaled_eval_3[1]

        if ( setup[4] > 0) {
            scaled_eval_3[(n_eq + 2):(n_constraints + 1)] <- scaled_eval_3[(n_eq + 2):(n_constraints + 1)] - param_trials[1:n_ineq, 3]
        }

        if (n_constraints > 0) {
            scaled_eval_3[2:(n_constraints + 1)] <- scaled_eval_3[2:(n_constraints + 1)] - augmented_jacobian %*% param_trials[, 3] + constraint_offset
            step_objectives[3] <- scaled_eval_3[1] - t(lagrange_mults) %*% scaled_eval_3[2:(n_constraints + 1)] + rho * .vnorm(scaled_eval_3[2:(n_constraints + 1) ])^2
        }

        step_interval_width <- 1
        while (step_interval_width > tol) {
            line_search_steps[2] <- (line_search_steps[1] + line_search_steps[3]) / 2
            param_trials[, 2] <- (1 - line_search_steps[2]) * p + line_search_steps[2] * working_params
            tmpv <- param_trials[(n_ineq + 1):n_pic, 2] * scaling_factors[(n_constraints + 2):(n_constraints + n_pars + 1)]
            fun_value <- .solnp_fun(tmpv)
            eq_value <- .solnp_eqfun(tmpv)
            ineq_value <- .solnp_ineqfun(tmpv)
            ctmp <- get(".solnp_nfn", envir =  .env)
            assign(".solnp_nfn", ctmp + 1, envir = .env)
            scaled_eval_2 <- c(fun_value, eq_value, ineq_value) / scaling_factors[1:(n_constraints + 1)]
            step_objectives[2] <- scaled_eval_2[1]
            if (setup[4] > 0) {
                scaled_eval_2[(n_eq + 2):(n_constraints + 1)] <- scaled_eval_2[(n_eq + 2):(n_constraints + 1)] - param_trials[1:n_ineq , 2]
            }
            if (n_constraints > 0) {
                scaled_eval_2[2:(n_constraints + 1)] <- scaled_eval_2[2:(n_constraints + 1)] - augmented_jacobian %*% param_trials[, 2] + constraint_offset
                step_objectives[2] <- scaled_eval_2[1] - t(lagrange_mults) %*% scaled_eval_2[2:(n_constraints + 1)] + rho * .vnorm(scaled_eval_2[2:(n_constraints + 1) ])^2
            }
            max_step_obj <- max(step_objectives)

            if (max_step_obj < j) {
                min_step_obj <- min(step_objectives)
                step_interval_width <- tol * (max_step_obj - min_step_obj) / (j - max_step_obj)
            }

            condif1 <- step_objectives[2] >= step_objectives[1]
            condif2 <- step_objectives[1] <= step_objectives[3] && step_objectives[2] < step_objectives[1]
            condif3 <- step_objectives[2] <  step_objectives[1] && step_objectives[1] > step_objectives[3]

            if (condif1) {
                step_objectives[3] <- step_objectives[2]
                scaled_eval_3 <- scaled_eval_2
                line_search_steps[3] <- line_search_steps[2]
                param_trials[, 3] <- param_trials[, 2]
            }

            if (condif2) {
                step_objectives[3] <- step_objectives[2]
                scaled_eval_3 <- scaled_eval_2
                line_search_steps[3] <- line_search_steps[2]
                param_trials[, 3] <- param_trials[, 2]
            }

            if (condif3) {
                step_objectives[1] <- step_objectives[2]
                scaled_eval_1 <- scaled_eval_2
                line_search_steps[1] <- line_search_steps[2]
                param_trials[, 1] <- param_trials[, 2]
            }

            if (step_interval_width >= tol) {
                step_interval_width <- line_search_steps[3] - line_search_steps[1]
            }

        }

        delta_position <- p
        delta_gradient <- grad_augmented_obj
        status_flag <- 1
        min_step_obj <- min(step_objectives)
        if (j <= min_step_obj) {
            maxit <- minit
        }
        reduce <- (j - min_step_obj) / (1 + abs(j))

        if (reduce < tol) {
            maxit <- minit
        }

        condif1 <- step_objectives[1] <  step_objectives[2]
        condif2 <- step_objectives[3] <  step_objectives[2] && step_objectives[1] >= step_objectives[2]
        condif3 <- step_objectives[1] >= step_objectives[2] && step_objectives[3] >= step_objectives[2]

        if (condif1) {
            j <- step_objectives[1]
            p <- param_trials[, 1]
            scaled_eval <- scaled_eval_1
        }

        if (condif2) {
            j <- step_objectives[3]
            p <- param_trials[, 3]
            scaled_eval <- scaled_eval_3
        }

        if (condif3) {
            j <- step_objectives[2]
            p <- param_trials[, 2]
            scaled_eval <- scaled_eval_2
        }
    }

    p <- p * scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]  # unscale the parameter vector
    if (n_constraints > 0) {
        y <- scaling_factors[1] * y / scaling_factors[2:(n_constraints + 1)] # unscale the lagrange multipliers
    }
    augmented_hessian <- scaling_factors[1] * augmented_hessian / (scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1) ] %*% t(scaling_factors[(n_eq + 2):(n_constraints + n_pars + 1)]))

    if (reduce > tol) {
        if (trace) .solver_warnings("M3")
    }

    ans = list(p = p, y = y, augmented_hessian = augmented_hessian, lambda = lambda)
    return(ans)
}

