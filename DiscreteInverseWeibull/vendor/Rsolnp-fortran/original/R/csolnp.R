#' Nonlinear optimization using augmented Lagrange method (C++ version)
#'
#' @param pars an numeric vector of decision variables (length n).
#' @param fn the objective function (must return a scalar).
#' @param gr an optional function for computing the analytic gradient of the function (must return
#' a vector of length n).
#' @param eq_fn an optional function for calculating equality constraints.
#' @param eq_b a vector of the equality bounds (if eq_fn provided).
#' @param eq_jac an optional function for computing the analytic jacobian of the equality.
#' function (a matrix with number of columns n and number of rows the same length as the number of equalities).
#' @param ineq_fn an optional function for calculating inequality constraints.
#' @param ineq_lower the lower bounds for the inequality (must be finite)
#' @param ineq_upper the upper bounds for the inequalitiy (must be finite)
#' @param ineq_jac an optional function for computing the analytic jacobian of the inequality (a matrix
#' with number of columns n and number of rows the same length as the number of inequalities).
#' @param lower lower bounds for the parameters. This is strictly required.
#' @param upper upper bounds for the parameters. This is strictly required.
#' @param control a list of solver control parameters (see details).
#' @param ... additional arguments passed to the supplied functions (common to all functions supplied).
#' @param use_r_version (logical) used for debugging and validation. Uses the R version of the solver
#' rather than the C++ version. Will be deprecated in future releases.
#' @returns A list with the following slot:
#' \describe{
#'   \item{pars}{ The parameters at the optimal solution found.}
#'   \item{objective}{ The value of the objective at the optimal solution found.}
#'   \item{objective_history}{ A vector of objective values obtained at each outer iteration.}
#'   \item{out_iterations}{The number of outer iterations used to arrive at the solution.}
#'   \item{convergence}{The convergence code (0 = converged).}
#'   \item{message}{The convergence message.}
#'   \item{kkt_diagnostics}{A list of optimal solution diagnostics.}
#'   \item{lagrange}{ The vector of Lagrange multipliers at the optimal solution found.}
#'   \item{n_eval}{ The number of function evaluations.}
#'   \item{elapsed}{ The time taken to find a solution.}
#'   \item{hessian}{ The Hessian at the optimal solution.}
#' }
#' @details
#' The optimization problem solved by \code{csolnp} is formulated as:
#' \deqn{
#' \begin{aligned}
#' \min_{x \in \mathbb{R}^n} \quad & f(x) \\
#' \text{s.t.} \quad & g(x) = b \\
#' & h_l \le h(x) \le h_u\\
#' & x_l \le x \le x_u\\
#' \end{aligned}
#' }
#'
#' where \eqn{f(x)} is the objective function, \eqn{g(x)} is the vector of equality constraints
#' with target value \eqn{b}, \eqn{h(x)} is the vector of inequality constraints bounded
#' by \eqn{h_l} and \eqn{h_u}, with parameter bounds \eqn{x_l} and \eqn{x_u}. Internally,
#' inequality constraints are converted into equality constraints using slack variables
#' and solved using an augmented Lagrangian approach.
#' This function is based on the original R code, but converted to C++, making use of
#' \code{Rcpp} and \code{RcppArmadillo}.
#' Additionally, it allows the user to pass in analytic gradient and Jacobians, else
#' finite differences using functions from the \code{numDeriv} package are used.
#'
#' The control list consists of the following options:
#' \describe{
#'   \item{rho}{Numeric. Initial penalty parameter for the augmented Lagrangian. Controls the weight given to constraint violation in the objective. Default is \code{1}.}
#'   \item{max_iter}{Integer. Maximum number of major (outer) iterations allowed. Default is \code{400}.}
#'   \item{min_iter}{Integer. Maximum number of minor (inner) iterations (per major iteration) for the quadratic subproblem solver. Default is \code{800}.}
#'   \item{tol}{Numeric. Convergence tolerance for both feasibility (constraint violation) and optimality (change in objective). The algorithm terminates when changes fall below this threshold. Default is \code{1e-8}.}
#'   \item{trace}{Integer If \code{1}, prints progress, \code{2} includes diagnostic information during optimization. Default is \code{0}.}
#' }
#'
#' Tracing information provides the following:
#' \describe{
#'   \item{Iter}{The current major iteration number.}
#'   \item{Obj}{The value of the objective function \eqn{f(x)} at the current iterate.}
#'   \item{||Constr||}{The norm of the current constraint violation, summarizing how well all constraints (equality and inequality) are satisfied. Typically the Euclidean or infinity norm.}
#'   \item{RelObj}{The relative change in the objective function value compared to the previous iteration, i.e., \eqn{|f_k - f_{k-1}| / max(1, |f_{k-1}|)}.}
#'   \item{Step}{The norm of the parameter update taken in this iteration, i.e., \eqn{||x_k - x_{k-1}||}.}
#'   \item{Penalty}{The current value of the penalty parameter (\eqn{\rho}) in the augmented Lagrangian. This parameter is adaptively updated to balance objective minimization and constraint satisfaction.}
#' }
#' @rdname csolnp
#' @author Alexios Galanos
#' @export
csolnp <- function(pars, fn, gr = NULL, eq_fn = NULL, eq_b = NULL, eq_jac = NULL,
                   ineq_fn = NULL, ineq_lower = NULL, ineq_upper = NULL,
                   ineq_jac = NULL, lower = NULL, upper = NULL,
                   control = list(), use_r_version = FALSE, ...)
{
  # Start timer for performance measurement
  if (use_r_version) {
    return(.rsolnp(pars, fn = fn, gr = gr, eq_fn = eq_fn, eq_b = eq_b, eq_jac = eq_jac,
                   ineq_fn = ineq_fn, ineq_lower = ineq_lower, ineq_upper = ineq_upper,
                   ineq_jac = ineq_jac, lower = lower, upper = upper, control = control, ...))
  }
  start_time <- Sys.time()
  # Store original parameter names
  parameter_names <- names(pars)
  # Indicator vector for problem characteristics
  # [1] Number of parameters (np)
  # [4] Has inequality constraints
  # [5] Number of inequality constraints
  # [7] Has equality constraints
  # [8] Number of equality constraints
  # [10] Has parameter upper/lower bounds
  # [11] Has any bounds (parameter or inequality)
  #
  if (is.null(lower)) {
    lower <- -1000 * abs(pars)
    warning("\nlower values are NULL. Setting to -1000 x abs(pars).")
  }

  if (is.null(upper)) {
    upper <- 1000 * abs(pars)
    warning("\nupper values are NULL. Setting to 1000 x abs(pars).")
  }

  # if pars are at their limits adjust
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

  opt_state <- list()
  opt_state$major_iteration_count <- major_iteration_count
  opt_state$max_major_iterations <- max_major_iterations
  opt_state$penalty_param <- penalty_param
  opt_state$min_iter <- max_minor_iterations
  opt_state$tol <- tol
  opt_state$ftol <- 1e-8
  opt_state$trace <- trace
  opt_state$n_eq <- n_eq
  opt_state$n_ineq <- n_ineq
  opt_state$num_parameters <- num_parameters
  opt_state$total_constraints <- total_constraints
  opt_state$problem_indicators <- problem_indicators
  opt_state$scaled_eval <- scaled_eval
  opt_state$augmented_parameters <- augmented_parameters
  opt_state$lagrange_mults <- lagrange_mults
  opt_state$augmented_hessian <- augmented_hessian
  opt_state$lambda <- lambda
  opt_state$lower_tmp <- lower
  opt_state$upper_tmp <- upper
  opt_state$ineq_lower <- ineq_lower
  opt_state$ineq_upper <- ineq_upper
  opt_state$all_bounds <- all_bounds
  opt_state$solnp_fun <- objective_fun
  opt_state$solnp_gradfun <- gradient_fun
  opt_state$solnp_eqfun <- eq_f
  opt_state$solnp_ineqfun <- ineq_f
  opt_state$solnp_eqjac <- eq_j
  opt_state$solnp_ineqjac <- ineq_j
  opt_state$previous_objective_value <- previous_objective_value
  opt_state$status_vector <- status_vector
  opt_state$historical_objective_values <- historical_objective_values
  result <- .csolnp(opt_state)
  # Calculate elapsed time and retrieve final function evaluation count
  elapsed_time <- Sys.time() - start_time

  optimized_parameters <- result$parameters

  if (problem_indicators[4] == 1) {
    initial_ineq_guess <- optimized_parameters[1:n_ineq]
  }
  optimized_parameters <- optimized_parameters[(n_ineq + 1):(n_ineq + num_parameters)]
  lagrange_mults <- result$lagrange_mults
  augmented_hessian <- result$augmented_hessian
  lambda <- result$lambda
  status_vector <- result$status_vector
  historical_objective_values <- result$historical_objective_values
  # Assign original names to optimized parameters
  names(optimized_parameters) <- parameter_names
  n_function_eval <- result$n_fun_evaluations
  major_iteration_count <- result$major_iteration_count
  convergence <- result$convergence
  convergence_message <- switch(as.character(convergence),
                                "0" = "Solution converged within tolerance.",
                                "1" = "Exiting after maximum number of iterations, tolerance not achieved.",
                                "2" = "Solution not reliable.")
  if (trace > 1) {
      print(convergence_message)
  }

  # Prepare results list
  results <- list(
    pars = optimized_parameters,
    convergence = result$convergence,
    message = convergence_message,
    objective = result$best_objective,
    objective_history = as.numeric(historical_objective_values),
    lagrange = lagrange_mults,
    hessian = augmented_hessian,
    ineq_initial_values = initial_ineq_guess,
    n_eval = n_function_eval,
    outer_iterations = major_iteration_count,
    kkt_diagnostics = result$kkt_diagnostics,
    elapsed_time = elapsed_time
  )
  return(results)
}


#' Summarize KKT Condition Diagnostics
#'
#' Given a list of KKT diagnostic statistics (stationarity, feasibility, complementarity, etc.),
#' this function prints a clear summary indicating which KKT conditions are satisfied at a specified tolerance.
#'
#' @param kkt A named list containing numeric entries for KKT diagnostics.
#'        Required names are \code{"kkt_stationarity"}, \code{"eq_violation"},
#'        \code{"ineq_violation"}, \code{"dual_feas_violation"}, \code{"compl_slackness"}.
#' @param tol Numeric tolerance for considering a condition as "satisfied". Default is \code{1e-8}.
#'
#' @return An object of class "solnp_kkt_summary" (a data.frame with columns: condition, value, status, tol).
#' @examples
#' kkt <- list(
#'   kkt_stationarity = 5.828909e-06,
#'   eq_violation = 0,
#'   ineq_violation = 0,
#'   dual_feas_violation = 0.4380053,
#'   compl_slackness = 0
#' )
#' kkt_diagnose(kkt, tol = 1e-8)
#'
#' @export
kkt_diagnose <- function(kkt, tol = 1e-8) {
  stopifnot(is.list(kkt))
  required_names <- c(
    "kkt_stationarity",
    "eq_violation",
    "ineq_violation",
    "dual_feas_violation",
    "compl_slackness"
  )
  missing <- setdiff(required_names, names(kkt))
  if (length(missing) > 0)
    stop("KKT diagnostic list is missing required element(s): ", paste(missing, collapse = ", "))

  nulls <- vapply(kkt[required_names], is.null, logical(1))
  if (any(nulls))
    stop("KKT diagnostic contains NULL for: ", paste(required_names[nulls], collapse = ", "))

  vals <- as.numeric(unlist(kkt[required_names]))
  na_idx <- which(is.na(vals))
  if (length(na_idx) > 0)
    stop("KKT diagnostic contains NA for: ", paste(required_names[na_idx], collapse = ", "))

  res <- data.frame(
    condition = c("stationarity", "equality_violation", "inequality_violation", "dual_feasibility", "complementarity"),
    value = as.numeric(unlist(kkt[required_names])),
    status = ifelse(abs(as.numeric(unlist(kkt[required_names]))) <= tol, "OK", "FAIL"),
    stringsAsFactors = FALSE
  )
  attr(res, "tolerance") <- tol
  class(res) <- c("solnp_kkt_summary", class(res))
  return(res)
}

#' @export
print.solnp_kkt_summary <- function(x, ...) {
  cat("KKT summary (tol =", attr(x, "tolerance"), "):\n")
  for (i in seq_len(nrow(x))) {
    cat(" -", x$condition[i], ":",
        if (x$status[i] == "OK") "OK" else "FAIL",
        sprintf("(|value| = %.2e)", abs(x$value[i])), "\n")
  }
  invisible(x)
}
