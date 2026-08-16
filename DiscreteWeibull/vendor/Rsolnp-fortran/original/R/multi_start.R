#' Multi-start version of csolnp
#'
#' Runs the csolnp solver from multiple diverse feasible starting values and
#' returns the best solution found.
#'
#' This function automates the process of generating multiple feasible starting
#' points (using lower, upper, and constraint information), runs `csolnp` from
#' each, and returns the solution with the lowest objective value. It is useful
#' for problems where local minima are a concern or the objective surface is
#' challenging.
#'
#' @param fn the objective function (must return a scalar).
#' @param gr an optional function for computing the analytic gradient of the
#' function (must return a vector of length n).
#' @param eq_fn an optional function for calculating equality constraints.
#' @param eq_b a vector of the equality bounds (if eq_fn provided).
#' @param eq_jac an optional function for computing the analytic Jacobian of the
#' equality function (a matrix with number of columns n and number of rows
#' equal to the number of equalities).
#' @param ineq_fn an optional function for calculating inequality constraints.
#' @param ineq_lower the lower bounds for the inequality constraints (must be finite).
#' @param ineq_upper the upper bounds for the inequality constraints (must be finite).
#' @param ineq_jac an optional function for computing the analytic Jacobian of the
#' inequality function (a matrix with number of columns n and number of rows equal
#' to the number of inequalities).
#' @param lower lower bounds for the parameters. This is strictly required.
#' @param upper upper bounds for the parameters. This is strictly required.
#' @param control a list of solver control parameters (see details).
#' @param n_candidates integer. The number of initial feasible candidate points
#' to generate for multi-start optimization. Default is 20.
#' @param penalty numeric. The penalty parameter used when projecting to feasibility
#' for candidate generation. Default is 1e4.
#' @param eq_tol Numeric. Tolerance for equality constraint violation
#' (default is 1e-6). Candidate solutions with \code{kkt_diagnostics\$eq_violation}
#' less than or equal to this value are considered feasible with respect
#' to equality constraints.
#' @param ineq_tol Numeric. Tolerance for inequality constraint violation
#' (default is 1e-6). Candidate solutions with \code{kkt_diagnostics\$ineq_violation}
#' less than or equal to this value are considered feasible with respect to
#' inequality constraints.
#' @param seed an optional random seed used to initialize the random number
#' generator for the random samples.
#' @param return_all logical. Whether to return all solutions as a list. This
#' may be useful for debugging.
#' @param ... additional arguments passed to the supplied functions
#' (common to all functions supplied).
#' @details
#' \strong{Candidate Generation:}
#'
#' The \code{generate_feasible_starts} approach creates a diverse set of initial
#' parameter vectors (candidates) that are feasible with respect to box and
#' (optionally) nonlinear constraints. The process is as follows:
#'
#' \enumerate{
#'   \item For each candidate, a random point is sampled inside the parameter
#'   box constraints (\code{lower} and \code{upper}) but a small distance away
#'   from the boundaries, to avoid numerical issues. This is achieved by a helper
#'   function that applies a user-specified buffer (\code{eps}).
#'   \item If nonlinear inequality constraints (\code{ineq_fn}) are provided,
#'   each sampled point is projected towards the feasible region using a
#'   fast penalized minimization. This step does not solve the feasibility problem
#'   exactly, but quickly produces a point that satisfies the constraints to
#'   within a specified tolerance, making it suitable as a starting point for
#'   optimization.
#'   \item If only box constraints are present, the sampled point is used
#'   directly as a feasible candidate.
#'   \item The set of feasible candidates is ranked by the objective with lower
#'   values considered better. This allows prioritization of candidates that start
#'   closer to optimality.
#' }
#'
#' This method efficiently creates a diverse set of robust initial values, improving
#' the chances that multi-start optimization will identify the global or a
#' high-quality local solution, especially in the presence of non-convexities or
#' challenging constraint boundaries.

#' \strong{Solution Selection:}
#' For each candidate starting point, \code{csolnp_ms} runs the \code{csolnp} solver
#' and collects the resulting solutions and their associated KKT diagnostics. If
#' equality or inequality constraints are present, candidate solutions are first
#' filtered to retain only those for which the maximum violation of
#' equality (\code{kkt_diagnostics\$eq_violation}) and/or
#' inequality (\code{kkt_diagnostics\$ineq_violation}) constraints are less than
#' or equal to user-specified tolerances (\code{eq_tol} and \code{ineq_tol}).
#' Among the feasible solutions (those satisfying all constraints within tolerance),
#' the solution with the lowest objective value is selected and returned as the
#' best result. If no candidate fully satisfies the constraints, the solution with
#' the smallest total constraint violation is returned, with a warning issued to
#' indicate that strict feasibility was not achieved. This two-stage selection
#' process ensures that the final result is both feasible (when possible) and
#' optimally minimizes the objective function among all feasible candidates.
#' @return
#' A list containing the best solution found by multi-start, with elements analogous
#' to those returned by \code{csolnp}. If \code{return_all} is TRUE, then a list of
#' all solutions is returned instead.
#' @seealso \code{\link{csolnp}}
#' @rdname csolnp_ms
#' @author Alexios Galanos
#' @export
csolnp_ms <- function(fn, gr = NULL, eq_fn = NULL, eq_b = NULL, eq_jac = NULL,
                      ineq_fn = NULL, ineq_lower = NULL, ineq_upper = NULL, ineq_jac = NULL,
                      lower = NULL, upper = NULL, control = list(), n_candidates = 20,
                      penalty = 1e4, eq_tol = 1e-6, ineq_tol = 1e-6, seed = NULL,
                      return_all = FALSE, ...)
{
    init <- generate_feasible_starts(n_candidates, fn, lower, upper, ineq_fn, ineq_lower, ineq_upper, eq_fn, eq_b,
                                     maxit = 100, penalty = penalty, eps = 1e-4, seed = seed)
    n <- NROW(init)
    sol <- future_lapply(1:n, function(i) {
        out <- try(csolnp(init[i,], fn = fn, gr = gr, eq_fn = eq_fn, eq_b = eq_b, eq_jac = eq_jac,
                          ineq_fn = ineq_fn, ineq_lower = ineq_lower, ineq_upper = ineq_upper, ineq_jac = ineq_jac,
                          lower = lower, upper = upper, control = control, ...), silent = TRUE)
        if (inherits(out, 'try-error')) return(NULL)
        return(out)
    }, future.seed = TRUE, future.packages = "Rsolnp")
    sol <- eval(sol)
    if (return_all) {
        return(sol)
    } else {
        # Remove NULL solutions up front
        sol <- sol[!sapply(sol, is.null)]
        if (length(sol) == 0) stop("No successful csolnp runs.")

        # Convergence filter as before
        convergence <- sapply(sol, function(x) if (is.null(x)) 10 else as.numeric(x$convergence))
        if (!any(convergence < 2)) {
            warning("No solution achieved convergence < 2. Results may be unreliable.")
        } else {
            sol <- sol[which(convergence < 2)]
        }
        # Get constraint violations
        has_eq <- !is.null(eq_fn)
        has_ineq <- !is.null(ineq_fn)

        eq_violation <- sapply(sol, function(x)
            if (!is.null(x$kkt_diagnostics$eq_violation)) as.numeric(x$kkt_diagnostics$eq_violation) else NA_real_)
        ineq_violation <- sapply(sol, function(x)
            if (!is.null(x$kkt_diagnostics$ineq_violation)) as.numeric(x$kkt_diagnostics$ineq_violation) else NA_real_)

        # Identify feasible candidates
        if (has_eq && has_ineq) {
            is_feasible <- (eq_violation <= eq_tol) & (ineq_violation <= ineq_tol)
        } else if (has_eq) {
            is_feasible <- (eq_violation <= eq_tol)
        } else if (has_ineq) {
            is_feasible <- (ineq_violation <= ineq_tol)
        } else {
            is_feasible <- rep(TRUE, length(sol)) # No constraints, all feasible
        }

        # Among feasible, choose lowest objective
        obj <- sapply(sol, function(x) if (is.null(x)) Inf else as.numeric(x$objective))

        if (any(is_feasible, na.rm=TRUE)) {
            idx <- which(is_feasible)
            min_index <- idx[which.min(obj[idx])]
            sol_best <- sol[[min_index]]
        } else {
            # No feasible solution, pick one with smallest total violation
            warning("No candidate fully satisfies constraints within tolerance. Returning least infeasible solution.")
            total_violation <- eq_violation
            if (has_eq && has_ineq) total_violation <- eq_violation + ineq_violation
            if (!has_eq && has_ineq) total_violation <- ineq_violation
            if (has_eq && !has_ineq) total_violation <- eq_violation
            min_index <- which.min(total_violation)
            sol_best <- sol[[min_index]]
        }
        return(sol_best)
    }
}
