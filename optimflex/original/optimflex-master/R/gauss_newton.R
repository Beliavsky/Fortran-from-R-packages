#' Gauss-Newton Optimization
#'
#' @description
#' Implements a full-featured Gauss-Newton algorithm for non-linear optimization, 
#' specifically optimized for Structural Equation Modeling (SEM).
#'
#' @details
#' \code{gauss_newton} is a specialized optimization algorithm for least-squares 
#' and Maximum Likelihood problems where the objective function can be 
#' expressed as a sum of squared residuals.
#' 
#' \bold{Scaling and SEM Consistency:}
#' To ensure consistent simulation results and standard error (SE) calculations, 
#' this implementation adjusts the Gradient \eqn{(2J^T r)} and the Approximate 
#' Hessian \eqn{(2J^T J)} to match the scale of the Maximum Likelihood (ML) 
#' fitting function \eqn{F_{ML}}. This alignment is critical when calculating 
#' asymptotic covariance matrices using the formula \eqn{\frac{2}{n} H^{-1}}.
#'
#' \bold{Comparison with Newton-Raphson:}
#' Unlike \code{newton_raphson} or \code{modified_newton}, which require the full 
#' second-order Hessian, Gauss-Newton approximates the Hessian using the 
#' Jacobian of the residuals. This is computationally more efficient and 
#' provides a naturally positive-semidefinite approximation, though a ridge 
#' adjustment is still provided for numerical stability.
#'
#' \bold{Ridge Adjustment Strategy:}
#' The function includes a "Ridge Rescue" mechanism. If the approximate Hessian 
#' is singular or poorly conditioned for Cholesky decomposition, it iteratively 
#' adds a diagonal ridge \eqn{(\tau I)} until numerical stability is achieved.
#'
#' @references
#' \itemize{
#'    \item Nocedal, J., & Wright, S. J. (2006). \emph{Numerical Optimization}. Springer.
#'    \item Bollen, K. A. (1989). \emph{Structural Equations with Latent Variables}. Wiley.
#' }
#'
#' @param start Numeric vector. Starting values for the optimization parameters.
#' @param objective Function. The objective function to minimize.
#' @param residual Function (optional). Function that returns the residuals vector.
#' @param gradient Function (optional). Gradient of the objective function.
#' @param hessian Function (optional). Hessian matrix of the objective function.
#' @param jac Function (optional). Jacobian matrix of the residuals.
#' @param control List. Control parameters including convergence flags:
#'    \itemize{
#'      \item \code{use_abs_f}: Logical. Use absolute change in objective for convergence.
#'      \item \code{use_rel_f}: Logical. Use relative change in objective for convergence.
#'      \item \code{use_abs_x}: Logical. Use absolute change in parameters for convergence.
#'      \item \code{use_rel_x}: Logical. Use relative change in parameters for convergence.
#'      \item \code{use_grad}: Logical. Use gradient norm for convergence.
#'      \item \code{use_posdef}: Logical. Verify positive definiteness at convergence.
#'      \item \code{use_pred_f}: Logical. Record predicted objective decrease.
#'      \item \code{use_pred_f_avg}: Logical. Record average predicted decrease.
#'      \item \code{diff_method}: String. Method for numerical differentiation.
#'    }
#' @param ... Additional arguments passed to objective, gradient, and Hessian functions.
#'
#' @return A list containing optimization results and iteration metadata.
#' @export
#' @examples
#' # Simple quadratic function optimization
#' quad <- function(x) (x[1] - 2)^2 + (x[2] + 1)^2
#' res <- gauss_newton(start = c(0, 0), objective = quad)
#' print(res$par)
gauss_newton <- function(
    start,
    objective,
    residual       = NULL,
    gradient       = NULL,
    hessian        = NULL,
    jac            = NULL,
    control        = list(),
    ...
) {
  
  # ---------- 1. Default Configuration (Synced with Suite) ----------
  ctrl0 <- list(
    # Convergence and recording flags
    use_abs_f       = FALSE,
    use_rel_f       = FALSE,
    use_abs_x       = FALSE,
    use_rel_x       = TRUE,
    use_grad        = TRUE,
    use_posdef      = TRUE,
    use_pred_f      = FALSE,
    use_pred_f_avg  = FALSE,
    
    # Algorithm parameters
    max_iter        = 1000L,
    tol_abs_f       = 1e-6,
    tol_rel_f       = 1e-6,
    tol_abs_x       = 1e-6,
    tol_rel_x       = 1e-6,
    tol_grad        = 1e-4,
    tol_pred_f      = 1e-4,
    tol_pred_f_avg  = 1e-4,
    wolfe_c1        = 1e-4,
    ls_alpha0       = 1.0,
    ls_max_steps    = 30L,
    ridge_offset    = 1e-4,
    ridge_mult      = 10.0,   
    ridge_max_tries = 8L,     
    diff_method      = "forward" 
  )
  ctrl <- utils::modifyList(ctrl0, control)
  ctrl$diff_method <- match.arg(ctrl$diff_method, c("forward", "central", "richardson"))
  
  if (ctrl$diff_method == "richardson") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' required.")
  }
  
  # ---------- 2. Internal Helpers ----------
  eval_obj <- function(z) as.numeric(objective(z, ...))[1]
  
  # Shared convergence test for the main loop (4.4) and the line-search-failure branch.
  # Returns only whether the enabled stopping criteria are met; the positive-definiteness
  # check and status assignment stay at the call sites. Reading every flag from ctrl here
  # means changing the convergence criteria affects both call sites identically.
  check_convergence <- function(g_inf, f, f_old, x, x_old, pred_dec, pred_dec_avg, cur_it) {
    res_conv <- TRUE
    if (ctrl$use_grad)                   res_conv <- res_conv && (g_inf <= ctrl$tol_grad)
    if (ctrl$use_abs_f && !is.na(f_old)) res_conv <- res_conv && (abs(f - f_old) <= ctrl$tol_abs_f)
    if (ctrl$use_rel_f && !is.na(f_old)) res_conv <- res_conv && (abs((f - f_old) / max(1, abs(f_old))) <= ctrl$tol_rel_f)
    if (ctrl$use_abs_x && cur_it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) <= ctrl$tol_abs_x)
    if (ctrl$use_rel_x && cur_it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) / max(1, max(abs(x_old)))) <= ctrl$tol_rel_x
    if (isTRUE(ctrl$use_pred_f))      res_conv <- res_conv && (is.finite(pred_dec) && pred_dec <= ctrl$tol_pred_f)
    if (isTRUE(ctrl$use_pred_f_avg))  res_conv <- res_conv && (is.finite(pred_dec_avg) && pred_dec_avg <= ctrl$tol_pred_f_avg)
    res_conv
  }

  grad_func <- if (!is.null(gradient)) {
    function(z) as.numeric(gradient(z, ...))
  } else if (ctrl$diff_method == "richardson") {
    function(z) as.numeric(numDeriv::grad(objective, z, method = "Richardson", ...))
  } else {
    function(z) fast_grad(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  hess_func <- if (!is.null(hessian)) {
    function(z) hessian(z, ...)
  } else if (ctrl$diff_method == "richardson") {
    function(z) numDeriv::hessian(objective, z, method = "Richardson", ...)
  } else {
    function(z) fast_hess(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  jac_func <- if (!is.null(jac)) {
    function(z) jac(z, ...)
  } else if (!is.null(residual)) {
    if (ctrl$diff_method == "richardson") {
      function(z) numDeriv::jacobian(residual, z, method = "Richardson", ...)
    } else {
      function(z) fast_jac(residual, z, diff_method = ctrl$diff_method, ...)
    }
  } else {
    NULL
  }
  
  solve_with_ridge <- function(H, rhs) {
    H <- 0.5 * (H + t(H)) 
    L <- try(chol(H), silent = TRUE)
    if (!inherits(L, "try-error")) {
      return(list(step = backsolve(L, backsolve(L, rhs, transpose = TRUE)), H_used = H, tau = 0, chol_ok = TRUE))
    }
    tau <- max(ctrl$ridge_offset, 1e-12)
    H_mod <- H
    for (j in seq_len(ctrl$ridge_max_tries)) {
      diag(H_mod) <- diag(H) + tau
      L <- try(chol(H_mod), silent = TRUE)
      if (!inherits(L, "try-error")) {
        return(list(step = backsolve(L, backsolve(L, rhs, transpose = TRUE)), H_used = H_mod, tau = tau, chol_ok = TRUE))
      }
      tau <- tau * ctrl$ridge_mult
    }
    diag(H_mod) <- diag(H) + tau
    list(step = tryCatch(solve(H_mod, rhs), error = function(e) rep(NA, length(rhs))), H_used = H_mod, tau = tau, chol_ok = FALSE)
  }
  
  # ---------- 3. Initialization ----------
  x <- as.numeric(start); n_par <- length(x)
  
  # start_clock defined here for CPU time calculation
  start_clock <- proc.time()
  
  f <- tryCatch(eval_obj(x), error = function(e) NA_real_)
  it <- 0L; x_old <- x; f_old <- NA_real_; converged <- FALSE; status <- "running"
  H_curr <- NULL; H_eval <- NULL; g_inf <- NA_real_; pred_dec <- NA_real_; pred_dec_avg <- NA_real_
  
  if (!is.finite(f)) {
    status <- "objective_error_at_start"
  } else if (is.null(jac_func)) {
    status <- "jacobian_unavailable"
  } else {
    # Scale Gradient (2*J'r) to match F_ML first-order derivatives
    get_g <- function(curr_x, curr_J) {
      if (!is.null(gradient)) return(grad_func(curr_x))
      if (!is.null(residual)) return(as.numeric(2 * crossprod(curr_J, as.numeric(residual(curr_x, ...)))))
      return(grad_func(curr_x))
    }
    
    J <- jac_func(x)
    g <- get_g(x, J)
    
    # ---------- 4. Main Optimization Loop ----------
    tryCatch({
      repeat {
        if (it >= ctrl$max_iter) { status <- "iteration_limit_reached"; break }
        it <- it + 1L; g_inf <- max(abs(g), na.rm = TRUE)
        
        # 4.1) Jacobian and Approximate Hessian (2 * J'J)
        if (is.null(J) || any(!is.finite(J))) { status <- "jacobian_error"; break }
        
        # Calculate approx_hessian as 2*J'J to match Newton-Raphson scale
        H_curr <- 2 * crossprod(J)
        
        # 4.2) Compute Step with Ridge Rescue
        sol <- solve_with_ridge(H_curr, rhs = -g)
        step <- as.numeric(sol$step); H_used <- sol$H_used
        if (any(!is.finite(step))) { status <- "step_failed"; break }
        
        # 4.3) Predicted Decrease
        if (isTRUE(ctrl$use_pred_f) || isTRUE(ctrl$use_pred_f_avg)) {
          Hp <- as.numeric(H_used %*% step)
          pred_dec <- as.numeric(-(sum(g * step) + 0.5 * sum(step * Hp)))
          pred_dec_avg <- pred_dec / n_par
        }
        
        # 4.4) Convergence Verification (shared check_convergence; all ctrl flags honored)
        res_conv <- check_convergence(g_inf, f, f_old, x, x_old, pred_dec, pred_dec_avg, it)
        
        if (res_conv && it > 1L) {
          converged <- TRUE; status <- "converged"; break
        }
        
        # 4.5) Backtracking Line Search (Armijo condition)
        dphi0 <- sum(g * step); alpha <- ctrl$ls_alpha0; ls_ok <- FALSE
        for (ls_it in seq_len(ctrl$ls_max_steps)) {
          xi <- x + alpha * step; fi <- eval_obj(xi)
          if (is.finite(fi) && fi <= f + ctrl$wolfe_c1 * alpha * dphi0) {
            x_old <- x; f_old <- f; x <- xi; f <- fi; ls_ok <- TRUE; break
          }
          alpha <- alpha * 0.5
        }
        if (!ls_ok) {
          # Line search found no acceptable step. Re-check convergence with the same
          # criteria as 4.4 (via check_convergence) so user-configured stopping rules are
          # honored: a line search can fail simply because we are at a stationary point.
          # The finalization PD check still guards against non-minimizing points.
          if (check_convergence(g_inf, f, f_old, x, x_old, pred_dec, pred_dec_avg, it)) {
            converged <- TRUE; status <- "converged"; break
          }
          status <- "line_search_failed"; break
        }
        
        # Prepare for next iteration
        J <- jac_func(x)
        g <- get_g(x, J)
        
        # Post-line-search convergence check (handles exact solutions, e.g., quadratics)
        g_inf_new <- max(abs(g), na.rm = TRUE)
        if (ctrl$use_grad && g_inf_new <= ctrl$tol_grad) {
          g_inf <- g_inf_new
          converged <- TRUE; status <- "converged"; break
        }
      }
    }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
  }
  
  # ---------- 5. Finalization & Mandatory PD Check ----------
  if (converged) {
    H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
    Hess_pd <- if (!is.null(H_eval)) is_pd_fast(H_eval) else FALSE
    
    if (isTRUE(ctrl$use_posdef) && !Hess_pd) {
      converged <- FALSE
      status <- "converged_but_not_positive_definite"
    }
  } else {
    Hess_pd <- FALSE
  }
  
  H_final <- if (!is.null(H_eval)) H_eval else if (!is.null(H_curr)) H_curr else NA_real_
  final_clock <- proc.time() - start_clock
  
  list(par = x, objective = f, converged = converged, status = status, iter = it,
       cpu_time = as.numeric(final_clock[1] + final_clock[2]), elapsed_time = as.numeric(final_clock[3]),
       max_grad = as.numeric(g_inf), Hess_is_pd = Hess_pd, Hessian = H_final, approx_hessian = H_curr,
       pred_dec = pred_dec, pred_dec_avg = pred_dec_avg)
}