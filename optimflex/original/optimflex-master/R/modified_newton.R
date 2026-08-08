#' Modified Newton-Raphson Optimization
#'
#' @description
#' Implements an optimized Newton-Raphson algorithm for non-linear optimization 
#' featuring dynamic ridge adjustment and backtracking line search.
#'
#' @details
#' \code{modified_newton} is a line search optimization algorithm that utilizes 
#' second-order curvature information (the Hessian matrix) to find the minimum 
#' of an objective function.
#' 
#' \bold{Modified Newton vs. Trust-Region:}
#' Unlike the \code{dogleg} and \code{double_dogleg} functions which use a 
#' Trust-Region approach to constrain the step size, this function uses a 
#' \bold{Line Search} approach. It first determines the Newton direction 
#' (the solution to \eqn{H \Delta x = -g}) and then performs a backtracking line 
#' search to find a step length \eqn{\alpha} that satisfies the sufficient decrease 
#' condition (Armijo condition).
#'
#' \bold{Dynamic Ridge Adjustment:}
#' If the Hessian matrix \eqn{H} is not positive definite (making it unsuitable for 
#' Cholesky decomposition), the algorithm applies a dynamic ridge adjustment. 
#' A diagonal matrix \eqn{\tau I} is added to the Hessian, where \eqn{\tau} is 
#' increased until the matrix becomes positive definite. This ensures the 
#' search direction always remains a descent direction.
#'
#' \bold{Differentiation Methods:}
#' The function allows for independent selection of differentiation methods for 
#' the gradient and Hessian:
#' \itemize{
#'    \item \code{forward}: Standard forward-difference numerical differentiation.
#'    \item \code{central}: Central-difference (more accurate but slower).
#'    \item \code{complex}: Complex-step differentiation (highly accurate for gradients).
#'    \item \code{richardson}: Richardson extrapolation via the \code{numDeriv} package.
#' }
#'
#' @param start Numeric vector. Starting values for the optimization parameters.
#' @param objective Function. The objective function to minimize.
#' @param gradient Function (optional). Gradient of the objective function.
#' @param hessian Function (optional). Hessian matrix of the objective function.
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
#'      \item \code{grad_diff}: String. Method for gradient differentiation.
#'      \item \code{hess_diff}: String. Method for Hessian differentiation.
#'    }
#' @param ... Additional arguments passed to objective, gradient, and Hessian functions.
#'
#' @return A list containing optimization results and iteration metadata.
#' @export
#' @examples
#' # Simple quadratic function optimization
#' quad <- function(x) (x[1] - 2)^2 + (x[2] + 1)^2
#' res <- modified_newton(start = c(0, 0), objective = quad)
#' print(res$par)
modified_newton <- function(
    start,
    objective,
    gradient       = NULL,
    hessian        = NULL,
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
    grad_diff       = "forward",
    hess_diff       = "forward"
  )
  ctrl <- utils::modifyList(ctrl0, control)
  
  ctrl$grad_diff <- match.arg(ctrl$grad_diff, c("forward", "central", "richardson", "complex"))
  ctrl$hess_diff <- match.arg(ctrl$hess_diff, c("forward", "central", "richardson"))
  
  if (any(c(ctrl$grad_diff, ctrl$hess_diff) %in% c("richardson", "complex"))) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' required.")
  }
  
  # ---------- 2. Internal Helpers ----------
  eval_obj <- function(z) as.numeric(objective(z, ...))[1]
  
  # Shared convergence test for the main loop (4.3) and the line-search-failure branch.
  # Returns only whether the enabled stopping criteria are met; posdef checks and status
  # assignment stay at the call sites. All criteria are read from ctrl.
  check_convergence <- function(g_inf, f, f_old, x, x_old, pred_dec, pred_dec_avg, cur_it) {
    res_conv <- TRUE
    if (ctrl$use_grad) res_conv <- res_conv && (g_inf <= ctrl$tol_grad)
    if (ctrl$use_abs_f && !is.na(f_old)) res_conv <- res_conv && (abs(f - f_old) <= ctrl$tol_abs_f)
    if (ctrl$use_rel_f && !is.na(f_old)) {
      res_conv <- res_conv && (abs((f - f_old) / max(1, abs(f_old))) <= ctrl$tol_rel_f)
    }
    if (ctrl$use_abs_x && cur_it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) <= ctrl$tol_abs_x)
    if (ctrl$use_rel_x && cur_it > 1L) {
      res_conv <- res_conv && (max(abs(x - x_old)) / max(1, max(abs(x_old))) <= ctrl$tol_rel_x)
    }
    if (isTRUE(ctrl$use_pred_f)) res_conv <- res_conv && (is.finite(pred_dec) && pred_dec <= ctrl$tol_pred_f)
    if (isTRUE(ctrl$use_pred_f_avg)) res_conv <- res_conv && (is.finite(pred_dec_avg) && pred_dec_avg <= ctrl$tol_pred_f_avg)
    res_conv
  }
  
  grad_func <- if (!is.null(gradient)) {
    function(z) as.numeric(gradient(z, ...))
  } else if (ctrl$grad_diff == "richardson") {
    function(z) as.numeric(numDeriv::grad(objective, z, method = "Richardson", ...))
  } else if (ctrl$grad_diff == "complex") {
    function(z) as.numeric(numDeriv::grad(objective, z, method = "complex", ...))
  } else {
    function(z) fast_grad(objective, z, diff_method = ctrl$grad_diff, ...)
  }
  
  hess_func <- if (!is.null(hessian)) {
    function(z) hessian(z, ...)
  } else if (ctrl$hess_diff == "richardson") {
    function(z) numDeriv::hessian(objective, z, method = "Richardson", ...)
  } else {
    function(z) fast_hess(objective, z, diff_method = ctrl$hess_diff, ...)
  }
  
  # ---------- 3. Initialization ----------
  x <- as.numeric(start); n <- length(x); start_clock <- proc.time()
  f <- tryCatch(eval_obj(x), error = function(e) NA_real_)
  it <- 0L; x_old <- x; f_old <- f; converged <- FALSE; status <- "running"
  H_last <- NULL; is_pd <- FALSE; g_inf <- NA_real_
  pred_dec <- NA_real_; pred_dec_avg <- NA_real_
  
  if (!is.finite(f)) { 
    status <- "objective_error_at_start" 
  } else {
    g <- grad_func(x)
    
    # ---------- 4. Main Optimization Loop ----------
    tryCatch({
      repeat {
        if (it >= ctrl$max_iter) { status <- "iteration_limit_reached"; break }
        it <- it + 1L; g_inf <- max(abs(g), na.rm = TRUE)
        
        # 4.1) Single Hessian Evaluation: Forced symmetry
        H <- hess_func(x); H <- 0.5 * (H + t(H)); H_last <- H
        
        # 4.2) Step Calculation with Dynamic Ridge Adjustment
        # H_used records the matrix actually factorized to solve for the step (H itself
        # when positive definite, otherwise the ridge-adjusted H_mod). It is used for the
        # predicted-decrease quadratic model below, consistent with gauss_newton().
        R <- try(chol(H), silent = TRUE)
        if (!inherits(R, "try-error")) {
          is_pd <- TRUE
          step <- backsolve(R, forwardsolve(t(R), -g)); H_used <- H
        } else {
          is_pd <- FALSE
          tau <- ctrl$ridge_offset
          repeat {
            H_mod <- H + diag(tau, n)
            R <- try(chol(H_mod), silent = TRUE)
            if (!inherits(R, "try-error")) {
              step <- backsolve(R, forwardsolve(t(R), -g)); H_used <- H_mod
              break
            }
            tau <- tau * 10
            if (tau > 1e6) { status <- "ridge_failed"; break }
          }
        }
        
        if (status == "ridge_failed") break
        gTp <- sum(g * step) 
        
        # 4.2b) Predicted Decrease of the Newton step (objective-scale quadratic model
        # -(g'p + 0.5 p'H p), using the matrix actually factorized, H_used). The step is
        # known before the line search, so all criteria are tested together in 4.3.
        if (isTRUE(ctrl$use_pred_f) || isTRUE(ctrl$use_pred_f_avg)) {
          pred_dec <- as.numeric(-(gTp + 0.5 * sum(step * (H_used %*% step)))); pred_dec_avg <- pred_dec / n
        }
        
        # 4.3) Convergence Verification (shared check_convergence, all ctrl flags honored)
        res_conv <- check_convergence(g_inf, f, f_old, x, x_old, pred_dec, pred_dec_avg, it)
        
        if (res_conv) {
          if (isTRUE(ctrl$use_posdef)) {
            if (is_pd) { converged <- TRUE; status <- "converged"; break } else { res_conv <- FALSE }
          } else { converged <- TRUE; status <- "converged"; break }
        }
        
        # 4.4) Backtracking Line Search (Armijo condition)
        alpha <- ctrl$ls_alpha0; ls_ok <- FALSE
        for (ls_it in seq_len(ctrl$ls_max_steps)) {
          xi <- x + alpha * step; fi <- eval_obj(xi)
          if (is.finite(fi) && fi <= f + ctrl$wolfe_c1 * alpha * gTp) {
            x_old <- x; f_old <- f; x <- xi; f <- fi; ls_ok <- TRUE; break
          }
          alpha <- alpha * 0.5 
        }
        if (!ls_ok) {
          # Re-check convergence at the current point before failing: a line search can
          # fail simply because g ~ 0 (nothing left to decrease). Same criteria as 4.3,
          # but a non-PD point here is a genuine failure rather than convergence.
          res_conv <- check_convergence(g_inf, f, f_old, x, x_old, pred_dec, pred_dec_avg, it)
          if (res_conv) {
            if (isTRUE(ctrl$use_posdef)) {
              if (is_pd) { converged <- TRUE; status <- "converged"; break }
            } else { converged <- TRUE; status <- "converged"; break }
          }
          status <- "line_search_failed"; break
        }
        
        g <- grad_func(x)
        
        # Post-line-search convergence check (handles exact solutions, e.g., quadratics)
        g_inf_new <- max(abs(g), na.rm = TRUE)
        if (ctrl$use_grad && g_inf_new <= ctrl$tol_grad) {
          g_inf <- g_inf_new
          if (isTRUE(ctrl$use_posdef)) {
            H_last <- tryCatch(hess_func(x), error = function(e) NULL)
            is_pd <- is_pd_fast(H_last)
            if (is_pd) { converged <- TRUE; status <- "converged"; break }
          } else { converged <- TRUE; status <- "converged"; break }
        }
      }
    }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
  }
  
  # ---------- 5. Final Reporting ----------
  final_clock <- proc.time() - start_clock
  list(
    par = x, objective = f, converged = converged, status = status, iter = it,
    cpu_time = as.numeric(final_clock[1] + final_clock[2]), 
    elapsed_time = as.numeric(final_clock[3]),
    max_grad = as.numeric(g_inf), Hess_is_pd = is_pd, Hessian = H_last,
    pred_dec = pred_dec, pred_dec_avg = pred_dec_avg
  )
}