#' Broyden-Fletcher-Goldfarb-Shanno (BFGS) Optimization
#'
#' @description
#' Implements the damped BFGS Quasi-Newton algorithm with a Strong Wolfe line search 
#' for non-linear optimization, specifically tailored for SEM.
#'
#' @details
#' \code{bfgs} is a Quasi-Newton method that maintains an approximation of the 
#' inverse Hessian matrix. It is widely considered the most robust and 
#' efficient member of the Broyden family of optimization methods.
#' 
#' \bold{BFGS vs. DFP:}
#' While both \code{bfgs} and \code{dfp} update the inverse Hessian using 
#' rank-two formulas, BFGS is generally more tolerant of inaccuracies in the 
#' line search. This implementation uses the Sherman-Morrison formula to 
#' update the inverse Hessian directly, avoiding the need for matrix inversion 
#' at each step.
#'
#' \bold{Strong Wolfe Line Search:}
#' To maintain the positive definiteness of the Hessian approximation and 
#' ensure global convergence, this algorithm employs a Strong Wolfe line search. 
#' This search identifies a step length \eqn{\alpha} that satisfies both sufficient 
#' decrease (Armijo condition) and the curvature condition.
#'
#' \bold{Damping for Non-Convexity:}
#' In Structural Equation Modeling (SEM), objective functions often exhibit 
#' non-convex regions. When \code{use_damped = TRUE}, Powell's damping 
#' strategy is applied to the update vectors to preserve the positive 
#' definiteness of the Hessian approximation even when the curvature condition 
#' is not naturally met.
#'
#' @references
#' \itemize{
#'    \item Nocedal, J., & Wright, S. J. (2006). \emph{Numerical Optimization}. Springer.
#'    \item Fletcher, R. (1987). \emph{Practical Methods of Optimization}. Wiley.
#' }
#'
#' @param start Numeric vector. Starting values for the optimization parameters.
#' @param objective Function. The objective function to minimize.
#' @param gradient Function (optional). Gradient of the objective function.
#' @param hessian Function (optional). Hessian matrix of the objective function.
#' @param lower Numeric vector. Lower bounds for box constraints.
#' @param upper Numeric vector. Upper bounds for box constraints.
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
#' 
#' @examples
#' # Simple quadratic function optimization
#' quad <- function(x) (x[1] - 2)^2 + (x[2] + 1)^2
#' res <- bfgs(start = c(0, 0), objective = quad)
#' print(res$par)
bfgs <- function(
    start,
    objective,
    gradient       = NULL,
    hessian        = NULL,
    lower          = -Inf,
    upper          = Inf,
    control        = list(),
    ...
) {
  
  # ---------- 1. Configuration (Synced with Optimization Suite) ----------
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
    max_iter        = 10000L,
    tol_abs_f       = 1e-6,
    tol_rel_f       = 1e-6,
    tol_abs_x       = 1e-6,
    tol_rel_x       = 1e-6,
    tol_grad        = 1e-4,
    tol_pred_f      = 1e-4,
    tol_pred_f_avg  = 1e-4,
    wolfe_c1        = 1e-4,
    wolfe_c2        = 0.9,
    ls_alpha0       = 1.0,
    ls_max_steps    = 30L,      
    zoom_max_steps  = 25L,
    curvature_eps   = 1e-12,
    Hinv_init_diag  = 1.0,
    diff_method      = "forward", 
    use_damped      = TRUE,      
    damp_phi        = 0.2        
  )
  ctrl <- utils::modifyList(ctrl0, control)
  ctrl$diff_method <- match.arg(ctrl$diff_method, c("forward", "central", "richardson"))
  
  if (ctrl$diff_method == "richardson") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' required.")
  }
  
  # ---------- 2. Internal Helpers ----------
  eval_obj <- function(z) as.numeric(objective(z, ...))[1]
  
  # Shared convergence test for the main loop (5.2) and the line-search-failure branch.
  # Returns whether the enabled gradient / parameter / function-value criteria are met;
  # the predicted-decrease tests need the accepted step, so they stay in 5.6b. The PD check
  # and status assignment stay at the call sites. Reading the flags here means changing the
  # convergence criteria affects both call sites identically.
  check_convergence <- function(g_inf, f, f_old, x, x_old, cur_it) {
    res_conv <- TRUE
    if (ctrl$use_grad) res_conv <- res_conv && (g_inf <= ctrl$tol_grad)
    if (ctrl$use_abs_f && !is.na(f_old)) res_conv <- res_conv && (abs(f - f_old) <= ctrl$tol_abs_f)
    if (ctrl$use_rel_f && !is.na(f_old)) res_conv <- res_conv && (abs((f - f_old) / max(1, abs(f_old))) <= ctrl$tol_rel_f)
    if (ctrl$use_abs_x && cur_it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) <= ctrl$tol_abs_x)
    if (ctrl$use_rel_x && cur_it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) / max(1, max(abs(x_old)))) <= ctrl$tol_rel_x
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
  
  # ---------- 3. Strong Wolfe Line Search ----------
  strong_wolfe_ls <- function(x, f, g, p, dphi0) {
    c1 <- ctrl$wolfe_c1; c2 <- ctrl$wolfe_c2; phi0 <- f
    
    zoom <- function(alo, ahi, flo, glo) {
      for (z_it in seq_len(ctrl$zoom_max_steps)) {
        a <- 0.5 * (alo + ahi); xa <- x + a * p; fa <- eval_obj(xa)
        if (fa > phi0 + c1 * a * dphi0 || fa >= flo) {
          ahi <- a
        } else {
          ga <- grad_func(xa); dphi_a <- sum(ga * p)
          if (abs(dphi_a) <= -c2 * dphi0) return(list(ok=TRUE, alpha=a, x=xa, f=fa, g=ga, status="wolfe"))
          if (dphi_a * (ahi - alo) >= 0) ahi <- alo
          alo <- a; flo <- fa; glo <- ga
        }
        if (abs(ahi - alo) < 1e-15) break 
      }
      list(ok=FALSE, status="zoom_failed")
    }
    
    a_prev <- 0.0; f_prev <- phi0; g_prev <- g; a <- ctrl$ls_alpha0
    for (ls_it in seq_len(ctrl$ls_max_steps)) {
      xa <- x + a * p; fa <- eval_obj(xa)
      if (fa > phi0 + c1 * a * dphi0 || (ls_it > 1L && fa >= f_prev)) return(zoom(a_prev, a, f_prev, g_prev))
      ga <- grad_func(xa); dphi_a <- sum(ga * p)
      if (abs(dphi_a) <= -c2 * dphi0) return(list(ok=TRUE, alpha=a, x=xa, f=fa, g=ga, status="wolfe"))
      if (dphi_a >= 0) return(zoom(a, a_prev, fa, ga))
      a_prev <- a; f_prev <- fa; g_prev <- ga; a <- a * 2.0
    }
    list(ok=FALSE, status="line_search_failed")
  }
  
  # ---------- 4. Optimization Initialization ----------
  x <- as.numeric(start); n <- length(x)
  
  # CRITICAL: start_clock must be defined here for final cpu_time reporting
  start_clock <- proc.time()
  
  f <- tryCatch(eval_obj(x), error = function(e) NA_real_)
  it <- 0L; x_old <- x; f_old <- NA_real_; converged <- FALSE; status <- "running"
  Hinv <- diag(ctrl$Hinv_init_diag, n); g_inf <- NA_real_; H_eval <- NULL
  pred_dec <- NA_real_; pred_dec_avg <- NA_real_
  
  if (!is.finite(f)) { 
    status <- "objective_error_at_start" 
  } else {
    g <- grad_func(x)
    
    # ---------- 5. Main Loop ----------
    tryCatch({
      repeat {
        if (it >= ctrl$max_iter) { status <- "iteration_limit_reached"; break }
        it <- it + 1L; g_inf <- max(abs(g), na.rm = TRUE)
        
        # 5.1) Compute Search Direction & Reset if not descent
        p <- as.numeric(-Hinv %*% g); dphi0 <- sum(g * p)
        if (dphi0 >= 0) {
          Hinv <- diag(ctrl$Hinv_init_diag, n); p <- -g; dphi0 <- sum(g * p)
        }
        
        # 5.2) Convergence Verification
        # Gradient / parameter / function-value tests use the current point (before the line
        # search), via the shared check_convergence so all flags are honored. The predicted-
        # decrease tests need the step the line search selects, so they are evaluated after
        # the line search (see 5.6b) and AND-combined with res_conv_pre.
        res_conv_pre <- check_convergence(g_inf, f, f_old, x, x_old, it)
        
        # When neither predicted-decrease test is active, the convergence decision is
        # complete here, so finalize immediately (this also skips an unnecessary line
        # search when already converged).
        if (!isTRUE(ctrl$use_pred_f) && !isTRUE(ctrl$use_pred_f_avg)) {
          if (res_conv_pre && it > 1) {
            if (isTRUE(ctrl$use_posdef)) {
              H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
              if (is_pd_fast(H_eval)) { status <- "converged"; converged <- TRUE; break } else { res_conv_pre <- FALSE }
            } else { status <- "converged"; converged <- TRUE; break }
          }
        }
        
        # 5.3) Line Search
        ls <- strong_wolfe_ls(x, f, g, p, dphi0)
        if (!isTRUE(ls$ok)) {
          # Line search found no acceptable step. Re-check convergence with the same criteria
          # as 5.2 (via check_convergence) so user-configured stopping rules are honored: a
          # line search can fail simply because we are at a stationary point. A non-positive-
          # definite point here is a genuine failure rather than convergence.
          if (check_convergence(g_inf, f, f_old, x, x_old, it)) {
            if (isTRUE(ctrl$use_posdef)) {
              H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
              if (is_pd_fast(H_eval)) { status <- "converged"; converged <- TRUE; break }
            } else { status <- "converged"; converged <- TRUE; break }
          }
          status <- ls$status; break
        }
        
        # 5.4) Update Parameters
        alpha_final <- ls$alpha; x_new <- ls$x; f_new <- ls$f; g_new <- ls$g
        s <- x_new - x; y <- g_new - g; sy <- sum(s * y)
        
        # 5.4a) Post-line-search convergence check (handles exact solutions, e.g., quadratics)
        g_inf_new <- max(abs(g_new), na.rm = TRUE)
        if (ctrl$use_grad && g_inf_new <= ctrl$tol_grad) {
          x <- x_new; f <- f_new; g <- g_new; g_inf <- g_inf_new
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
            if (is_pd_fast(H_eval)) { status <- "converged"; converged <- TRUE; break }
          } else { status <- "converged"; converged <- TRUE; break }
        }
        update_ok <- FALSE; y_star <- y; sy_star <- sy
        B_s_approx <- -alpha_final * g; sBs <- as.numeric(sum(s * B_s_approx))
        
        # 5.5) Powell's Damping Strategy
        if (isTRUE(ctrl$use_damped)) {
          if (is.finite(sBs) && sBs > ctrl$curvature_eps) {
            if (sy < ctrl$damp_phi * sBs) {
              theta <- ((1 - ctrl$damp_phi) * sBs) / (sBs - sy)
              y_star <- theta * y + (1 - theta) * B_s_approx; sy_star <- sum(s * y_star)
            }
          }
          if (is.finite(sy_star) && sy_star > ctrl$curvature_eps) { y <- y_star; sy <- sy_star; update_ok <- TRUE }
        } else { if (is.finite(sy) && sy > ctrl$curvature_eps) update_ok <- TRUE }
        
        # 5.6) Predicted Decrease
        if (isTRUE(ctrl$use_pred_f) || isTRUE(ctrl$use_pred_f_avg)) {
          pred_dec <- as.numeric(-(sum(g * s) + 0.5 * sBs)); pred_dec_avg <- pred_dec / n
        }
        
        # 5.6b) Predicted-decrease convergence tests (use the step selected by the line
        # search, i.e. the same pred_dec that is reported). AND-combined with the
        # gradient/parameter/function-value result from 5.2 (res_conv_pre).
        if (isTRUE(ctrl$use_pred_f) || isTRUE(ctrl$use_pred_f_avg)) {
          res_conv <- res_conv_pre
          if (isTRUE(ctrl$use_pred_f))     res_conv <- res_conv && (is.finite(pred_dec) && pred_dec <= ctrl$tol_pred_f)
          if (isTRUE(ctrl$use_pred_f_avg)) res_conv <- res_conv && (is.finite(pred_dec_avg) && pred_dec_avg <= ctrl$tol_pred_f_avg)
          if (res_conv && it > 1) {
            if (isTRUE(ctrl$use_posdef)) {
              H_eval <- tryCatch(hess_func(x_new), error = function(e) NULL)
              if (is_pd_fast(H_eval)) {
                x_old <- x; f_old <- f; x <- x_new; f <- f_new; g <- g_new
                status <- "converged"; converged <- TRUE; break
              }
            } else {
              x_old <- x; f_old <- f; x <- x_new; f <- f_new; g <- g_new
              status <- "converged"; converged <- TRUE; break
            }
          }
        }
        
        # 5.7) Inverse Hessian Update (Sherman-Morrison Formula)
        if (update_ok) {
          # Added a small epsilon (1e-16) to sy to prevent division by zero or Inf
          rho <- 1 / (sy + 1e-16) 
          I_mat <- diag(n)
          
          # Initial scaling: Adjust the initial Hinv to better match the true Hessian scale
          Hy <- as.numeric(Hinv %*% y)
          yHy <- as.numeric(crossprod(y, Hy))
          if (it == 1L && is.finite(yHy) && yHy > 1e-12) {
            Hinv <- Hinv * (sy / yHy)
          }
          
          # Sherman-Morrison Update: Update Hinv directly to avoid matrix inversion
          # Use double-sided multiplication for better numerical stability
          H_left  <- (I_mat - rho * (s %*% t(y)))
          H_right <- (I_mat - rho * (y %*% t(s)))
          Hinv <- H_left %*% Hinv %*% H_right + rho * (s %*% t(s))
          
          # Force symmetry: Prevent accumulation of rounding errors over many iterations
          Hinv <- 0.5 * (Hinv + t(Hinv))
        }
        
        x_old <- x; f_old <- f; x <- x_new; f <- f_new; g <- g_new
      }
    }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
  }
  
  # ---------- 6. Final Status & Output Construction ----------
  if (is.null(H_eval)) H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
  H_final <- if (!is.null(H_eval)) H_eval else tryCatch(solve(Hinv), error = function(e) matrix(NA, n, n))
  
  final_clock <- proc.time() - start_clock
  list(
    par          = x, 
    objective    = f, 
    converged    = converged, 
    status       = status, 
    iter         = it, 
    cpu_time     = as.numeric(final_clock[1] + final_clock[2]), 
    elapsed_time = as.numeric(final_clock[3]),
    max_grad     = g_inf, 
    Hess_is_pd   = is_pd_fast(H_final), 
    Hessian      = H_final, 
    approx_hinv  = Hinv, 
    pred_dec     = pred_dec, 
    pred_dec_avg = pred_dec_avg
  )
}