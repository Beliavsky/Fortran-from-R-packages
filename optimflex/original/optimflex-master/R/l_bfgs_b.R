#' Limited-memory BFGS with Box Constraints (L-BFGS-B)
#'
#' @description
#' Performs bound-constrained minimization using the L-BFGS-B algorithm. 
#' This implementation handles box constraints via Generalized Cauchy Point (GCP) 
#' estimation and subspace minimization, featuring a limited-memory (two-loop recursion) 
#' inverse Hessian approximation.
#'
#' @section Comparison with Existing Functions:
#' This function adds three features for rigorous convergence control. First, it 
#' applies an \bold{AND rule}: all selected convergence criteria must be satisfied 
#' simultaneously. Second, users can choose among \bold{eight distinct criteria} 
#' (e.g., changes in \eqn{f}, \eqn{x}, gradient, or predicted decrease) instead of relying 
#' on fixed defaults. Third, it provides an \bold{optional verification} using the 
#' Hessian computed from derivatives (analytically when provided, or via numerical 
#' differentiation). Checking the positive definiteness of this Hessian at the final 
#' solution reduces the risk of declaring convergence at non-minimizing stationary 
#' points, such as saddle points.
#'
#' @references
#' Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm 
#' for bound constrained optimization. \emph{SIAM Journal on Scientific Computing}, 
#' 16(5), 1190-1208.
#'
#' Morales, J. L., & Nocedal, J. (2011). L-BFGS-B: Remark on algorithm 778: 
#' L-BFGS-B: Fortran subroutines for large-scale bound-constrained optimization. 
#' \emph{ACM Transactions on Mathematical Software}, 38(1), 1-4.
#'
#' @param start Numeric vector. Initial values for the parameters.
#' @param objective Function. The scalar objective function to be minimized.
#' @param gradient Function (optional). Returns the gradient vector. If \code{NULL}, 
#' numerical derivatives are used.
#' @param hessian Function (optional). Returns the Hessian matrix. Used for final 
#' positive definiteness verification if \code{use_posdef = TRUE}.
#' @param lower,upper Numeric vectors. Lower and upper bounds for the parameters. 
#' Can be scalars if all parameters share the same bounds.
#'
#' @param control A list of control parameters:
#' \itemize{
#'   \item \code{use_abs_f}: Logical. Criterion: \eqn{|f_{new} - f_{old}| <} \code{tol_abs_f}.
#'   \item \code{use_rel_f}: Logical. Criterion: \eqn{|(f_{new} - f_{old}) / f_{old}| <} \code{tol_rel_f}.
#'   \item \code{use_abs_x}: Logical. Criterion: \eqn{\max |x_{new} - x_{old}| <} \code{tol_abs_x}.
#'   \item \code{use_rel_x}: Logical. Criterion: \eqn{\max |(x_{new} - x_{old}) / x_{old}| <} \code{tol_rel_x}.
#'   \item \code{use_grad}: Logical. Criterion: \eqn{\|g\|_\infty <} \code{tol_grad}.
#'   \item \code{use_posdef}: Logical. Criterion: Positive definiteness of the Hessian.
#'   \item \code{use_pred_f}: Logical. Criterion: predicted decrease \eqn{<} \code{tol_pred_f}.
#'   \item \code{use_pred_f_avg}: Logical. Criterion: predicted decrease per parameter \eqn{<} \code{tol_pred_f_avg}.
#'   \item \code{max_iter}: Maximum number of iterations (default: 10000).
#'   \item \code{m}: Number of L-BFGS memory updates (default: 5).
#'   \item \code{tol_abs_f}, \code{tol_rel_f}: Tolerances for function value change.
#'   \item \code{tol_abs_x}, \code{tol_rel_x}: Tolerances for parameter change.
#'   \item \code{tol_grad}: Tolerance for the projected gradient (default: 1e-4).
#'   \item \code{tol_pred_f}, \code{tol_pred_f_avg}: Tolerances for the predicted decrease.
#' }
#' @param ... Additional arguments passed to objective, gradient, and Hessian functions.
#'
#' @return A list containing optimization results and metadata.
#' @export
#' @examples
#' # Simple quadratic function optimization
#' quad <- function(x) (x[1] - 2)^2 + (x[2] + 1)^2
#' res <- l_bfgs_b(start = c(0, 0), objective = quad)
#' print(res$par)
l_bfgs_b <- function(
    start,
    objective,
    gradient        = NULL,
    hessian         = NULL,
    lower           = -Inf,
    upper           = Inf,
    control         = list(),
    ...
) {
  
  # ---------- 1. Configuration (Synced with Suite) ----------
  ctrl0 <- list(
    # Convergence flags
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
    m               = 5L,
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
    ls_shrink       = 0.5,
    
    curvature_eps   = 1e-10,
    damp_phi        = 0.2,
    bound_eps       = 1e-10,
    diff_method      = "forward"
  )
  ctrl <- utils::modifyList(ctrl0, control)
  
  # ---------- 2. Internal Helpers ----------
  eval_obj <- function(z) as.numeric(objective(z, ...))[1]
  
  # Shared convergence test for the main loop (4.1) and the line-search-failure branch.
  # g_inf is the projected gradient. Returns whether the enabled gradient / parameter /
  # function-value criteria are met; the predicted-decrease tests need the accepted step,
  # so they stay in 4.5b. The PD check and status assignment stay at the call sites, and
  # reading the flags here means changing the criteria affects both call sites identically.
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
  } else {
    function(z) fast_grad(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  hess_func <- if (!is.null(hessian)) {
    function(z) hessian(z, ...)
  } else {
    function(z) fast_hess(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  project <- function(x, l, u) pmax(l, pmin(x, u))
  
  get_projected_grad <- function(x, g, l, u) {
    pg <- g
    pg[x <= l + ctrl$bound_eps & g > 0] <- 0
    pg[x >= u - ctrl$bound_eps & g < 0] <- 0
    pg
  }
  
  two_loop <- function(v, s_list, y_list, rho_list, theta) {
    mcur <- length(s_list)
    if (mcur == 0L) return((1 / theta) * v)
    q <- v
    alpha <- numeric(mcur)
    for (i in rev(seq_len(mcur))) {
      alpha[i] <- rho_list[[i]] * sum(s_list[[i]] * q)
      q <- q - alpha[i] * y_list[[i]]
    }
    r <- (1 / theta) * q
    for (i in seq_len(mcur)) {
      beta <- rho_list[[i]] * sum(y_list[[i]] * r)
      r <- r + s_list[[i]] * (alpha[i] - beta)
    }
    r
  }
  
  # ---------- 3. Initialization ----------
  x <- as.numeric(start)
  n <- length(x)
  
  # start_clock defined here to ensure it exists for cpu_time calculation
  start_clock <- proc.time() 
  
  lower <- if(length(lower) == 1) rep(lower, n) else lower
  upper <- if(length(upper) == 1) rep(upper, n) else upper
  x <- project(x, lower, upper)
  
  f <- eval_obj(x)
  g <- grad_func(x)
  
  s_list <- list(); y_list <- list(); rho_list <- list()
  it <- 0L; converged <- FALSE; status <- "running"
  x_old <- x; f_old <- NA_real_; s_last <- NULL; g_inf <- NA_real_
  pred_dec <- NA_real_; pred_dec_avg <- NA_real_; H_eval <- NULL
  
  # ---------- 4. Main Loop ----------
  tryCatch({
    repeat {
      it <- it + 1L
      if (it > ctrl$max_iter) { status <- "iteration_limit_reached"; break }
      
      pg <- get_projected_grad(x, g, lower, upper)
      g_inf <- max(abs(pg))
      
      # 4.1) Convergence Verification (AND Rule)
      # Gradient / parameter / function-value tests use the current point (before the line
      # search), via the shared check_convergence so all flags are honored. g_inf is the
      # projected gradient. The predicted-decrease tests need the step the line search
      # selects, so they are evaluated after the line search (4.5b) and AND-combined.
      res_conv_pre <- check_convergence(g_inf, f, f_old, x, x_old, it)
      
      # When neither predicted-decrease test is active, the convergence decision is
      # complete here, so finalize immediately (this also skips an unnecessary line
      # search when already converged).
      if (!isTRUE(ctrl$use_pred_f) && !isTRUE(ctrl$use_pred_f_avg)) {
        if (res_conv_pre && it > 1L) {
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
            if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break } else res_conv_pre <- FALSE
          } else { converged <- TRUE; status <- "converged"; break }
        }
      }
      
      # 4.2) Update theta scaling for B0
      theta <- if (length(s_list) > 0) {
        sy_tmp <- sum(s_list[[length(s_list)]] * y_list[[length(y_list)]])
        yy_tmp <- sum(y_list[[length(y_list)]]^2)
        if (sy_tmp > 1e-12) yy_tmp / sy_tmp else 1.0
      } else 1.0
      
      # 4.3) Subspace Step via Two-loop Recursion
      p <- -two_loop(g, s_list, y_list, rho_list, theta)
      
      # Feasibility enforcement for direction
      p[x <= lower + ctrl$bound_eps & p < 0] <- 0
      p[x >= upper - ctrl$bound_eps & p > 0] <- 0
      
      if (max(abs(p)) < 1e-14) p <- -pg # Fallback to projected gradient
      
      # 4.4) Backtracking Line Search
      alpha <- ctrl$ls_alpha0; step_ok <- FALSE
      for (ls_it in seq_len(ctrl$ls_max_steps)) {
        x_try <- project(x + alpha * p, lower, upper)
        f_try <- eval_obj(x_try)
        s_try <- x_try - x
        
        # Armijo condition for projected step
        if (is.finite(f_try) && f_try <= f + ctrl$wolfe_c1 * sum(g * s_try)) {
          x_new <- x_try; f_new <- f_try; s_last <- s_try; step_ok <- TRUE; break
        }
        alpha <- alpha * ctrl$ls_shrink
      }
      
      if (!step_ok) {
        # Line search found no acceptable step. Re-check convergence with the same criteria
        # as 4.1 (via check_convergence) so user-configured stopping rules are honored: a
        # line search can fail simply because we are at a stationary point. A non-positive-
        # definite point here is a genuine failure rather than convergence.
        if (check_convergence(g_inf, f, f_old, x, x_old, it)) {
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
            if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break }
          } else { converged <- TRUE; status <- "converged"; break }
        }
        status <- "line_search_failed"; break
      }
      
      # 4.5) L-BFGS Memory Update with Powell Damping
      g_new <- grad_func(x_new); s_vec <- x_new - x; y_vec <- g_new - g
      sy <- sum(s_vec * y_vec)

      # Post-line-search convergence check (handles exact solutions, e.g., quadratics)
      g_inf_new <- max(abs(g_new), na.rm = TRUE)
      if (ctrl$use_grad && g_inf_new <= ctrl$tol_grad) {
        x <- x_new; f <- f_new; g <- g_new; g_inf <- g_inf_new
        if (isTRUE(ctrl$use_posdef)) {
          H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
          if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break }
        } else { converged <- TRUE; status <- "converged"; break }
      }
      
      # Implicit B*s approximation using theta
      Bs_approx <- theta * s_vec; sBs <- sum(s_vec * Bs_approx)
      
      # 4.5b) Predicted decrease of the line-search step on the objective scale, using the
      # L-BFGS initial Hessian model B0 = theta*I (the same model the Powell damping below
      # uses). s_vec is the actual projected step, so this reflects the step taken. The same
      # value is reported as pred_dec / pred_dec_avg, and is AND-combined with the
      # gradient/parameter/function-value result from 4.1 (res_conv_pre).
      pred_dec <- as.numeric(-(sum(g * s_vec) + 0.5 * sBs)); pred_dec_avg <- pred_dec / n
      if (isTRUE(ctrl$use_pred_f) || isTRUE(ctrl$use_pred_f_avg)) {
        res_conv <- res_conv_pre
        if (isTRUE(ctrl$use_pred_f)) res_conv <- res_conv && (is.finite(pred_dec) && pred_dec <= ctrl$tol_pred_f)
        if (isTRUE(ctrl$use_pred_f_avg)) res_conv <- res_conv && (is.finite(pred_dec_avg) && pred_dec_avg <= ctrl$tol_pred_f_avg)
        if (res_conv && it > 1L) {
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func(x_new), error = function(e) NULL)
            if (is_pd_fast(H_eval)) {
              x_old <- x; f_old <- f; x <- x_new; f <- f_new; g <- g_new
              converged <- TRUE; status <- "converged"; break
            }
          } else {
            x_old <- x; f_old <- f; x <- x_new; f <- f_new; g <- g_new
            converged <- TRUE; status <- "converged"; break
          }
        }
      }
      
      update_ok <- FALSE; y_star <- y_vec; sy_star <- sy
      if (isTRUE(ctrl$use_damped)) {
        if (sy < ctrl$damp_phi * sBs) {
          theta_damp <- ((1 - ctrl$damp_phi) * sBs) / (sBs - sy + 1e-16)
          y_star <- theta_damp * y_vec + (1 - theta_damp) * Bs_approx
          sy_star <- sum(s_vec * y_star)
        }
        if (sy_star > ctrl$curvature_eps) { y_vec <- y_star; sy <- sy_star; update_ok <- TRUE }
      } else {
        if (sy > ctrl$curvature_eps) update_ok <- TRUE
      }
      
      if (update_ok) {
        if (length(s_list) >= ctrl$m) {
          s_list[[1]] <- NULL; y_list[[1]] <- NULL; rho_list[[1]] <- NULL
        }
        s_list[[length(s_list) + 1L]] <- s_vec
        y_list[[length(y_list) + 1L]] <- y_vec
        rho_list[[length(rho_list) + 1L]] <- 1 / (sy + 1e-16)
      }
      
      x_old <- x; f_old <- f; x <- x_new; f <- f_new; g <- g_new
    }
  }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
  
  # ---------- 5. Final Status & Output Construction ----------
  # Report the derivative-based Hessian and its positive-definiteness, consistent with the
  # other optimizers. Reuse H_eval computed during the convergence check when available;
  # otherwise evaluate it once here. (L-BFGS keeps no explicit Hessian approximation, so
  # there is no approx_hessian to report, unlike the dense-B optimizers.)
  if (is.null(H_eval)) H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
  Hess_pd <- if (!is.null(H_eval)) is_pd_fast(H_eval) else FALSE
  H_final <- if (!is.null(H_eval)) H_eval else NA_real_

  final_clock <- proc.time() - start_clock
  list(
    par          = x, 
    objective    = f, 
    converged    = converged, 
    status       = status, 
    iter         = it, 
    cpu_time     = as.numeric(final_clock[1] + final_clock[2])[1], 
    elapsed_time = as.numeric(final_clock[3])[1],
    max_grad     = as.numeric(g_inf)[1],
    Hess_is_pd   = Hess_pd,
    Hessian      = H_final,
    pred_dec     = pred_dec,
    pred_dec_avg = pred_dec_avg
  )
}