#' Levenberg-Marquardt Optimization
#'
#' @description
#' Local optimizer based on the Levenberg-Marquardt algorithm. The search
#' direction blends the Newton and steepest-descent directions through an
#' adaptively damped, diagonally inflated Hessian, which keeps the local
#' quadratic model positive-definite and yields a descent direction at every
#' iteration. The function follows the common optimflex interface, so it shares
#' the same arguments, control parameters, convergence criteria, box-constraint
#' handling, and return structure as the other optimizers in the package.
#'
#' @details
#' \bold{Search direction.}
#' At each iteration the step solves \eqn{(\tilde{H}) p = -g}, where \eqn{g} is
#' the gradient and \eqn{\tilde{H}} is the Hessian (or its approximation) with an
#' inflated diagonal in the Fletcher (1971) form:
#' \deqn{\tilde{H}_{ii} = H_{ii} + da\,((1-ga)\,|H_{ii}| + ga\,tr),}
#' with \eqn{tr} the mean of \eqn{|H_{ii}|}. The damping \eqn{da} (Marquardt's
#' \eqn{\lambda}) and the curvature-mixing factor \eqn{ga} are increased until the
#' inflated matrix admits a Cholesky factorization, guaranteeing a positive-definite
#' model and a descent direction. A small \eqn{da} makes the step close to a Newton
#' step; a large \eqn{da} drives it toward steepest descent. When a step improves
#' the objective, \eqn{da} is decreased so the next step is more Newton-like; when
#' it does not, the damping is increased and a backtracking line search shrinks the
#' step length until the objective decreases.
#'
#' \bold{Hessian source (two modes).}
#' \code{hessian_update = "exact"} recomputes the Hessian at every accepted step
#' (only when an analytic Hessian is supplied). \code{hessian_update = "bfgs"}
#' (the default) refines a Hessian approximation B with damped BFGS rank-two
#' updates, using a supplied analytic Hessian (if any) only as the starting B.
#'
#' \bold{Convergence.}
#' Convergence is assessed with the package's standard criteria, combined with an
#' AND rule: gradient norm, parameter and objective-function stability, and an
#' optional positive-definiteness check at the candidate solution.
#'
#' @references
#' \itemize{
#'    \item Marquardt, D. W. (1963). An Algorithm for Least-Squares Estimation of
#'          Nonlinear Parameters. SIAM Journal.
#'    \item Fletcher, R. (1971). A modified Marquardt subroutine for nonlinear
#'          least squares.
#'    \item Philipps V., Hejblum B.P., Prague M., Commenges D., Proust-Lima C.
#'          (2021). Robust and Efficient Optimization Using a Marquardt-Levenberg
#'          Algorithm with R Package marqLevAlg. The R Journal.
#'    \item Nocedal, J., & Wright, S. J. (2006). Numerical Optimization. Springer.
#' }
#'
#' @param start Numeric vector. Starting values for the optimization parameters.
#' @param objective Function. The objective function to minimize.
#' @param gradient Function (optional). Gradient of the objective function.
#' @param hessian Function (optional). Hessian matrix of the objective function.
#' @param gn_hessian Function (optional). Gauss-Newton curvature (e.g. 2 * t(J) %*% J)
#'   used as the iteration curvature B in place of \code{hessian}, recomputed at every
#'   accepted step (damped Gauss-Newton). When supplied, the positive-definiteness check
#'   at convergence still uses \code{hessian} (the observed Hessian), so a Gauss-Newton
#'   curvature - which is positive definite by construction - does not make that check
#'   vacuous. If \code{hessian} is omitted, the check falls back to finite differences
#'   of the objective.
#' @param residual Function (optional). Residual vector r(x); together with \code{jac}, activates a Gauss-Newton least-squares mode.
#' @param jac Function (optional). Jacobian of the residual, J(x) = d r / d x.
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
#'      \item \code{hessian_update}: "bfgs" (default) or "exact".
#'    }
#' @param ... Additional arguments passed to objective, gradient, and Hessian functions.
#'
#' @return A list containing optimization results and iteration metadata.
#' @export
#' @examples
#' # Simple quadratic function optimization
#' quad <- function(x) (x[1] - 2)^2 + (x[2] + 1)^2
#' res <- levenberg_marquardt(start = c(0, 0), objective = quad)
#' print(res$par)
levenberg_marquardt <- function(
    start,
    objective,
    gradient = NULL,
    hessian  = NULL,
    gn_hessian = NULL,
    residual = NULL,
    jac      = NULL,
    lower    = -Inf,
    upper    = Inf,
    control  = list(),
    ...
) {
  # ---------- 1. Configuration ----------
  ctrl0 <- list(
    use_abs_f       = FALSE,
    use_rel_f       = FALSE,
    use_abs_x       = FALSE,
    use_rel_x       = TRUE,
    use_grad        = TRUE,
    use_posdef      = TRUE,
    use_pred_f      = FALSE,
    use_pred_f_avg  = FALSE,
    
    max_iter        = 10000L,
    tol_abs_f       = 1e-6,
    tol_rel_f       = 1e-6,
    tol_abs_x       = 1e-6,
    tol_rel_x       = 1e-6,
    tol_grad        = 1e-4,
    tol_pred_f      = 1e-4,
    tol_pred_f_avg  = 1e-4,
    
    H_init_diag     = 1.0,
    diff_method     = "forward",
    hessian_update  = "bfgs",
    use_damped      = TRUE,
    damp_phi        = 0.2,
    
    # Fletcher diagonal-inflation damping controls
    da_init         = 1e-2,   # initial Marquardt damping (lambda)
    ga_init         = 0.01,   # initial curvature-mixing factor (eta)
    da_factor       = 5.0,    # multiplicative inflation factor
    da_min          = 1e-7,   # floor for da before a successful step
    gonfle_max      = 10L,    # max diagonal inflations per iteration
    
    # backtracking line search controls
    ls_max          = 20L,    # max backtracking steps
    ls_shrink       = 0.5,    # step-length shrink factor
    ls_min_step     = 1e-14   # smallest step length before giving up
  )
  ctrl <- utils::modifyList(ctrl0, control)
  ctrl$hessian_update <- match.arg(ctrl$hessian_update, c("bfgs", "exact"))

  # Least-squares mode is used only when both residual and jac are supplied.
  ls_mode <- !is.null(residual) && !is.null(jac)
  
  # ---------- 2. Internal Helpers ----------
  eval_obj <- function(z) as.numeric(objective(z, ...))[1]
  
  grad_func <- if (!is.null(gradient)) {
    function(z) as.numeric(gradient(z, ...))
  } else {
    function(z) fast_grad(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  # Iteration curvature B: prefer an explicit Gauss-Newton curvature when supplied,
  # otherwise the supplied Hessian, otherwise finite differences of the objective.
  hess_func <- if (!is.null(gn_hessian)) {
    function(z) gn_hessian(z, ...)
  } else if (!is.null(hessian)) {
    function(z) hessian(z, ...)
  } else {
    function(z) fast_hess(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  # Hessian used for the positive-definiteness check. Always the observed Hessian
  # (a supplied analytic Hessian, else finite differences) - never gn_hessian, which is
  # positive definite by construction and would make the check vacuous.
  hess_func_pd <- if (!is.null(hessian)) {
    function(z) hessian(z, ...)
  } else {
    function(z) fast_hess(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  project <- function(z, l, u) pmax(l, pmin(z, u))
  
  # Fletcher diagonal inflation of a symmetric matrix M.
  # Increases da (and ga, after a few tries) until chol() succeeds, then returns
  # the Cholesky factor R (upper triangular) of the inflated matrix together with
  # the da/ga that worked. tr is the mean of the absolute diagonal of M.
  inflate_and_chol <- function(M, da, ga) {
    n_ <- nrow(M)
    diagM <- diag(M)
    tr <- mean(abs(diagM))
    if (!is.finite(tr) || tr == 0) tr <- 1
    
    build <- function(da_, ga_) {
      Mi <- M
      add <- ifelse(diagM != 0,
                    da_ * ((1 - ga_) * abs(diagM) + ga_ * tr),
                    da_ * ga_ * tr)
      diag(Mi) <- diagM + add
      Mi
    }
    
    ncount <- 0L
    repeat {
      Mi <- build(da, ga)
      R <- tryCatch(chol(0.5 * (Mi + t(Mi))), error = function(e) NULL)
      if (!is.null(R)) {
        return(list(R = R, da = da, ga = ga, ok = TRUE))
      }
      ncount <- ncount + 1L
      if (ncount <= 3L || ga >= 1) {
        da <- da * ctrl$da_factor
      } else {
        ga <- ga * ctrl$da_factor
        if (ga > 1) ga <- 1
      }
      if (ncount >= ctrl$gonfle_max) {
        # Last resort: a strongly diagonally dominant fallback that is PD.
        Mi <- build(da, ga)
        diag(Mi) <- diag(Mi) + (abs(min(diag(Mi))) + 1) 
        R <- tryCatch(chol(0.5 * (Mi + t(Mi))), error = function(e) NULL)
        return(list(R = R, da = da, ga = ga, ok = !is.null(R)))
      }
    }
  }
  
  # ---------- 3. Initialization ----------
  x <- project(as.numeric(start), lower, upper)
  n <- length(x)
  
  # start_clock defined here for CPU time calculation
  start_clock <- proc.time()
  
  f <- tryCatch(eval_obj(x), error = function(e) NA_real_)
  it <- 0L; converged <- FALSE; status <- "running"
  x_old <- x; f_old <- NA_real_
  
  # Hessian-source mode flags.
  # use_exact_hess: an iteration curvature was supplied - either a Gauss-Newton curvature
  #   (gn_hessian) or a Hessian - and is used to initialize B. The positive-definiteness
  #   check at convergence always uses the observed Hessian (hessian / finite differences),
  #   never gn_hessian.
  # update_exact: recompute the iteration curvature exactly at every accepted step. True
  #   when a Gauss-Newton curvature is supplied, or when a Hessian is supplied AND
  #   hessian_update == "exact". Otherwise (the default "bfgs") B is refined with damped
  #   BFGS rank-two updates from the supplied starting Hessian.
  use_exact_hess <- !is.null(gn_hessian) || !is.null(hessian)
  update_exact <- use_exact_hess && (!is.null(gn_hessian) || identical(ctrl$hessian_update, "exact"))
  
  # Initialize Hessian (approximation) B.
  B <- if (use_exact_hess) {
    B_init <- tryCatch(hess_func(x), error = function(e) diag(ctrl$H_init_diag, n))
    B_init <- 0.5 * (B_init + t(B_init))
    if (!is_pd_fast(B_init)) {
      ev <- eigen(B_init, symmetric = TRUE, only.values = TRUE)$values
      shift <- max(abs(min(ev)) + 1e-6, max(abs(ev)) * 1e-7)
      B_init + diag(shift, n)
    } else B_init
  } else {
    diag(ctrl$H_init_diag, n)
  }
  
  H_eval <- NULL; g_inf <- NA_real_
  pred_dec <- NA_real_; pred_dec_avg <- NA_real_
  scaled_B <- FALSE   # whether the one-time self-scaling of B has been applied (BFGS mode)
  da <- ctrl$da_init   # current Marquardt damping, carried across iterations
  
  # ---------- 4. Main Loop ----------
  if (!is.finite(f)) {
    status <- "objective_error_at_start"
  } else if (ls_mode) {
    # ---------- Least-squares mode (residual + jac): damped Gauss-Newton via QR ----------
    # Activated only when both residual and jac are supplied. g = 2 J'r, and the damped
    # Gauss-Newton step solves (2 J'J + da I) p = -g as a linear least-squares problem in J
    # augmented with sqrt(da) I, so 2 J'J is not assembled for the solve. The Marquardt
    # damping da is decreased after a successful step and increased after a failed one,
    # with the same convergence tests, line search, and positive-definiteness check.
    r <- as.numeric(residual(x, ...))
    J <- matrix(as.numeric(jac(x, ...)), nrow = length(r))
    g <- 2 * as.numeric(crossprod(J, r))

    tryCatch({
      repeat {
        if (it >= ctrl$max_iter) { status <- "iteration_limit_reached"; break }
        it <- it + 1L
        g_inf <- max(abs(g), na.rm = TRUE)

        # Damped Gauss-Newton step via augmented QR of [sqrt(2) J; sqrt(da) I].
        p_step <- tryCatch({
          A_aug <- rbind(sqrt(2) * J, sqrt(da) * diag(n))
          qr.solve(A_aug, c(-sqrt(2) * r, rep(0, n)))
        }, error = function(e) NULL)
        if (is.null(p_step) || any(!is.finite(p_step))) {
          current_pred_dec <- 0
        } else {
          Jp <- as.numeric(J %*% p_step)
          current_pred_dec <- as.numeric(-(sum(g * p_step) + sum(Jp^2)))
        }

        pred_dec <- current_pred_dec
        pred_dec_avg <- current_pred_dec / n
        res_conv <- TRUE
        if (ctrl$use_grad) res_conv <- res_conv && (g_inf <= ctrl$tol_grad)
        if (ctrl$use_abs_f && !is.na(f_old)) res_conv <- res_conv && (abs(f - f_old) <= ctrl$tol_abs_f)
        if (ctrl$use_rel_f && !is.na(f_old)) res_conv <- res_conv && (abs((f - f_old) / max(1, abs(f_old))) <= ctrl$tol_rel_f)
        if (ctrl$use_abs_x && it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) <= ctrl$tol_abs_x)
        if (ctrl$use_rel_x && it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) / max(1, max(abs(x_old)))) <= ctrl$tol_rel_x
        if (isTRUE(ctrl$use_pred_f)) res_conv <- res_conv && (is.finite(pred_dec) && pred_dec <= ctrl$tol_pred_f)
        if (isTRUE(ctrl$use_pred_f_avg)) res_conv <- res_conv && (is.finite(pred_dec_avg) && pred_dec_avg <= ctrl$tol_pred_f_avg)

        if (res_conv && it > 1L) {
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
            if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break } else res_conv <- FALSE
          } else { converged <- TRUE; status <- "converged"; break }
        }

        # Step acceptance with backtracking line search.
        accepted <- FALSE; x_try <- x; f_try <- f
        if (!is.null(p_step) && any(p_step != 0)) {
          step_len <- 1.0; ls <- 0L
          repeat {
            x_cand <- project(x + step_len * p_step, lower, upper)
            f_cand <- tryCatch(eval_obj(x_cand), error = function(e) NA_real_)
            if (is.finite(f_cand) && f_cand < f) { x_try <- x_cand; f_try <- f_cand; accepted <- TRUE; break }
            ls <- ls + 1L; step_len <- step_len * ctrl$ls_shrink
            if (ls >= ctrl$ls_max || step_len < ctrl$ls_min_step) break
          }
        }
        actual_red <- f - f_try

        if (accepted && actual_red > 0) {
          # Successful step: decrease damping (toward a more Newton-like step), with a floor.
          da <- if (da < ctrl$da_min) ctrl$da_min else da / (ctrl$da_factor + 2)
          x_old <- x; f_old <- f; x <- x_try; f <- f_try
          r <- as.numeric(residual(x, ...)); J <- matrix(as.numeric(jac(x, ...)), nrow = length(r))
          g <- 2 * as.numeric(crossprod(J, r))

          # Post-step convergence check (handles exact solutions).
          g_inf_new <- max(abs(g), na.rm = TRUE)
          if (ctrl$use_grad && g_inf_new <= ctrl$tol_grad) {
            g_inf <- g_inf_new
            if (isTRUE(ctrl$use_posdef)) {
              H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
              if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break }
            } else { converged <- TRUE; status <- "converged"; break }
          }
        } else {
          # Failed step: increase damping so the next model is more conservative.
          da <- da * ctrl$da_factor
          if (da > 1e12) { status <- "damping_too_large"; break }
        }
      }
    }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
    B <- 2 * crossprod(J)   # Gauss-Newton Hessian for the approx_hessian return field
  } else {
    g <- grad_func(x)
    
    tryCatch({
      repeat {
        if (it >= ctrl$max_iter) { status <- "iteration_limit_reached"; break }
        it <- it + 1L
        
        # 4.1) Free Variables Identification (for Box Constraints)
        is_free <- !((x <= lower + 1e-10 & g > 0) | (x >= upper - 1e-10 & g < 0))
        free_idx <- which(is_free); nfree <- length(free_idx)
        
        if (nfree > 0L) {
          g_f <- g[free_idx]; g_inf <- max(abs(g_f), na.rm = TRUE)
          B_f <- B[free_idx, free_idx, drop = FALSE]
          
          # 4.2) Damped step: inflate the diagonal (Fletcher) until Cholesky
          #      succeeds, then solve (inflated B_f) p = -g_f via that factor.
          infl <- inflate_and_chol(B_f, da, ctrl$ga_init)
          R_f <- infl$R
          ga_used <- infl$ga
          da_used <- infl$da
          
          if (is.null(R_f)) {
            # Could not build a PD model; treat as a failed step (shrink via da).
            p_f <- rep(0, nfree)
            current_pred_dec <- 0
          } else {
            # solve (R'R) p = -g_f
            p_f <- backsolve(R_f, forwardsolve(t(R_f), -g_f))
            current_pred_dec <- as.numeric(-(sum(g_f * p_f) + 0.5 * sum(p_f * (B_f %*% p_f))))
          }
        } else {
          g_inf <- 0; current_pred_dec <- 0
          p_f <- numeric(0); ga_used <- ctrl$ga_init; da_used <- da
        }
        
        # 4.3) Convergence Check
        # The full set of convergence criteria is combined with an AND rule. All
        # tests use the current point (before the step is taken).
        pred_dec <- current_pred_dec
        pred_dec_avg <- current_pred_dec / n
        res_conv <- TRUE
        if (ctrl$use_grad) res_conv <- res_conv && (g_inf <= ctrl$tol_grad)
        if (ctrl$use_abs_f && !is.na(f_old)) res_conv <- res_conv && (abs(f - f_old) <= ctrl$tol_abs_f)
        if (ctrl$use_rel_f && !is.na(f_old)) res_conv <- res_conv && (abs((f - f_old) / max(1, abs(f_old))) <= ctrl$tol_rel_f)
        if (ctrl$use_abs_x && it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) <= ctrl$tol_abs_x)
        if (ctrl$use_rel_x && it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) / max(1, max(abs(x_old)))) <= ctrl$tol_rel_x
        if (isTRUE(ctrl$use_pred_f)) res_conv <- res_conv && (is.finite(pred_dec) && pred_dec <= ctrl$tol_pred_f)
        if (isTRUE(ctrl$use_pred_f_avg)) res_conv <- res_conv && (is.finite(pred_dec_avg) && pred_dec_avg <= ctrl$tol_pred_f_avg)
        
        if (res_conv && it > 1L) {
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
            if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break } else res_conv <- FALSE
          } else { converged <- TRUE; status <- "converged"; break }
        }
        
        # 4.4) Step Acceptance with backtracking line search
        # Try the full step; if it does not reduce f, backtrack by halving until it
        # does or the step length becomes too small.
        p_full <- rep(0, n); if (nfree > 0L) p_full[free_idx] <- p_f
        
        step_len <- 1.0
        accepted <- FALSE
        x_try <- x; f_try <- f
        if (nfree > 0L && any(p_full != 0)) {
          ls <- 0L
          repeat {
            x_cand <- project(x + step_len * p_full, lower, upper)
            f_cand <- tryCatch(eval_obj(x_cand), error = function(e) NA_real_)
            if (is.finite(f_cand) && f_cand < f) {
              x_try <- x_cand; f_try <- f_cand; accepted <- TRUE
              break
            }
            ls <- ls + 1L
            step_len <- step_len * ctrl$ls_shrink
            if (ls >= ctrl$ls_max || step_len < ctrl$ls_min_step) break
          }
        }
        
        actual_red <- f - f_try
        
        if (accepted && actual_red > 0) {
          # Successful step: decrease damping (toward a more Newton-like step), with a floor.
          da <- if (da_used < ctrl$da_min) ctrl$da_min else da_used / (ctrl$da_factor + 2)
          
          g_new <- grad_func(x_try); s <- x_try - x; y <- g_new - g
          
          # Self-scale the initial B on the first BFGS update
          # (Nocedal & Wright 2006, eq. 6.20; direct-Hessian form B0 <- (y'y / s'y) I).
          if (!use_exact_hess && !scaled_B) {
            yy0 <- sum(y * y); sy0 <- sum(s * y)
            if (is.finite(yy0) && yy0 > 1e-12 && is.finite(sy0) && sy0 > 1e-12) {
              B <- diag(yy0 / sy0, n)
              scaled_B <- TRUE
            }
          }
          
          # Damped BFGS Update
          Bs <- as.numeric(B %*% s); sBs <- sum(s * Bs); sy <- sum(s * y)
          update_ok <- FALSE; y_star <- y; sy_star <- sy
          if (isTRUE(ctrl$use_damped)) {
            if (sy < ctrl$damp_phi * sBs) {
              theta_damp <- ((1 - ctrl$damp_phi) * sBs) / (sBs - sy)
              y_star <- theta_damp * y + (1 - theta_damp) * Bs; sy_star <- sum(s * y_star)
            }
            if (sy_star > 1e-12) { y <- y_star; sy <- sy_star; update_ok <- TRUE }
          } else { if (sy > 1e-12) update_ok <- TRUE }
          
          if (update_exact) {
            B_new <- tryCatch(hess_func(x_try), error = function(e) B)
            B <- 0.5 * (B_new + t(B_new))
          } else {
            if (update_ok) {
              B <- B - (Bs %*% t(Bs)) / (sBs + 1e-16) + (y %*% t(y)) / (sy + 1e-16)
              B <- 0.5 * (B + t(B))
            }
          }
          
          x_old <- x; f_old <- f; x <- x_try; f <- f_try; g <- g_new
          
          # Post-step convergence check (handles exact solutions, e.g., quadratics)
          g_inf_new <- max(abs(g_new), na.rm = TRUE)
          if (ctrl$use_grad && g_inf_new <= ctrl$tol_grad) {
            g_inf <- g_inf_new
            if (isTRUE(ctrl$use_posdef)) {
              H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
              if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break }
            } else { converged <- TRUE; status <- "converged"; break }
          }
        } else {
          # Failed step: increase damping so the next model is more conservative.
          da <- da_used * ctrl$da_factor
          if (da > 1e12) { status <- "damping_too_large"; break }
        }
      }
    }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
  }
  
  # ---------- 5. Final Reporting ----------
  if (is.null(H_eval)) H_eval <- tryCatch(hess_func(x), error = function(e) NULL)
  final_clock <- proc.time() - start_clock
  
  list(
    par              = x,
    objective        = f,
    converged        = converged,
    status           = status,
    iter             = it,
    cpu_time         = as.numeric(final_clock[1] + final_clock[2]),
    elapsed_time     = as.numeric(final_clock[3]),
    max_grad         = as.numeric(g_inf),
    Hessian          = H_eval,
    approx_hessian   = B,
    pred_dec         = pred_dec,
    pred_dec_avg     = pred_dec_avg
  )
}