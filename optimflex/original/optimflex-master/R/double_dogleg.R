#' Double Dogleg Trust-Region Optimization
#'
#' @description
#' Implements the Double Dogleg Trust-Region algorithm for non-linear optimization.
#'
#' @details
#' This function implements the Double Dogleg method within a Trust-Region framework, 
#' primarily based on the work of Dennis and Mei (1979). 
#' 
#' \bold{Trust-Region vs. Line Search:}
#' While Line Search methods (like BFGS) first determine a search direction and then 
#' find an appropriate step length, Trust-Region methods define a neighborhood 
#' around the current point (the trust region with radius \eqn{\Delta}) where a local 
#' quadratic model is assumed to be reliable. The algorithm then finds a step that 
#' minimizes this model within the radius. This approach is generally more robust, 
#' especially when the Hessian is not positive definite.
#'
#' \bold{Powell's Dogleg vs. Double Dogleg:}
#' Powell's original Dogleg method (1970) constructs a trajectory consisting of 
#' two line segments: one from the current point to the Cauchy point, and another 
#' from the Cauchy point to the Newton point. The "Double Dogleg" modification 
#' by Dennis and Mei (1979) introduces an intermediate "bias" point (\eqn{p_W}) 
#' along the Newton direction.
#' \itemize{
#'   \item \bold{Cauchy Point (\eqn{p_C}):} The minimizer of the quadratic model along 
#'         the steepest descent direction.
#'   \item \bold{Newton Point (\eqn{p_N}):} The minimizer of the quadratic model \eqn{(B^{-1}g)}.
#'   \item \bold{Double Dogleg Point (\eqn{p_W}):} A point defined as \eqn{\gamma \cdot p_N}, 
#'         where \eqn{\gamma} is a scaling factor (bias) that ensures the path stays 
#'         closer to the Newton direction while maintaining monotonic descent in 
#'         the model.
#' }
#' This modification allows the algorithm to perform more like a Newton method 
#' earlier in the optimization process compared to the standard Dogleg.
#'
#' \bold{BFGS curvature as a Cholesky factor:}
#' When \code{hessian_update == "bfgs"} the iteration curvature is maintained as the
#' upper-triangular Cholesky factor \eqn{R} of \eqn{B} (so \eqn{B = R^T R}) and refined by
#' a factored (dual) damped BFGS update: a rank-one Cholesky update by
#' \eqn{y / \sqrt{s^T y}} followed by a rank-one Cholesky downdate by
#' \eqn{B s / \sqrt{s^T B s}}. The Newton point is obtained from two triangular solves,
#' and positive definiteness is structural (no eigenvalue safeguard needed on the BFGS path).
#'
#' @references
#' \itemize{
#'   \item Dennis, J. E., & Mei, H. H. (1979). Two New Unconstrained Optimization 
#'         Algorithms which use Function and Gradient Values. 
#'         \emph{Journal of Optimization Theory and Applications}, 28(4), 453-482.
#'   \item Powell, M. J. D. (1970). A Hybrid Method for Nonlinear Equations. 
#'         \emph{Numerical Methods for Nonlinear Algebraic Equations}.
#'   \item Nocedal, J., & Wright, S. J. (2006). \emph{Numerical Optimization}. Springer.
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
#' @param control List. Control parameters including convergence flags starting with 'use_'.
#'   \itemize{
#'     \item \code{use_abs_f}: Logical. Use absolute change in objective for convergence.
#'     \item \code{use_rel_f}: Logical. Use relative change in objective for convergence.
#'     \item \code{use_abs_x}: Logical. Use absolute change in parameters for convergence.
#'     \item \code{use_rel_x}: Logical. Use relative change in parameters for convergence.
#'     \item \code{use_grad}: Logical. Use gradient norm for convergence.
#'     \item \code{use_posdef}: Logical. Verify positive definiteness at convergence.
#'     \item \code{use_pred_f}: Logical. Record predicted objective decrease.
#'     \item \code{use_pred_f_avg}: Logical. Record average predicted decrease.
#'   }
#' @param ... Additional arguments passed to objective, gradient, and Hessian functions.
#'
#' @return A list containing optimization results and iteration metadata.
#' @export
#' @examples
#' # Simple quadratic function optimization
#' quad <- function(x) (x[1] - 2)^2 + (x[2] + 1)^2
#' res <- double_dogleg(start = c(0, 0), objective = quad)
#' print(res$par)
double_dogleg <- function(
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
    initial_delta   = 1.0, 
    delta_max       = 100.0, 
    dd_bias         = 0.8, 
    rho_accept      = 0.1, 
    rho_expand      = 0.75,
    delta_shrink    = 0.25, 
    delta_expand    = 2.0, 
    H_init_diag     = 1.0, 
    ridge_eps       = 1e-9, 
    diff_method      = "forward",
    hessian_update   = "bfgs",
    use_damped      = TRUE, 
    damp_phi        = 0.2
  )
  ctrl <- utils::modifyList(ctrl0, control)
  ctrl$diff_method <- match.arg(ctrl$diff_method, c("forward", "central", "richardson"))
  ctrl$hessian_update <- match.arg(ctrl$hessian_update, c("bfgs", "exact"))
  
  # Least-squares mode is used only when both residual and jac are supplied.
  ls_mode <- !is.null(residual) && !is.null(jac)
  
  if (ctrl$diff_method == "richardson") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' required.")
  }
  
  # ---------- 2. Internal Helpers ----------
  eval_obj <- function(z) as.numeric(objective(z, ...))[1]
  
  grad_func <- if (!is.null(gradient)) {
    function(z) as.numeric(gradient(z, ...))
  } else if (ctrl$diff_method == "richardson") {
    function(z) as.numeric(numDeriv::grad(objective, z, method = "Richardson", ...))
  } else {
    function(z) fast_grad(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  # Iteration curvature B: prefer an explicit Gauss-Newton curvature when supplied,
  # otherwise the supplied Hessian, otherwise finite differences of the objective.
  hess_func <- if (!is.null(gn_hessian)) {
    function(z) gn_hessian(z, ...)
  } else if (!is.null(hessian)) {
    function(z) hessian(z, ...)
  } else if (ctrl$diff_method == "richardson") {
    function(z) numDeriv::hessian(objective, z, method = "Richardson", ...)
  } else {
    function(z) fast_hess(objective, z, diff_method = ctrl$diff_method, ...)
  }
  # Hessian used for the positive-definiteness check and the reported H_eval. Always the
  # observed Hessian (a supplied analytic Hessian, else finite differences) - never
  # gn_hessian, which is positive definite by construction and would make the check vacuous.
  hess_func_pd <- if (!is.null(hessian)) {
    function(z) hessian(z, ...)
  } else if (ctrl$diff_method == "richardson") {
    function(z) numDeriv::hessian(objective, z, method = "Richardson", ...)
  } else {
    function(z) fast_hess(objective, z, diff_method = ctrl$diff_method, ...)
  }
  
  # Cholesky helpers for the factored (dual) BFGS curvature.
  # Convention: upper-triangular R with B = R^T R, matching base R's chol().
  # Positive-definite Cholesky with a shift fallback (mirrors the eigenvalue safeguard).
  chol_psd <- function(M) {
    M <- 0.5 * (M + t(M))
    R <- tryCatch(chol(M), error = function(e) NULL)
    if (!is.null(R)) return(R)
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    shift <- max(abs(min(ev)) + 1e-6, max(abs(ev)) * 1e-7)
    chol(M + diag(shift, nrow(M)))
  }
  # Rank-one Cholesky update: returns R~ with R~^T R~ = R^T R + x x^T.
  chol_update <- function(R, x) {
    x <- as.numeric(x); m <- length(x)
    for (k in seq_len(m)) {
      Rkk <- R[k, k]
      rr  <- sqrt(Rkk * Rkk + x[k] * x[k])
      cc  <- rr / Rkk
      ss  <- x[k] / Rkk
      R[k, k] <- rr
      if (k < m) {
        idx <- (k + 1):m
        R[k, idx] <- (R[k, idx] + ss * x[idx]) / cc
        x[idx]    <- cc * x[idx] - ss * R[k, idx]
      }
    }
    R
  }
  # Rank-one Cholesky downdate: returns R~ with R~^T R~ = R^T R - x x^T, or NULL if not PD.
  chol_downdate <- function(R, x) {
    x <- as.numeric(x); m <- length(x)
    for (k in seq_len(m)) {
      Rkk <- R[k, k]
      dd  <- Rkk * Rkk - x[k] * x[k]
      if (!is.finite(dd) || dd <= 0) return(NULL)
      rr <- sqrt(dd)
      cc <- rr / Rkk
      ss <- x[k] / Rkk
      R[k, k] <- rr
      if (k < m) {
        idx <- (k + 1):m
        R[k, idx] <- (R[k, idx] - ss * x[idx]) / cc
        x[idx]    <- cc * x[idx] - ss * R[k, idx]
      }
    }
    R
  }
  
  # ---------- 3. Initialization ----------
  x <- pmax(lower, pmin(as.numeric(start), upper))
  n <- length(x)
  
  start_clock <- proc.time() 
  f <- tryCatch(eval_obj(x), error = function(e) NA_real_)
  
  it <- 0L; converged <- FALSE; status <- "running"
  x_old <- x; f_old <- NA_real_; delta <- ctrl$initial_delta
  
  # Initialize Hessian approximation
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
  scaled_B <- FALSE   # whether the one-time self-scaling of B has been applied
  R <- NULL           # Cholesky factor of the iteration curvature (B = R^T R), set in the general loop
  
  # ---------- 4. Main Loop ----------
  if (!is.finite(f)) {
    status <- "objective_error_at_start"
  } else if (ls_mode) {
    # ---------- Least-squares mode (residual + jac): Gauss-Newton double dogleg via QR ----------
    # Activated only when both residual and jac are supplied. The Gauss-Newton point is the
    # least-squares solution of min || -r - J p || obtained from J (2 J'J is not assembled
    # for the solve); g = 2 J'r, the curvature is 2 J'J, and the double-dogleg bias point
    # pW = gamma * pN uses gamma = dd_bias * ||g||^4 / ((g' B g)(g' B^{-1} g)) with
    # g' B^{-1} g = -g' pN. The path is restricted to the trust region of radius delta.
    r <- as.numeric(residual(x, ...))
    J <- matrix(as.numeric(jac(x, ...)), nrow = length(r))
    g <- 2 * as.numeric(crossprod(J, r))
    tryCatch({
      repeat {
        if (it >= ctrl$max_iter) { status <- "iteration_limit_reached"; break }
        it <- it + 1L
        g_inf <- max(abs(g), na.rm = TRUE)
        
        # Gauss-Newton point via QR; small ridge fallback if rank-deficient.
        pN <- tryCatch(qr.solve(J, -r), error = function(e) {
          A <- rbind(J, sqrt(ctrl$ridge_eps) * diag(n)); qr.solve(A, c(-r, rep(0, n)))
        })
        gnorm <- sqrt(sum(g^2)); Jg <- as.numeric(J %*% g); gBg <- 2 * sum(Jg^2)
        alpha_c <- if (gBg > 1e-15) (gnorm^2) / gBg else delta / max(gnorm, 1e-12)
        pC <- -alpha_c * g
        
        ghinvg <- sum(g * (-pN))                                   # g' B^{-1} g = -g' pN
        gamma <- if (ghinvg > 1e-15 && gBg > 1e-15) ctrl$dd_bias * (gnorm^4 / (gBg * ghinvg)) else 1.0
        gamma <- max(alpha_c, min(1.0, gamma))
        pW <- gamma * pN
        
        nPN <- sqrt(sum(pN^2)); nPC <- sqrt(sum(pC^2)); nPW <- sqrt(sum(pW^2))
        if (nPN <= delta) {
          p <- pN
        } else if (nPC >= delta) {
          p <- (delta / nPC) * pC
        } else if (nPW <= delta) {
          d <- (pN - pW); aa <- sum(d^2); bb <- 2 * sum(pW * d); cc <- sum(pW^2) - delta^2
          tau <- (-bb + sqrt(max(0, bb^2 - 4 * aa * cc))) / (2 * (aa + 1e-16))
          p <- pW + tau * d
        } else {
          d <- (pW - pC); aa <- sum(d^2); bb <- 2 * sum(pC * d); cc <- sum(pC^2) - delta^2
          tau <- (-bb + sqrt(max(0, bb^2 - 4 * aa * cc))) / (2 * (aa + 1e-16))
          p <- pC + tau * d
        }
        Jp <- as.numeric(J %*% p)
        current_pred_dec <- as.numeric(-(sum(g * p) + sum(Jp^2)))
        
        pred_dec <- current_pred_dec; pred_dec_avg <- current_pred_dec / n
        res_conv <- TRUE
        if (ctrl$use_grad) res_conv <- res_conv && (g_inf <= ctrl$tol_grad)
        if (ctrl$use_abs_f && !is.na(f_old)) res_conv <- res_conv && (abs(f - f_old) <= ctrl$tol_abs_f)
        if (ctrl$use_rel_f && !is.na(f_old)) res_conv <- res_conv && (abs((f - f_old) / max(1, abs(f_old))) <= ctrl$tol_rel_f)
        if (ctrl$use_abs_x && it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) <= ctrl$tol_abs_x)
        if (ctrl$use_rel_x && it > 1L) res_conv <- res_conv && (max(abs(x - x_old)) / max(1, max(abs(x_old)))) <= ctrl$tol_rel_x
        if (isTRUE(ctrl$use_pred_f)) res_conv <- res_conv && (is.finite(pred_dec) && pred_dec <= ctrl$tol_pred_f)
        if (isTRUE(ctrl$use_pred_f_avg)) res_conv <- res_conv && (is.finite(pred_dec_avg) && pred_dec_avg <= ctrl$tol_pred_f_avg)
        
        if (res_conv && it > 1) {
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
            if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break } else res_conv <- FALSE
          } else { converged <- TRUE; status <- "converged"; break }
        }
        
        x_try <- pmax(lower, pmin(x + p, upper)); f_try <- eval_obj(x_try)
        actual_red <- f - f_try
        rho <- if (is.finite(current_pred_dec) && current_pred_dec > 1e-15) actual_red / current_pred_dec else 0
        
        if (rho > ctrl$rho_accept && actual_red > 0) {
          x_old <- x; f_old <- f; x <- x_try; f <- f_try
          r <- as.numeric(residual(x, ...)); J <- matrix(as.numeric(jac(x, ...)), nrow = length(r))
          g <- 2 * as.numeric(crossprod(J, r))
          g_inf_new <- max(abs(g), na.rm = TRUE)
          if (ctrl$use_grad && g_inf_new <= ctrl$tol_grad) {
            g_inf <- g_inf_new
            if (isTRUE(ctrl$use_posdef)) {
              H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
              if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break }
            } else { converged <- TRUE; status <- "converged"; break }
          }
          if (rho > ctrl$rho_expand) delta <- min(ctrl$delta_max, ctrl$delta_expand * delta)
        } else {
          delta <- ctrl$delta_shrink * delta
          if (delta < 1e-14) { status <- "radius_too_small"; break }
        }
      }
    }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
    B <- 2 * crossprod(J)   # Gauss-Newton Hessian for the approx_hessian return field
  } else {
    g <- grad_func(x)
    R <- chol_psd(B)        # maintain the curvature as its Cholesky factor (B = R^T R)
    
    tryCatch({
      repeat {
        if (it >= ctrl$max_iter) { status <- "iteration_limit_reached"; break }
        it <- it + 1L
        
        # 4.1) Free Variables Identification (for Box Constraints)
        is_free <- !((x <= lower + 1e-10 & g > 0) | (x >= upper - 1e-10 & g < 0))
        free_idx <- which(is_free)
        nfree <- length(free_idx)
        
        if (nfree > 0L) {
          g_f <- g[free_idx]
          g_inf <- max(abs(g_f), na.rm = TRUE)
          # Cholesky factor of the free-subspace curvature: use R directly when every
          # variable is free, otherwise factor the free principal submatrix
          # B_f = R[,free]' R[,free] (positive definite because B is).
          R_f <- if (nfree == n) R else chol_psd(crossprod(R[, free_idx, drop = FALSE]))
          
          # 4.2) Double Dogleg Subproblem Step (solved through the Cholesky factor)
          pN_f <- as.numeric(backsolve(R_f, forwardsolve(t(R_f), -g_f)))
          
          gnorm <- sqrt(sum(g_f^2))
          Rg_f  <- as.numeric(R_f %*% g_f)
          gBg   <- sum(Rg_f^2)                                      # g' B_f g = ||R_f g||^2
          alpha_c <- if (gBg > 1e-15) (gnorm^2) / gBg else delta / max(gnorm, 1e-12)
          pC_f <- -alpha_c * g_f
          
          ghinvg <- sum(g_f * (-pN_f))
          gamma <- if (ghinvg > 1e-15 && gBg > 1e-15) ctrl$dd_bias * (gnorm^4 / (gBg * ghinvg)) else 1.0
          gamma <- max(alpha_c, min(1.0, gamma))
          pW_f <- gamma * pN_f
          
          nPN <- sqrt(sum(pN_f^2)); nPC <- sqrt(sum(pC_f^2)); nPW <- sqrt(sum(pW_f^2))
          
          # Interpolation logic for step selection
          if (nPN <= delta) {
            p_f <- pN_f
          } else if (nPC >= delta) {
            p_f <- (delta / nPC) * pC_f
          } else if (nPW <= delta) {
            d <- (pN_f - pW_f); aa <- sum(d^2); bb <- 2 * sum(pW_f * d); cc <- sum(pW_f^2) - delta^2
            tau <- (-bb + sqrt(max(0, bb^2 - 4 * aa * cc))) / (2 * (aa + 1e-16))
            p_f <- pW_f + tau * d
          } else {
            d <- (pW_f - pC_f); aa <- sum(d^2); bb <- 2 * sum(pC_f * d); cc <- sum(pC_f^2) - delta^2
            tau <- (-bb + sqrt(max(0, bb^2 - 4 * aa * cc))) / (2 * (aa + 1e-16))
            p_f <- pC_f + tau * d
          }
          Rp_f <- as.numeric(R_f %*% p_f)
          current_pred_dec <- as.numeric(-(sum(g_f * p_f) + 0.5 * sum(Rp_f^2)))
        } else {
          g_inf <- 0; current_pred_dec <- 0
        }
        
        # 4.3) Convergence Verification
        # All convergence criteria are combined with an AND rule, matching gauss_newton().
        # The predicted-decrease tests reuse current_pred_dec from 4.2 (the quadratic-model
        # predicted decrease of the double-dogleg step on the free subspace); the same value
        # is reported as pred_dec / pred_dec_avg.
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
        
        if (res_conv && it > 1) {
          if (isTRUE(ctrl$use_posdef)) {
            H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
            if (is_pd_fast(H_eval)) { converged <- TRUE; status <- "converged"; break } else res_conv <- FALSE
          } else { converged <- TRUE; status <- "converged"; break }
        }
        
        # 4.4) Step Acceptance & Trust Region Management
        p_full <- rep(0, n)
        if (nfree > 0L) p_full[free_idx] <- p_f
        x_try <- pmax(lower, pmin(x + p_full, upper))
        f_try <- eval_obj(x_try)
        actual_red <- f - f_try
        rho <- if (is.finite(current_pred_dec) && current_pred_dec > 1e-15) actual_red / current_pred_dec else 0
        
        if (rho > ctrl$rho_accept && actual_red > 0) {
          g_new <- grad_func(x_try); s <- x_try - x; y <- g_new - g
          
          # Curvature Update (maintained as the Cholesky factor R, with B = R^T R)
          if (update_exact) {
            B_new <- tryCatch(hess_func(x_try), error = function(e) crossprod(R))
            B_new <- 0.5 * (B_new + t(B_new))
            R <- chol_psd(B_new)   # chol_psd safeguards positive definiteness via a shift if needed
          } else {
            Rs  <- as.numeric(R %*% s); sBs <- sum(Rs^2); sy <- sum(s * y)
            Bs  <- as.numeric(crossprod(R, Rs))                     # B s = R^T (R s)
            
            # Initial Scaling
            # Self-scale B once on the first BFGS update (scale R by sqrt of the scalar);
            # skipped when an analytic Hessian seeded B, preserving the initial curvature scale.
            if (!use_exact_hess && !scaled_B && is.finite(sy) && sy > 1e-12) {
              y_norm_sq <- sum(y * y)
              if (y_norm_sq > 1e-12) { R <- R * sqrt(y_norm_sq / sy); scaled_B <- TRUE }
              Rs <- as.numeric(R %*% s); sBs <- sum(Rs^2); Bs <- as.numeric(crossprod(R, Rs))
            }
            
            # Powell's Damping Strategy
            update_ok <- FALSE
            y_star <- y; sy_star <- sy
            if (isTRUE(ctrl$use_damped)) {
              if (is.finite(sBs) && sBs > 1e-12) {
                if (sy < ctrl$damp_phi * sBs) {
                  theta <- ((1 - ctrl$damp_phi) * sBs) / (sBs - sy)
                  y_star <- theta * y + (1 - theta) * Bs
                  sy_star <- sum(s * y_star)
                }
              }
              if (is.finite(sy_star) && sy_star > 1e-12) { y <- y_star; sy <- sy_star; update_ok <- TRUE }
            } else {
              if (is.finite(sy) && sy > 1e-12) update_ok <- TRUE
            }
            
            # Factored (dual) damped BFGS update of the Cholesky factor:
            # rank-one update by y/sqrt(s'y), then rank-one downdate by Bs/sqrt(s'Bs),
            # so that R^T R becomes B - (Bs)(Bs)'/(s'Bs) + (y y')/(s'y). A supplied analytic
            # Hessian only seeds the initial B.
            if (update_ok) {
              a_vec <- y / sqrt(sy + 1e-16)
              b_vec <- Bs / sqrt(sBs + 1e-16)
              R_upd <- chol_update(R, a_vec)
              R_dd  <- chol_downdate(R_upd, b_vec)
              if (!is.null(R_dd)) {
                R <- R_dd
              } else {
                # Numerical fallback: assemble B^+ explicitly and refactor.
                B_cur  <- crossprod(R)
                B_plus <- B_cur - (Bs %*% t(Bs)) / (sBs + 1e-16) + (y %*% t(y)) / (sy + 1e-16)
                R <- chol_psd(0.5 * (B_plus + t(B_plus)))
              }
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
          
          if (rho > ctrl$rho_expand) delta <- min(ctrl$delta_max, ctrl$delta_expand * delta)
        } else {
          # Step Rejected: Shrink Trust Region
          delta <- ctrl$delta_shrink * delta
          if (delta < 1e-14) { status <- "radius_too_small"; break }
        }
      }
    }, error = function(e) { status <<- paste0("runtime_error: ", conditionMessage(e)) })
  }
  
  # ---------- 5. Final Status & Output Construction ----------
  if (is.null(H_eval)) H_eval <- tryCatch(hess_func_pd(x), error = function(e) NULL)
  
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
    Hess_is_pd       = is_pd_fast(H_eval), 
    Hessian          = H_eval, 
    approx_hessian   = if (!is.null(R)) crossprod(R) else B,
    pred_dec         = pred_dec,
    pred_dec_avg     = pred_dec_avg
  )
}