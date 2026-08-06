#' Regime-Switching Correlation (TVTP) — DEoptim + L-BFGS-B
#'
#' Estimates a multivariate correlation model with \emph{time-varying transition probabilities} (TVTP).
#' Regime persistence is modeled via a logistic link for \eqn{N=2} or a softmax for \eqn{N \ge 3}
#' on exogenous covariates \code{X}. For \eqn{N=2}, off-diagonals equal \eqn{(1-p_{ii,t})/(N-1)};
#' for \eqn{N \ge 3}, all entries are independently determined by the softmax. Correlations are
#' regime-specific and reconstructed from lower-triangular parameters. Optimization uses global search
#' (\pkg{DEoptim}) followed by local refinement (L-BFGS-B).
#'
#' @param N Integer. Number of regimes.
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals or returns
#'   (unit variance per column is assumed).
#' @param X Numeric matrix \eqn{T \times p} of exogenous covariates driving TVTP.
#'   Include an intercept column yourself if desired (no automatic intercept).
#' @param out_of_sample Logical. If \code{TRUE}, estimation uses the first 70% of rows; the remainder is
#'   held out (indices are fixed by a 70/30 split).
#' @param control Optional list for basic verbosity control, currently supporting
#'   \code{do_trace = TRUE} to print optimizer progress and \code{seed} for
#'   reproducibility (default = 123). Algorithmic hyperparameters are set internally.
#'
#' @returns A list with:
#' \describe{
#'   \item{transition_matrix}{Representative \eqn{N \times N} transition matrix evaluated at
#'     the in-sample covariate means (\code{colMeans} of the estimation window; the first 70\%
#'     when \code{out_of_sample = TRUE}). This is a point summary at the mean covariate, \emph{not} the
#'     time-average \eqn{\frac{1}{T}\sum_t P_t}; under the nonlinear logit/softmax link the
#'     two differ (Jensen's inequality). The full per-period dynamics are governed by \code{beta}.}
#'   \item{means}{\eqn{N \times K} zeros (placeholders; means are not estimated).}
#'   \item{volatilities}{\eqn{N \times K} ones (placeholders; vols are not estimated).}
#'   \item{correlations}{\eqn{N \times C} matrix of lower-triangular correlations, \eqn{C = K(K-1)/2}.}
#'   \item{covariances}{Array \eqn{K \times K \times N} of regime correlation matrices.}
#'   \item{log_likelihood}{In-sample log-likelihood (numeric scalar); full-sample when \code{out_of_sample = FALSE}.}
#'   \item{log_likelihood_oos}{OOS log-likelihood evaluated on the held-out 30\% (numeric scalar),
#'     or \code{NULL} when \code{out_of_sample = FALSE}.}
#'   \item{beta}{\eqn{N \times p} matrix of TVTP coefficients for \eqn{N=2};
#'     \eqn{N \times (N-1)p} for \eqn{N \ge 3}, packed row-wise.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{TVTP parameterization:} for \eqn{N=2}, diagonal entries use
#'         \code{plogis(X \%*\% beta[i, ])} and off-diagonals are equal;
#'         for \eqn{N \ge 3}, a softmax over \eqn{(N-1)} free beta vectors per row,
#'         with all entries independently determined.
#'   \item \strong{Bounds:} \eqn{\beta \in [-10, 10]} elementwise; correlations \eqn{\in (-1, 1)}.
#'   \item \strong{State ordering:} regimes are reordered by ascending mean correlation for identifiability.
#'   \item \strong{Numerical safeguards:} tiny ridge is added in inverses; non-PD proposals are penalized.
#' }
#'
#' @section Assumptions:
#' Inputs in \code{residuals} are treated as mean-zero with unit variance; only the correlation structure is estimated.
#'
#' @seealso \code{\link{rsdc_estimate}}, \code{\link{rsdc_hamilton}}
#' @references
#' \insertRef{DEoptim-JSS}{RSDC}
#'
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim plogis
#' @noRd
f_optim <- function(N, residuals, X, out_of_sample = FALSE, control = list()) {
  stopifnot(is.matrix(residuals))
  if (N < 2) stop("f_optim requires N >= 2.")

  con <- list(seed = 123, do_trace = FALSE, itermax = 500, NP = NULL,
              parallelType = 0, steptol = 50, maxit = 1000, compute_se = TRUE,
              cores = 1L, start = NULL)
  con[names(control)] <- control

  # TVTP needs a regime intercept to identify a baseline persistence; the package
  # does not add one. Warn if X has no constant column.
  if (!any(apply(X, 2, function(col) isTRUE(all(abs(col - col[1]) < 1e-12)))))
    warning("X has no constant/intercept column; include one (e.g. cbind(1, X)) ",
            "so each regime has a baseline transition intercept.")

  # Parameter setup
  K <- ncol(residuals)
  p <- ncol(X)
  # N=2: one logistic beta vector per regime; N=3: (N-1) softmax beta vectors per regime row
  n_beta <- if (N == 2) N * p else N * (N - 1L) * p
  n_rho <- N * (K * (K - 1) / 2)
  rho_indices <- (n_beta + 1):(n_beta + n_rho)

  # Bounds (beta: [-10, 10] for all N; rho: (-1, 1))
  bounds <- list(
    lower = c(rep(-10, n_beta), rep(-1, n_rho)),
    upper = c(rep(10, n_beta), rep(1, n_rho))
  )

  # Out-of-sample handling (fixed indexing)
  if (out_of_sample) {
    is_cut        <- round(0.7 * nrow(residuals))
    if (is_cut < 2L || nrow(residuals) - is_cut < 1L)
      stop("Sample too small for a 70/30 out-of-sample split.")
    residuals_is  <- residuals[1:is_cut, , drop = FALSE]
    residuals_oos <- residuals[(is_cut + 1):nrow(residuals), , drop = FALSE]
    X_is  <- X[1:is_cut, , drop = FALSE]
    X_oos <- X[(is_cut + 1):nrow(X), , drop = FALSE]
    y_optim <- residuals_is
    X_optim <- X_is
  } else {
    y_optim <- residuals
    X_optim <- X
  }

  # Global search (skipped when a warm start is supplied). The RNG is seeded
  # only when the stochastic search actually runs, so warm-started refits do
  # not reset the caller's RNG stream (essential for the bootstrap, which
  # interleaves simulation and warm-started re-estimation).
  if (!is.null(con$start)) {
    if (length(con$start) != length(bounds$lower))
      stop("control$start must have length ", length(bounds$lower), " for this model.")
    de_best <- con$start
    # A warm start can be a relabeled optimum lying outside the default box:
    # after the identifiability re-ordering, re-referenced softmax betas can
    # reach [-20, 20]. L-BFGS-B silently projects an infeasible start onto the
    # box, which can quietly degrade the refit (multi-start, bootstrap); widen
    # the box just enough to keep the supplied start feasible instead.
    bounds$lower <- pmin(bounds$lower, de_best)
    bounds$upper <- pmax(bounds$upper, de_best)
  } else {
    set.seed(con$seed)
    # Reparameterized global search (1.7-0): DEoptim explores partial
    # correlations, where the whole box is feasible, instead of the natural
    # correlations, where nearly no random draw is positive definite for
    # K >~ 5. Same objective as the refinement below. The returned point is
    # already L-BFGS-B refined; the optim() call below only polishes it.
    gs <- .rsdc_theta_search(
      function(par) rsdc_likelihood(par, y = y_optim, exog = X_optim, K = K, N = N),
      method = "tvtp", K = K, N = N, p = p,
      lower = bounds$lower, upper = bounds$upper, con = con)
    de_best <- gs$par
    # A softmax-head candidate can sit outside the default box (see the
    # warm-start branch above): widen just enough to keep it feasible.
    bounds$lower <- pmin(bounds$lower, de_best)
    bounds$upper <- pmax(bounds$upper, de_best)
  }

  # Optim refinement (fixed residual reference)
  result_optim <- optim(
    par = de_best,
    fn = rsdc_likelihood,
    method = "L-BFGS-B",
    lower = bounds$lower,
    upper = bounds$upper,
    y = y_optim,
    exog = X_optim,
    K = K,
    N = N,
    control = list(maxit = con$maxit)
  )
  result_optim$convergence <- .rsdc_polish_code(
    result_optim$convergence, result_optim$value,
    rsdc_likelihood(de_best, y = y_optim, exog = X_optim, K = K, N = N))

  # Parameter processing
  optim_params <- result_optim$par
  C_per_state <- K * (K - 1) / 2

  # Generalized state ordering
  corr_means <- sapply(1:N, function(i) {
    mean(optim_params[rho_indices][(1:C_per_state) + (i - 1)*C_per_state])
  })
  sorted_states <- order(corr_means)

  # beta_dim: parameters per regime row in the flat vector
  beta_dim <- if (N == 2) p else (N - 1L) * p

  if (!all(sorted_states == 1:N)) {
    if (N == 2L) {
      beta_reordered <- unlist(lapply(sorted_states, function(s)
        optim_params[((s - 1L) * beta_dim + 1L):(s * beta_dim)]))
    } else {
      # N=3: beta columns must also be permuted and shifted for the new reference.
      # new_beta[new_row_i, new_col_j] = old_beta[s_i, s_j] - old_beta[s_i, s_N]
      # where s_i = sorted_states[i], s_j = sorted_states[j], s_N = sorted_states[N].
      ref_col <- sorted_states[N]
      beta_reordered <- unlist(lapply(sorted_states, function(s) {
        old_betas <- lapply(seq_len(N), function(col) {
          if (col < N)
            optim_params[((s - 1L) * beta_dim + (col - 1L) * p + 1L):((s - 1L) * beta_dim + col * p)]
          else rep(0, p)
        })
        ref_vec <- old_betas[[ref_col]]
        unlist(lapply(seq_len(N - 1L), function(j)
          old_betas[[sorted_states[j]]] - ref_vec))
      }))
    }
    optim_params <- c(
      beta_reordered,
      unlist(lapply(sorted_states, function(s)
        optim_params[rho_indices][((s - 1L) * C_per_state + 1L):(s * C_per_state)]))
    )
  }

  # Parameter extraction
  beta <- matrix(optim_params[1:n_beta], nrow = N, byrow = TRUE)
  rho_matrix <- matrix(optim_params[rho_indices], nrow = N, byrow = TRUE)

  # Covariance matrix reconstruction
  sigma_array <- array(dim = c(K, K, N))
  for (i in 1:N) {
    cor_mat <- diag(K)
    cor_mat[lower.tri(cor_mat)] <- rho_matrix[i, ]
    cor_mat[upper.tri(cor_mat)] <- t(cor_mat)[upper.tri(cor_mat)]
    sigma_array[,,i] <- cor_mat
  }

  # Transition matrix evaluated at the in-sample covariate means. X_optim is the
  # estimation window (X_is when out_of_sample = TRUE, full X otherwise), so this
  # summary does not leak held-out OOS covariates that beta was never fit on.
  avg_X <- colMeans(X_optim)
  P <- matrix(nrow = N, ncol = N)
  if (N == 2) {
    for (i in 1:N) {
      p_ii <- plogis(avg_X %*% beta[i, ])
      P[i, i] <- p_ii
      P[i, -i] <- (1 - p_ii) / (N - 1)
    }
  } else {
    # N=3: softmax
    for (i in 1:N) {
      logits <- c(sapply(seq_len(N - 1L), function(j)
        sum(avg_X * beta[i, ((j - 1L) * p + 1L):(j * p)])), 0)
      exp_l <- exp(logits - max(logits))
      P[i, ] <- exp_l / sum(exp_l)
    }
  }

  ll_oos <- if (out_of_sample)
    -rsdc_likelihood(optim_params, y = residuals_oos, exog = X_oos, K = K, N = N)
  else NULL

  vcov <- if (isTRUE(con$compute_se))
    .rsdc_vcov(optim_params, rsdc_likelihood, y = y_optim, exog = X_optim, K = K, N = N)
  else NULL

  list(
    transition_matrix  = P,
    means              = matrix(0, N, K),
    volatilities       = matrix(1, N, K),
    correlations       = rho_matrix,
    covariances        = sigma_array,
    log_likelihood     = -result_optim$value,
    log_likelihood_oos = ll_oos,
    beta               = beta,
    par                = optim_params,
    vcov               = vcov,
    convergence        = result_optim$convergence,
    npar               = length(optim_params),
    nobs               = nrow(y_optim),
    residuals          = y_optim,
    X                  = X_optim
  )
}

#' Regime-Switching Correlation (Fixed Transition Matrix, no X)
#'
#' Estimates a multivariate correlation model with a \emph{time-invariant} (fixed) regime
#' transition matrix. No exogenous covariates are used; both the fixed transition probabilities
#' and the regime-specific correlation parameters are estimated by global search
#' (\pkg{DEoptim}) followed by local refinement (L-BFGS-B).
#'
#' @param N Integer. Number of regimes.
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals or returns
#'   (columns are treated as mean-zero with unit variance).
#' @param out_of_sample Logical. If \code{TRUE}, estimation uses the first 70\% of rows; the
#'   remainder is held out. (Split index is a fixed 70/30 cut.)
#' @param control Optional list for basic verbosity control, currently supporting
#'   \code{do_trace = TRUE} to print optimizer progress and \code{seed} for
#'   reproducibility (default = 123). Algorithmic hyperparameters are set internally.
#'
#' @returns A list with:
#' \describe{
#'   \item{transition_matrix}{Estimated \eqn{N \times N} fixed transition matrix.}
#'   \item{means}{\eqn{N \times K} zeros (placeholders; means are not estimated).}
#'   \item{volatilities}{\eqn{N \times K} ones (placeholders; vols are not estimated).}
#'   \item{correlations}{\eqn{N \times C} matrix of lower-triangular correlations,
#'     with \eqn{C = K(K-1)/2}.}
#'   \item{covariances}{Array \eqn{K \times K \times N} of regime correlation matrices.}
#'   \item{log_likelihood}{In-sample log-likelihood (numeric scalar); full-sample when \code{out_of_sample = FALSE}.}
#'   \item{log_likelihood_oos}{OOS log-likelihood evaluated on the held-out 30\% (numeric scalar),
#'     or \code{NULL} when \code{out_of_sample = FALSE}.}
#'   \item{beta}{\code{NULL} (no covariates in this specification).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Parameterization:} For \eqn{N=2}, the transition parameters are
#'         \eqn{\{p_{11}, p_{22}\}}; off-diagonals are \eqn{1 - p_{ii}}. For
#'         \eqn{N \ge 3}, each row carries \eqn{N-1} free probabilities (packed
#'         row-wise) and the last column is the row complement; infeasible rows
#'         (free probabilities summing above 1) are penalized during optimization
#'         and, if present at the optimum, projected onto the simplex with a warning.
#'   \item \strong{Bounds:} \eqn{p_{ii} \in [0.01, 0.99]} and \eqn{\rho \in (-1, 1)}
#'     for numerical stability and identifiability.
#'   \item \strong{Correlation build:} Each regime's correlation matrix is reconstructed from
#'     its lower-triangular vector and symmetrized; non-PD proposals are penalized in the
#'     likelihood routine.
#'   \item \strong{State ordering:} regimes are reordered by ascending mean correlation for identifiability.
#' }
#'
#' @section Assumptions:
#' Inputs in \code{residuals} are treated as mean-zero with unit variance; only the correlation
#' structure and fixed transition probabilities are estimated.
#'
#' @seealso \code{\link{rsdc_estimate}} (wrapper), \code{\link{rsdc_hamilton}}, \code{\link{rsdc_likelihood}}
#'
#' @references
#' \insertRef{DEoptim-JSS}{RSDC}
#'
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim
#' @noRd
f_optim_noX <- function(N, residuals, out_of_sample = FALSE, control = list()) {
  stopifnot(is.matrix(residuals))
  if (N < 2) stop("f_optim_noX requires N >= 2.")

  con <- list(seed = 123, do_trace = FALSE, itermax = 500, NP = NULL,
              parallelType = 0, steptol = 50, maxit = 1000, compute_se = TRUE,
              cores = 1L, start = NULL)
  con[names(control)] <- control

  K <- ncol(residuals)
  n_obs <- nrow(residuals)
  n_trans <- N * (N - 1)
  C_per_state <- K * (K - 1) / 2
  n_rho <- N * C_per_state

  # Each free transition prob in [0.01, 0.99] for all N. For N=3 the implied third
  # prob (1 - p_i1 - p_i2) can be negative for some proposals; rsdc_likelihood rejects
  # those rows (returns 1e10). This keeps persistent diagonals (p_ii up to 0.99)
  # reachable instead of capping them at 0.485.
  trans_upper <- 0.99
  bounds <- list(
    lower = c(rep(0.01, n_trans), rep(-1, n_rho)),
    upper = c(rep(trans_upper, n_trans), rep(1, n_rho))
  )

  if (out_of_sample) {
    cut     <- round(0.7 * n_obs)
    if (cut < 2L || n_obs - cut < 1L)
      stop("Sample too small for a 70/30 out-of-sample split.")
    y_optim <- residuals[1:cut, , drop = FALSE]
    y_oos   <- residuals[(cut + 1):n_obs, , drop = FALSE]
  } else {
    y_optim <- residuals
    y_oos   <- NULL
  }

  # RNG seeded only when the stochastic global search runs (see f_optim).
  if (!is.null(con$start)) {
    if (length(con$start) != length(bounds$lower))
      stop("control$start must have length ", length(bounds$lower), " for this model.")
    de_best <- con$start
    # A warm start can be a relabeled optimum lying outside the default box:
    # after the identifiability re-ordering (N >= 3), a former row-complement
    # transition entry can fall below 0.01 (or above 0.99) while remaining a
    # valid probability. L-BFGS-B silently projects an infeasible start onto
    # the box, which can quietly degrade the refit (multi-start, bootstrap);
    # widen the box just enough to keep the supplied start feasible instead.
    bounds$lower <- pmin(bounds$lower, de_best)
    bounds$upper <- pmax(bounds$upper, de_best)
  } else {
    set.seed(con$seed)
    # Reparameterized global search (1.7-0): partial-correlation blocks and,
    # for N >= 3, a bounded-softmax transition head (the [0.01, 0.99] box on
    # the log-odds scale), so every point of the search box is a valid,
    # non-degenerate model. Same objective as the refinement below.
    gs <- .rsdc_theta_search(
      function(par) rsdc_likelihood(par, y = y_optim, exog = NULL, K = K, N = N),
      method = "noX", K = K, N = N,
      lower = bounds$lower, upper = bounds$upper, con = con)
    de_best <- gs$par
    # A softmax-head candidate can hold transition probabilities below the
    # 0.01 box edge: widen just enough, as in the warm-start branch above.
    bounds$lower <- pmin(bounds$lower, de_best)
    bounds$upper <- pmax(bounds$upper, de_best)
  }

  optim_result <- optim(
    par = de_best,
    fn = rsdc_likelihood,
    method = "L-BFGS-B",
    lower = bounds$lower,
    upper = bounds$upper,
    y = y_optim, exog = NULL, K = K, N = N,
    control = list(maxit = con$maxit)
  )
  optim_result$convergence <- .rsdc_polish_code(
    optim_result$convergence, optim_result$value,
    rsdc_likelihood(de_best, y = y_optim, exog = NULL, K = K, N = N))

  final_par <- optim_result$par

  # Regime reordering for identifiability (ascending mean correlation)
  trans_per_regime <- N - 1L
  rho_params <- final_par[(n_trans + 1):(n_trans + n_rho)]
  corr_means <- sapply(1:N, function(i) mean(rho_params[(1:C_per_state) + (i - 1L) * C_per_state]))
  sorted_states <- order(corr_means)
  if (!all(sorted_states == 1:N)) {
    orig_trans <- final_par[1:n_trans]
    if (N == 2L) {
      trans_reordered <- unlist(lapply(sorted_states, function(s)
        orig_trans[((s - 1L) * trans_per_regime + 1L):(s * trans_per_regime)]))
    } else {
      # N=3: reconstruct full P row, then select columns in new state order
      trans_reordered <- unlist(lapply(sorted_states, function(s) {
        free     <- orig_trans[((s - 1L) * trans_per_regime + 1L):(s * trans_per_regime)]
        full_row <- c(free, 1 - sum(free))
        full_row[sorted_states][seq_len(N - 1L)]
      }))
    }
    final_par <- c(
      trans_reordered,
      unlist(lapply(sorted_states, function(s)
        rho_params[((s - 1L) * C_per_state + 1L):(s * C_per_state)]))
    )
  }

  rho_matrix <- matrix(final_par[(n_trans + 1):(n_trans + n_rho)], nrow = N, byrow = TRUE)

  sigma_array <- array(dim = c(K, K, N))
  for (i in 1:N) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_matrix[i, ]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    sigma_array[,,i] <- R
  }

  # Build transition matrix. N=2 keeps the diagonal-stay parameterization
  # (param i = p_ii); for N >= 3 the free entries fill the first N-1 columns and
  # the last column is the row complement (generic for any N).
  trans_params <- final_par[1:n_trans]
  P <- matrix(0, N, N)
  infeasible_row <- FALSE
  if (N == 2) {
    P[1, 1] <- trans_params[1]
    P[1, 2] <- 1 - trans_params[1]
    P[2, 1] <- 1 - trans_params[2]
    P[2, 2] <- trans_params[2]
  } else {
    for (i in 1:N) {
      free <- trans_params[((i - 1L) * (N - 1L) + 1L):(i * (N - 1L))]
      row  <- c(free, 1 - sum(free))
      # A row can be infeasible (any entry < 0) when the optimizer has not reached
      # a feasible simplex, or after the identifiability re-ordering re-references
      # the free entries. Project onto the simplex (clamp negatives, renormalize)
      # so a valid stochastic matrix is always returned, and flag non-convergence.
      if (any(row < 0)) {
        infeasible_row <- TRUE
        row <- pmax(row, 0)
        row <- row / sum(row)
      }
      P[i, ] <- row
    }
  }
  if (infeasible_row) {
    # Re-write the transition block of the parameter vector from the projected
    # (valid) matrix so that par, log-likelihood, vcov and the OOS likelihood
    # below are all consistent with the returned transition_matrix.
    final_par[1:n_trans] <- as.vector(t(P[, seq_len(N - 1L), drop = FALSE]))
    optim_result$value <- rsdc_likelihood(final_par, y = y_optim, exog = NULL, K = K, N = N)
    # The projected vector is not a true optimum under the simplex constraints, so flag
    # non-convergence and suppress observed-information SEs (the Hessian at a projected,
    # non-stationary point is not a valid basis for inference).
    optim_result$convergence <- 99L
    warning("Transition matrix had an infeasible row (free probabilities summed > 1) ",
            "and was projected onto the simplex; the optimizer did not converge to a ",
            "feasible optimum. Standard errors are suppressed (convergence flagged). ",
            "Increase control$itermax or control$NP.")
  }

  # Warn when a free transition probability is pinned to its [0.01, 0.99] bound
  # (degenerate / near-absorbing regime; SEs there are unreliable).
  if (any(trans_params <= 0.01 + 1e-4 | trans_params >= trans_upper - 1e-4))
    warning("A transition probability is at the [0.01, 0.99] bound; the regime is ",
            "near-degenerate and its standard error may be unreliable.")

  ll_oos <- if (!is.null(y_oos))
    -rsdc_likelihood(final_par, y = y_oos, exog = NULL, K = K, N = N)
  else NULL

  vcov <- if (isTRUE(con$compute_se) && !infeasible_row)
    .rsdc_vcov(final_par, rsdc_likelihood, y = y_optim, exog = NULL, K = K, N = N)
  else NULL

  list(
    transition_matrix  = P,
    means              = matrix(0, N, K),
    volatilities       = matrix(1, N, K),
    correlations       = rho_matrix,
    covariances        = sigma_array,
    log_likelihood     = -optim_result$value,
    log_likelihood_oos = ll_oos,
    beta               = NULL,
    par                = final_par,
    vcov               = vcov,
    convergence        = optim_result$convergence,
    npar               = length(final_par),
    nobs               = nrow(y_optim),
    residuals          = y_optim,
    X                  = NULL
  )
}

#' Constant-Correlation Model (Single Regime)
#'
#' Estimates a single (time-invariant) correlation matrix assuming asset correlations
#' are constant over time. No regime switching or exogenous covariates are used.
#' Optimization uses a global stage (\pkg{DEoptim}) followed by local refinement
#' (L-BFGS-B via \code{stats::optim}).
#'
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals or
#'   returns (columns treated as mean-zero with unit variance).
#' @param out_of_sample Logical. If \code{TRUE}, estimation uses the first 70\% of rows; the
#'   remainder is held out. (Split index is a fixed 70/30 cut.)
#' @param control Optional list for basic verbosity control, currently supporting
#'   \code{do_trace = TRUE} to print optimizer progress and \code{seed} for
#'   reproducibility (default = 123). Algorithmic hyperparameters are set internally.
#'
#' @returns A list with:
#' \describe{
#'   \item{transition_matrix}{\eqn{1 \times 1} matrix equal to 1 (no switching).}
#'   \item{means}{\eqn{1 \times K} zeros (placeholders; means are not estimated).}
#'   \item{volatilities}{\eqn{1 \times K} ones (placeholders; vols are not estimated).}
#'   \item{correlations}{\eqn{1 \times C} vector of lower-triangular correlations,
#'     where \eqn{C = K(K-1)/2}.}
#'   \item{covariances}{Array \eqn{K \times K \times 1} with the full correlation matrix.}
#'   \item{log_likelihood}{In-sample log-likelihood (numeric scalar); full-sample when \code{out_of_sample = FALSE}.}
#'   \item{log_likelihood_oos}{OOS log-likelihood evaluated on the held-out 30\% (numeric scalar),
#'     or \code{NULL} when \code{out_of_sample = FALSE}.}
#'   \item{beta}{\code{NULL} (no covariates).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Parameterization:} the free parameters are the \eqn{C = K(K-1)/2}
#'         pairwise correlations stacked in \code{lower.tri} order.
#'   \item \strong{Bounds / PD checks:} elementwise bounds \eqn{\rho \in (-1, 1)};
#'         non–positive-definite proposals are penalized in the objective.
#'   \item \strong{Numerical safeguards:} a tiny ridge is added before inverting
#'         the correlation matrix to improve stability.
#'   \item \strong{Scaling:} since inputs are treated as unit-variance, only the
#'         correlation structure is estimated.
#' }
#'
#' @seealso \code{\link{rsdc_estimate}} (wrapper).
#'
#' @references
#' \insertRef{DEoptim-JSS}{RSDC}
#'
#' \insertRef{hamilton1989}{RSDC}
#'
#' \insertRef{pelletier2006regime}{RSDC}
#'
#' @note Uses \code{DEoptim::DEoptim()} then \code{stats::optim(method = "L-BFGS-B")}.
#'
#' @importFrom stats optim
#' @noRd
f_optim_const <- function(residuals, out_of_sample = FALSE, control = list()) {
  stopifnot(is.matrix(residuals))

  con <- list(seed = 123, do_trace = FALSE, itermax = 500, NP = NULL,
              steptol = 50, maxit = 1000, compute_se = TRUE, start = NULL)
  con[names(control)] <- control

  K <- ncol(residuals)
  n_rho <- K * (K - 1) / 2

  neg_loglik_const <- function(rho_vec, y, K) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_vec
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    # PD check
    eig <- tryCatch(eigen(R, only.values = TRUE)$values, error = function(e) NA_real_)
    if (any(!is.finite(eig)) || min(eig) < 1e-8) return(1e10)

    R_rdg  <- R + diag(1e-8, K)            # ridge used consistently for inverse and log-det
    inv_R  <- solve(R_rdg)
    logdet <- determinant(R_rdg, logarithm = TRUE)$modulus[1]
    quad   <- rowSums((y %*% inv_R) * y)
    # negative log-likelihood
    return(-(-0.5 * sum(logdet + quad + K * log(2 * pi))))
  }

  # Split (if requested)
  if (out_of_sample) {
    cut <- round(0.7 * nrow(residuals))
    if (cut < 2L || nrow(residuals) - cut < 1L)
      stop("Sample too small for a 70/30 out-of-sample split.")
    y_is  <- residuals[1:cut, , drop = FALSE]
    y_oos <- residuals[(cut + 1):nrow(residuals), , drop = FALSE]
  } else {
    y_is  <- residuals
    y_oos <- NULL
  }

  # Fit on IS only. RNG seeded only when the stochastic global search runs
  # (see f_optim).
  if (!is.null(con$start)) {
    if (length(con$start) != n_rho)
      stop("control$start must have length ", n_rho, " for the constant model.")
    de_best <- con$start
  } else {
    set.seed(con$seed)
    # Reparameterized global search (1.7-0): the partial-correlation box is
    # fully feasible, so the search needs no penalty rejections even at
    # large K. Same objective as the refinement below.
    gs <- .rsdc_theta_search(
      function(par) neg_loglik_const(par, y = y_is, K = K),
      method = "const", K = K, N = 1L,
      lower = rep(-1, n_rho), upper = rep(1, n_rho),
      con = list(itermax = con$itermax, NP = con$NP, steptol = con$steptol,
                 maxit = con$maxit, do_trace = con$do_trace))
    de_best <- gs$par
  }

  optim_result <- optim(
    par = de_best,
    fn = neg_loglik_const,
    method = "L-BFGS-B",
    lower = rep(-1, n_rho),
    upper = rep(1, n_rho),
    y = y_is, K = K,
    control = list(maxit = con$maxit)
  )
  optim_result$convergence <- .rsdc_polish_code(
    optim_result$convergence, optim_result$value,
    neg_loglik_const(de_best, y = y_is, K = K))

  rho_vec <- optim_result$par
  R <- diag(K)
  R[lower.tri(R)] <- rho_vec
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  array_R <- array(R, dim = c(K, K, 1))

  # IS log-likelihood (full-sample when out_of_sample = FALSE)
  ll_return <- -optim_result$value

  # OOS log-likelihood: evaluate IS-fitted params on held-out data
  ll_oos <- if (!is.null(y_oos))
    -neg_loglik_const(rho_vec, y = y_oos, K = K)
  else NULL

  vcov <- if (isTRUE(con$compute_se))
    .rsdc_vcov(rho_vec, neg_loglik_const, y = y_is, K = K)
  else NULL

  list(
    transition_matrix  = matrix(1, 1, 1),
    means              = matrix(0, 1, K),
    volatilities       = matrix(1, 1, K),
    correlations       = matrix(rho_vec, nrow = 1),
    covariances        = array_R,
    log_likelihood     = ll_return,
    log_likelihood_oos = ll_oos,
    beta               = NULL,
    par                = rho_vec,
    vcov               = vcov,
    convergence        = optim_result$convergence,
    npar               = length(rho_vec),
    nobs               = nrow(y_is),
    residuals          = y_is,
    X                  = NULL
  )
}

#' Estimate Regime-Switching or Constant Correlation Model (Wrapper)
#'
#' Unified front-end that dispatches to one of three estimators:
#' \itemize{
#'   \item \code{f_optim()} — TVTP specification (\code{method = "tvtp"}).
#'   \item \code{f_optim_noX()} — fixed transition matrix (\code{method = "noX"}).
#'   \item \code{f_optim_const()} — constant correlation, single regime (\code{method = "const"}).
#' }
#'
#' @param method Character. One of \code{"tvtp"}, \code{"noX"}, \code{"const"}.
#' @param residuals Numeric matrix \eqn{T \times K}. Typically standardized residuals/returns.
#' @param N Integer. Number of regimes. Ignored when \code{method = "const"}.
#' @param X Numeric matrix \eqn{T \times p} of exogenous covariates (required for \code{"tvtp"}).
#' @param out_of_sample Logical. If \code{TRUE}, a fixed 70/30 split is applied prior to estimation.
#' @param control Optional named list of estimation settings, forwarded to the
#'   backends and optimizers. Any subset of the following elements may be
#'   supplied (defaults in parentheses):
#'   \describe{
#'     \item{\code{start} (\code{NULL})}{Warm start; skips the \pkg{DEoptim}
#'       global search. Two forms are accepted. A \emph{numeric vector} in the
#'       objective's packing order (see \code{\link{rsdc_likelihood}}) refines
#'       a single point locally — e.g. re-examine an earlier optimum with
#'       \code{control = list(start = coef(fit0))}. An
#'       \code{\link{rsdc_starts}} \emph{object} runs a data-driven
#'       multi-start instead: each start is refined locally, the
#'       highest-likelihood fit is kept, and the fit carries
#'       \code{start_logliks} and \code{start_pars}. Since 1.7-0 the default
#'       global search is feasible at any \eqn{K} (partial-correlation
#'       parameterization); \code{rsdc_starts} remains a fast deterministic
#'       complement — when both routes agree, the optimum is cross-certified:
#'       \code{st <- rsdc_starts(y, N = 2, method = "noX")} then
#'       \code{control = list(start = st, cores = 4)}. A warm start disables
#'       \code{n_starts}.}
#'     \item{\code{n_starts} (\code{1})}{Number of independent global+local
#'       searches from seeds \code{seed, seed + 1, \dots}, keeping the
#'       highest-likelihood fit. Because \pkg{DEoptim} is a \emph{stochastic}
#'       global search, a single run can land in a suboptimal basin of a
#'       multimodal likelihood; independent searches make that failure both
#'       unlikely and \emph{visible}: with \code{n_starts >= 2} the fit
#'       carries \code{start_logliks}, and a warning is issued when the best
#'       value is not replicated by at least two searches.
#'       \code{control = list(n_starts = 4, cores = 4)} is a good default for
#'       multimodal surfaces (e.g. \code{"noX"} with \eqn{N \ge 3}) at the
#'       wall-clock cost of a single search.}
#'     \item{\code{cores} (\code{1})}{Parallel workers, via the base
#'       \pkg{parallel} package (forked on Unix, PSOCK cluster on Windows).
#'       With either multi-start form, the starts are fitted in parallel —
#'       \code{control = list(n_starts = 4, cores = 4)} — and the result is
#'       identical for any value, because each fit is deterministic given its
#'       start or seed; only the wall time changes.}
#'     \item{\code{seed} (\code{123})}{RNG seed of the stochastic \pkg{DEoptim}
#'       search: \code{control = list(seed = 42)}. Unused when a warm start
#'       skips the global search.}
#'     \item{\code{compute_se} (\code{TRUE})}{Compute observed-information
#'       standard errors at the optimum. Disable for speed inside loops
#'       (bootstraps, searches): \code{control = list(compute_se = FALSE)}.}
#'     \item{\code{itermax} (\code{500}), \code{NP} (\code{10 x} npar),
#'       \code{steptol} (\code{50})}{Budget of the \pkg{DEoptim} global
#'       search: maximum generations, population size, early stop after
#'       \code{steptol} stagnant generations. E.g. a quick exploratory fit:
#'       \code{control = list(itermax = 100, NP = 200)}.}
#'     \item{\code{maxit} (\code{1000})}{Maximum iterations of the local
#'       L-BFGS-B refinement (\code{\link[stats]{optim}}).}
#'     \item{\code{do_trace} (\code{FALSE})}{Print optimizer progress.}
#'   }
#'
#' @return An object of class \code{"rsdc_fit"}: a list with components
#' \describe{
#'   \item{\code{transition_matrix}}{Estimated transition matrix (\eqn{1 \times 1} for \code{"const"};
#'     for \code{"tvtp"}, evaluated at the in-sample covariate means).}
#'   \item{\code{correlations}}{Regime lower-triangular correlations (rows ordered by ascending mean).}
#'   \item{\code{covariances}}{Array \eqn{K \times K \times N} of regime correlation matrices.}
#'   \item{\code{log_likelihood}}{In-sample log-likelihood; full-sample when \code{out_of_sample = FALSE}.}
#'   \item{\code{log_likelihood_oos}}{OOS log-likelihood on held-out 30 percent, or \code{NULL}.}
#'   \item{\code{beta}}{TVTP coefficients (only for \code{"tvtp"}; \code{NULL} otherwise).}
#'   \item{\code{par}, \code{coefficients}}{Named vector of estimated parameters in objective order.}
#'   \item{\code{vcov}}{Observed-information variance-covariance at the MLE, or \code{NULL} if
#'     the Hessian is singular/unavailable (e.g. at a bound).}
#'   \item{\code{se}}{Named standard errors (\code{sqrt(diag(vcov))}), or \code{NULL}.}
#'   \item{\code{convergence}}{\code{optim} convergence code (0 = converged).}
#'   \item{\code{npar}, \code{nobs}}{Number of free parameters and observations used in estimation.}
#'   \item{\code{method}, \code{N}, \code{K}, \code{p}, \code{call}}{Fit metadata.}
#' }
#' The returned object supports the standard fitted-model generics —
#' \code{print}, \code{summary}, \code{coef}, \code{logLik} (so
#' \code{\link[stats]{AIC}}/\code{\link[stats]{BIC}} work), \code{nobs},
#' \code{vcov}, \code{confint}, \code{predict}, \code{simulate}, and
#' \code{plot} — documented together in \link[=rsdc_fit-methods]{Methods for
#' fitted RSDC models}.
#'
#' @details
#' \itemize{
#'   \item \strong{Method selection:} \code{match.arg()} validates \code{method}.
#'   \item \strong{Inputs:} \code{"tvtp"} requires non-NULL \code{X}; \code{N} is ignored for \code{"const"}.
#'   \item \strong{Split:} If \code{out_of_sample = TRUE}, the first 70% is used for fitting.
#' }
#'
#' @references
#' \insertRef{DEoptim-JSS}{RSDC} \cr
#'
#' \insertRef{hamilton1989}{RSDC} \cr
#'
#' \insertRef{pelletier2006regime}{RSDC}
#'
#' @examples
#' \donttest{
#' # Two persistent regimes: low (0.1) vs high (0.8) correlation
#' sim <- rsdc_simulate(n = 500, X = matrix(1, 500, 1),
#'                      beta = matrix(qlogis(0.9), 2, 1),
#'                      mu = matrix(0, 2, 2),
#'                      sigma = array(c(1, 0.1, 0.1, 1,
#'                                      1, 0.8, 0.8, 1), c(2, 2, 2)),
#'                      N = 2, seed = 2)
#' fit <- rsdc_estimate("noX", residuals = sim$observations, N = 2)
#' print(fit)
#' summary(fit)
#' c(AIC = AIC(fit), BIC = BIC(fit))
#' X <- cbind(1, scale(seq_len(nrow(sim$observations))))
#' fit_tvtp <- rsdc_estimate("tvtp", residuals = sim$observations, N = 2, X = X)
#' coef(fit_tvtp)
#' }
#'
#' @seealso \link[=rsdc_fit-methods]{Methods for fitted RSDC models},
#'   \code{\link{rsdc_hamilton}} and \code{\link{rsdc_likelihood}}.
#'
#' @export
rsdc_estimate <- function(method = c("tvtp", "noX", "const"),
                          residuals, N = 2, X = NULL,
                          out_of_sample = FALSE, control = list()) {
  method <- match.arg(method)
  # N must be a whole number: a fractional N propagates into the parameter
  # counts, every candidate then fails rsdc_likelihood's length check, and the
  # fit is returned with the infeasibility penalty reported as a log-likelihood.
  if (!is.null(N) && (length(N) != 1L || !is.finite(N) || N != as.integer(N)))
    stop("N must be a single whole number (got ", paste(format(N), collapse = ", "), ").")
  if (!is.null(N)) N <- as.integer(N)
  if (method != "const" && !is.null(N) && N < 2)
    stop("N must be at least 2. Use method='const' for a single-regime model.")

  if (!is.matrix(residuals)) stop("residuals must be a numeric matrix.")
  if (ncol(residuals) < 2L)
    stop("residuals must have at least 2 columns (K >= 2); the model estimates a ",
         "multivariate correlation structure.")
  if (any(!is.finite(residuals))) stop("residuals must not contain NA/NaN/Inf values.")
  if (!is.null(X) && any(!is.finite(as.matrix(X))))
    stop("X must not contain NA/NaN/Inf values.")
  if (!is.null(X) && nrow(as.matrix(X)) != nrow(residuals))
    stop(sprintf("X must have the same number of rows as residuals (nrow(X) = %d, nrow(residuals) = %d).",
                 nrow(as.matrix(X)), nrow(residuals)))

  if (nrow(residuals) < 2L)
    stop("residuals must have at least 2 rows (nrow(residuals) = ",
         nrow(residuals), "); a correlation model cannot be identified from ",
         "a single observation.")

  # Correlation-only model: columns are assumed mean-zero, unit-variance. Warn (do not
  # auto-transform) when inputs deviate materially so raw returns are not silently misspecified.
  if (is.matrix(residuals)) {
    col_sd <- apply(residuals, 2, stats::sd)
    if (any(is.finite(col_sd) & (col_sd < 0.5 | col_sd > 2)))
      warning("residuals columns are not ~unit-variance (sd outside [0.5, 2]); the model treats ",
              "inputs as standardized. Consider scale() so only the correlation structure is fit.")
  }

  cl <- match.call()
  p_x <- if (!is.null(X)) ncol(X) else NA_integer_

  if (method == "tvtp" && is.null(X)) stop("X must be provided for method = 'tvtp'")

  # Single fit at a given seed (compute_se controllable for the multi-start search).
  fit_once <- function(ctrl) {
    if (method == "tvtp")
      f_optim(N = N, residuals = residuals, X = X,
              out_of_sample = out_of_sample, control = ctrl)
    else if (method == "noX")
      f_optim_noX(N = N, residuals = residuals,
                  out_of_sample = out_of_sample, control = ctrl)
    else if (method == "const")
      f_optim_const(residuals = residuals, out_of_sample = out_of_sample,
                    control = ctrl)
    else stop("Unknown method: ", method)
  }

  # Multi-start: with control$n_starts > 1, re-run the global+local search from
  # several seeds, keep the highest-likelihood fit, and record the spread of
  # log-likelihoods so users can judge whether the optimum is stable (a guard
  # against local optima). A warm start (control$start) disables multi-start.
  # Alternatively, control$start may be an "rsdc_starts" object (data-driven
  # warm starts for high-dimensional problems, see rsdc_starts()): each start
  # is refined by the local optimizer (the global search is skipped) and the
  # highest-likelihood fit is kept. control$cores parallelizes whichever
  # multi-start form is active; the backends never see it there (no nested
  # parallelism). Since 1.7-0 a plain single fit ignores cores: parallelism
  # lives across the independent searches of the multi-start forms, not
  # inside DEoptim (see NEWS).
  n_starts <- if (!is.null(control$n_starts)) as.integer(control$n_starts) else 1L
  base_seed <- if (!is.null(control$seed)) control$seed else 123L
  cores <- if (!is.null(control$cores)) as.integer(control$cores) else 1L
  start_logliks <- NULL; start_pars <- NULL
  if (inherits(control$start, "rsdc_starts")) {
    st <- control$start
    if (!identical(st$method, method))
      stop("The rsdc_starts object was built for method = '", st$method,
           "', not '", method, "'.")
    if (method != "const" && st$N != N)
      stop("The rsdc_starts object was built for N = ", st$N, ", not N = ", N, ".")
    if (st$K != ncol(residuals))
      stop("The rsdc_starts object was built for K = ", st$K, " series, not ",
           ncol(residuals), ".")
    if (method == "tvtp" && !is.na(st$p) && st$p != ncol(X))
      stop("The rsdc_starts object was built for p = ", st$p,
           " covariates, not ", ncol(X), ".")
    keep_starts <- which(is.finite(st$loglik0))
    if (!length(keep_starts))
      stop("No feasible start in the rsdc_starts object (all loglik0 infinite).")
    search_ctrl <- control; search_ctrl$compute_se <- FALSE
    search_ctrl$n_starts <- NULL; search_ctrl$cores <- NULL
    fits <- .rsdc_lapply(keep_starts, function(i) {
      sc <- search_ctrl; sc$start <- st$starts[i, ]
      tryCatch(fit_once(sc), error = function(e) NULL)
    }, cores = cores)
    ok <- !vapply(fits, is.null, logical(1))
    if (!any(ok)) stop("All warm starts failed to converge.")
    fits <- fits[ok]
    start_logliks <- vapply(fits, function(f) f$log_likelihood, numeric(1))
    start_pars <- t(vapply(fits, function(f) f$par,
                           numeric(length(fits[[1]]$par))))
    best <- which.max(start_logliks)
    # Recompute the winning fit with the user's compute_se setting,
    # warm-starting the local optimizer at the winner's parameters.
    final_ctrl <- control; final_ctrl$n_starts <- NULL; final_ctrl$cores <- NULL
    final_ctrl$start <- fits[[best]]$par
    fit <- fit_once(final_ctrl)
    n_starts <- length(start_logliks)
  } else if (n_starts > 1L && is.null(control$start)) {
    seeds <- base_seed + seq_len(n_starts) - 1L
    search_ctrl <- control; search_ctrl$compute_se <- FALSE
    search_ctrl$n_starts <- NULL; search_ctrl$cores <- NULL
    fits <- .rsdc_lapply(seeds, function(s) {
      search_ctrl$seed <- s; fit_once(search_ctrl)
    }, cores = cores)
    start_logliks <- vapply(fits, function(f) f$log_likelihood, numeric(1))
    best <- which.max(start_logliks)
    # Replication (agreement) check: at a reproducible optimum, several
    # independent searches should reach the same maximum. A lone maximum on a
    # multimodal surface is indistinguishable from a lucky draw.
    if (sum(start_logliks >= start_logliks[best] - .RSDC_AGREE_TOL) < 2L)
      warning("The best log-likelihood was reached by only one of ", n_starts,
              " independent searches; the optimum may not be global. ",
              "Inspect fit$start_logliks and consider increasing ",
              "control$n_starts.")
    # Recompute the winning fit with the user's compute_se setting, warm-starting the
    # local optimizer at the winner's parameters (avoids repeating the global search).
    final_ctrl <- control; final_ctrl$n_starts <- NULL; final_ctrl$cores <- NULL
    final_ctrl$seed <- seeds[best]; final_ctrl$start <- fits[[best]]$par
    fit <- fit_once(final_ctrl)
  } else {
    fc <- control; fc$n_starts <- NULL
    fit <- fit_once(fc)   # single fit: cores is ignored since 1.7-0 (no DEoptim parallelType)
  }

  out <- .rsdc_finalize(fit, method = method, N = N, K = ncol(residuals), p = p_x, call = cl,
                        se_requested = is.null(control$compute_se) ||
                                       isTRUE(control$compute_se))
  if (!is.null(start_logliks)) {
    out$n_starts      <- n_starts
    out$start_logliks <- start_logliks
  }
  if (!is.null(start_pars)) {
    colnames(start_pars) <- names(out$coefficients)
    out$start_pars <- start_pars
  }
  out
}
