#' Hamilton Filter (Fixed P or TVTP)
#'
#' Runs the Hamilton (1989) filter for a multivariate regime-switching \emph{correlation} model.
#' Supports either a fixed (time-invariant) transition matrix \eqn{P} or time-varying transition
#' probabilities (TVTP) built from exogenous covariates \code{X}: a logistic link for \eqn{N=2}
#' or a softmax link for \eqn{N \ge 3}.
#' Returns filtered/smoothed regime probabilities and the log-likelihood.
#'
#' @param y Numeric matrix \eqn{T \times K} of observations (e.g., standardized residuals/returns).
#'   Columns are treated as mean-zero with unit variance; only the correlation structure is modeled.
#' @param X Optional numeric matrix \eqn{T \times p} of covariates for TVTP. Required if \code{beta} is supplied.
#' @param beta Optional numeric matrix of TVTP coefficients.
#'   For \eqn{N=2}: an \eqn{N \times p} matrix; row \eqn{i} gives the logistic coefficient
#'   vector so that \eqn{p_{ii,t} = \mathrm{logit}^{-1}(X_t^\top \beta_i)}.
#'   For \eqn{N \ge 3}: an \eqn{N \times (N-1)p} matrix packed row-wise; row \eqn{i}
#'   contains \eqn{N-1} consecutive length-\eqn{p} softmax vectors; the \eqn{N}-th
#'   logit is fixed at 0 (reference category).
#' @param rho_matrix Numeric matrix \eqn{N \times C} of regime correlation parameters, where
#'   \eqn{C = K(K-1)/2}. Each row is the lower-triangular part (by \code{lower.tri}) of a regime's
#'   correlation matrix.
#' @param K Integer. Number of observed series (columns of \code{y}).
#' @param N Integer. Number of regimes.
#' @param P Optional \eqn{N \times N} fixed transition matrix. Used only when \code{X} or \code{beta} is \code{NULL}.
#' @param xi_init Optional numeric vector of length \code{N}. If supplied, the
#'   forward filter is initialized from this vector (normalized internally to sum
#'   to 1) instead of the default uniform \eqn{1/N}. Intended for two-pass
#'   out-of-sample forecasting: run the filter on the in-sample period first,
#'   extract the terminal column of \code{filtered_probs}, and pass it here to
#'   initialize the out-of-sample filter run. Must be non-negative and finite.
#' @param engine Character; \code{"cpp"} (default) runs the filter/smoother in C++
#'   (\pkg{RcppArmadillo}), \code{"r"} the pure-R reference. Both give identical
#'   results (used for equivalence testing).
#'
#' @returns A list with:
#' \describe{
#'   \item{filtered_probs}{\eqn{N \times T} matrix of filtered probabilities
#'     \eqn{\Pr(S_t = j \mid \Omega_t)}.}
#'   \item{smoothed_probs}{\eqn{N \times T} matrix of smoothed probabilities
#'     \eqn{\Pr(S_t = j \mid \Omega_T)}.}
#'   \item{log_likelihood}{Scalar log-likelihood of the model given \code{y}.}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Correlation rebuild:} For regime \eqn{m}, a correlation matrix \eqn{R_m} is reconstructed
#'         from \code{rho_matrix[m, ]} (lower-triangular fill + symmetrization). Non-PD proposals are penalized.
#'   \item \strong{Transition dynamics:}
#'         \itemize{
#'           \item \emph{Fixed P:} If \code{X} or \code{beta} is missing, a constant \eqn{P} is used
#'                 (user-provided via \code{P}; otherwise uniform \eqn{1/N} rows).
#'           \item \emph{TVTP:} With \code{X} and \code{beta}:
#'                 for \eqn{N=2}, \eqn{p_{ii,t} = \mathrm{logit}^{-1}(X_t^\top \beta_i)} and
#'                 \eqn{p_{ij,t} = 1 - p_{ii,t}} (\eqn{j \ne i});
#'                 for \eqn{N \ge 3}, a softmax over \eqn{N-1} free logit vectors per row with
#'                 the \eqn{N}-th logit fixed at 0 (reference category); log-sum-exp stabilized.
#'         }
#'   \item \strong{Numerical safeguards:} A small ridge is added before inversion; if filtering
#'         degenerates at a time step, \code{log_likelihood = -Inf} is returned.
#' }
#'
#' @examples
#' set.seed(1)
#' T <- 50; K <- 3; N <- 2
#' y <- scale(matrix(rnorm(T * K), T, K), center = TRUE, scale = TRUE)
#'
#' # Example rho: two regimes with different average correlations
#' rho <- rbind(c(0.10, 0.05, 0.00),
#'              c(0.60, 0.40, 0.30))  # lower-tri order for K=3
#'
#' # Fixed-P filtering
#' Pfix <- matrix(c(0.9, 0.1,
#'                  0.2, 0.8), nrow = 2, byrow = TRUE)
#' out_fix <- rsdc_hamilton(y = y, X = NULL, beta = NULL,
#'                          rho_matrix = rho, K = K, N = N, P = Pfix)
#' str(out_fix$filtered_probs)
#'
#' # TVTP filtering (include an intercept yourself)
#' X <- cbind(1, scale(seq_len(T)))
#' beta <- rbind(c(1.2, 0.0),
#'               c(0.8, -0.1))
#' out_tvtp <- rsdc_hamilton(y = y, X = X, beta = beta,
#'                           rho_matrix = rho, K = K, N = N)
#' out_tvtp$log_likelihood
#'
#' @seealso \code{\link{rsdc_likelihood}} and \code{\link{rsdc_estimate}}.
#'
#' @references
#' \insertRef{hamilton1989}{RSDC}
#'
#' @note For \eqn{N=2}, TVTP uses a logistic link on the diagonal; off-diagonals follow by
#'   complement. For \eqn{N \ge 3}, a full softmax is used with \eqn{(N-1)} free logit
#'   vectors per row; all entries are independently determined.
#'
#'   X-timing: uses \code{X[t, ]} (contemporaneous) to form \eqn{P_t}, consistent with
#'   \code{\link{rsdc_simulate}}.
#'
#' @importFrom stats plogis
#' @export
rsdc_hamilton <- function(y, X = NULL, beta = NULL, rho_matrix, K, N, P = NULL,
                          xi_init = NULL, engine = c("cpp", "r")) {
  engine <- match.arg(engine)

  if (!is.matrix(y)) stop("y must be a numeric matrix.")
  if (any(!is.finite(y))) stop("y must not contain NA/NaN/Inf values.")
  if (!is.null(X) && !is.matrix(X)) stop("X must be a numeric matrix or NULL.")
  if (!is.null(X) && any(!is.finite(X))) stop("X must not contain NA/NaN/Inf values.")
  if (!is.null(X) && nrow(X) != nrow(y))
    stop(sprintf("X must have the same number of rows as y (nrow(X) = %d, nrow(y) = %d).",
                 nrow(X), nrow(y)))
  # TVTP needs both X and beta; supplying only one silently fell back to a fixed P.
  if (xor(is.null(X), is.null(beta)))
    stop("For time-varying transitions supply both X and beta; supply neither for a ",
         "fixed transition matrix.")
  if (!is.null(beta) && (!is.matrix(beta) || nrow(beta) != N)) stop("beta must be a matrix with N rows.")
  if (!is.null(beta) && !is.null(X)) {
    expected_ncol <- if (N == 2L) ncol(X) else (N - 1L) * ncol(X)
    if (ncol(beta) != expected_ncol)
      stop(sprintf("beta must have %d column(s) for N=%d and p=%d (got %d).",
                   expected_ncol, N, ncol(X), ncol(beta)))
  }
  if (!is.matrix(rho_matrix) || nrow(rho_matrix) != N || ncol(rho_matrix) != K*(K - 1)/2) {
    stop("rho_matrix must be of dimension N x (K*(K-1)/2).")
  }
  if (!is.null(P)) {
    if (!is.matrix(P) || !all(dim(P) == c(N, N)))
      stop("P must be an N x N matrix if provided.")
    if (any(!is.finite(P)) || any(P < 0) || any(abs(rowSums(P) - 1) > 1e-8))
      stop("P must be a valid row-stochastic matrix (finite, non-negative entries, ",
           "rows summing to 1).")
  }
  if (!is.null(xi_init) && (length(xi_init) != N || any(!is.finite(xi_init)) ||
                            any(xi_init < 0) || sum(xi_init) <= 0)) {
    stop("xi_init must be a finite, non-negative vector of length N with a positive sum.")
  }

  n_obs <- nrow(y)
  if (K != ncol(y)) stop(sprintf("K = %d but ncol(y) = %d.", K, ncol(y)))

  # Fast path: run the full filter + smoother in C++ (engine = "cpp", default).
  # The pure-R engine below is kept as the reference and for equivalence testing.
  if (engine == "cpp") {
    sig <- array(0, dim = c(K, K, N))
    for (m in 1:N) {
      Rm <- diag(K)
      Rm[lower.tri(Rm)] <- rho_matrix[m, ]
      Rm[upper.tri(Rm)] <- t(Rm)[upper.tri(Rm)]
      sig[, , m] <- Rm
    }
    tvtp  <- !is.null(X) && !is.null(beta)
    empty <- matrix(0, 0, 0)
    P0    <- if (tvtp) empty else if (!is.null(P)) P else matrix(1 / N, N, N)
    res <- rsdc_filter_full_cpp(
      as.matrix(y), sig, if (tvtp) 1L else 0L, P0,
      if (tvtp) as.matrix(X) else empty,
      if (tvtp) as.matrix(beta) else empty,
      if (is.null(xi_init)) numeric(0) else xi_init)
    if (!isTRUE(res$ok) || !is.finite(res$log_likelihood))
      return(list(filtered_probs = NULL, smoothed_probs = NULL,
                  log_likelihood = -Inf, loglik_t = NULL))
    return(list(filtered_probs = res$filtered_probs,
                smoothed_probs = res$smoothed_probs,
                log_likelihood = res$log_likelihood,
                loglik_t       = as.numeric(res$loglik_t)))
  }

  # Build regime correlation matrices
  sigma <- array(dim = c(K, K, N))
  for (m in 1:N) {
    R <- diag(K)
    R[lower.tri(R)] <- rho_matrix[m, ]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]

    # Check positive definiteness
    if (min(eigen(R)$values) < 1e-8)
      return(list(filtered_probs = NULL, smoothed_probs = NULL, log_likelihood = -Inf))

    sigma[,,m] <- R
  }

  # Transition probability matrices
  P_mats <- array(dim = c(N, N, n_obs))

  if (is.null(X) || is.null(beta)) {
    # No exogenous variables, use fixed transition matrix P
    if (is.null(P)) {
      # If no P matrix is provided, use equal probabilities
      for (t in 1:n_obs) {
        P_mats[,,t] <- matrix(1/N, N, N)
      }
    } else {
      # Use provided P matrix
      for (t in 1:n_obs) {
        P_mats[,,t] <- P
      }
    }
  } else {
    # Time-varying transition probabilities
    p_cov <- ncol(X)
    for (t in 1:n_obs) {
      for (i in 1:N) {
        if (N == 1) {
          P_mats[,,t] <- matrix(1)
          next
        } else if (N == 2) {
          # Logistic link: one free parameter vector per regime
          p_ii <- plogis(X[t,] %*% beta[i,])
          off_prob <- (1 - p_ii) / (N - 1)
          P_mats[i,,t] <- off_prob
          P_mats[i,i,t] <- p_ii
        } else {
          # Softmax: (N-1) free beta vectors per regime row; last column is reference (0)
          logits <- numeric(N)
          for (j in seq_len(N - 1)) {
            logits[j] <- sum(X[t, ] * beta[i, ((j - 1L) * p_cov + 1L):(j * p_cov)])
          }
          logits[N] <- 0
          exp_l <- exp(logits - max(logits))
          P_mats[i,,t] <- exp_l / sum(exp_l)
        }
      }
    }
  }

  # Multivariate normal log-densities
  log_densities <- matrix(nrow = N, ncol = n_obs)
  for (m in 1:N) {
    sigma_m   <- sigma[,,m]
    sigma_rdg <- sigma_m + diag(1e-8, K)   # ridge used consistently for inverse and log-det
    inv_sigma <- try(solve(sigma_rdg), silent = TRUE)
    if (inherits(inv_sigma, "try-error"))
      return(list(filtered_probs = NULL, smoothed_probs = NULL, log_likelihood = -Inf))

    log_det <- determinant(sigma_rdg, logarithm = TRUE)$modulus[1]
    centered <- as.matrix(y)
    inv_sigma <- as.matrix(inv_sigma)
    centered <- matrix(apply(centered, 2, as.numeric), nrow = n_obs)
    inv_sigma <- matrix(as.numeric(inv_sigma), nrow = nrow(inv_sigma))
    quad_form <- rowSums((centered %*% inv_sigma) * centered)
    log_densities[m,] <- -0.5*(log_det + quad_form + K*log(2*pi))
  }

  # Filtering and smoothing
  filtered <- matrix(0, N, n_obs)
  predicted <- matrix(0, N, n_obs)
  smoothed <- matrix(0, N, n_obs)
  loglik_t <- numeric(n_obs)                  # per-observation log-likelihood contributions
  xi <- if (is.null(xi_init)) rep(1/N, N) else xi_init / sum(xi_init)
  log_lik <- 0

  for (t in 1:n_obs) {
    P_t <- P_mats[,,t]
    predicted[,t] <- t(P_t) %*% xi

    # Log-space update: factor out the largest log-density to avoid exp() underflow.
    # exp(c_t) cancels in the normalization, so filtered probs are unchanged; the
    # log-likelihood increment log(sum_w) + c_t = log(sum predicted * exp(logdens)).
    ld    <- log_densities[, t]
    c_t   <- max(ld)
    w     <- predicted[, t] * exp(ld - c_t)
    sum_w <- sum(w)

    if (!is.finite(sum_w) || sum_w <= 0)
      return(list(filtered_probs = NULL, smoothed_probs = NULL,
                  log_likelihood = -Inf, loglik_t = NULL))

    filtered[,t] <- w / sum_w
    loglik_t[t]  <- log(sum_w) + c_t
    log_lik <- log_lik + loglik_t[t]
    xi <- filtered[,t]
  }

  # Smoothing
  smoothed[,n_obs] <- filtered[,n_obs]
  if (n_obs > 1L) {
    for (t in (n_obs - 1):1) {
      P_t1  <- P_mats[,, t + 1]
      ratio <- smoothed[, t + 1] / pmax(predicted[, t + 1], .Machine$double.eps)
      temp  <- filtered[, t] * (P_t1 %*% ratio)
      s <- sum(temp)
      if (!is.finite(s) || s < .Machine$double.eps)
        return(list(filtered_probs = NULL, smoothed_probs = NULL,
                    log_likelihood = -Inf, loglik_t = NULL))
      smoothed[, t] <- temp / s
    }
  }

  list(
    filtered_probs = filtered,
    smoothed_probs = smoothed,
    log_likelihood = log_lik,
    loglik_t       = loglik_t
  )
}

#' Negative Log-Likelihood for Regime-Switching Correlation Models
#'
#' Computes the negative log-likelihood for a multivariate \emph{correlation-only}
#' regime-switching model, with either a fixed (time-invariant) transition matrix or
#' time-varying transition probabilities (TVTP) driven by exogenous covariates.
#' Likelihood evaluation uses the Hamilton (1989) filter.
#'
#' @param params Numeric vector of model parameters packed as:
#' \itemize{
#'   \item \strong{No exogenous covariates} (\code{exog = NULL}): first
#'         \eqn{N(N-1)} transition parameters (for the fixed transition matrix),
#'         followed by \eqn{N \times K(K-1)/2} correlation parameters, stacked
#'         \emph{row-wise by regime} in \code{lower.tri} order.
#'   \item \strong{With exogenous covariates} (\code{exog} provided): first
#'         \eqn{N \times p} TVTP coefficients (\code{beta}) for \eqn{N=2}, or
#'         \eqn{N \times (N-1)p} coefficients for \eqn{N \ge 3} (row \eqn{i}
#'         holds \eqn{N-1} length-\eqn{p} softmax vectors; last logit = 0 reference),
#'         followed by \eqn{N \times K(K-1)/2} correlation parameters,
#'         stacked \emph{row-wise by regime} in \code{lower.tri} order.
#' }
#' @param y Numeric matrix \eqn{T \times K} of observations (e.g., standardized residuals).
#'          Columns are treated as mean-zero, unit-variance; only correlation is modeled.
#' @param exog Optional numeric matrix \eqn{T \times p} of exogenous covariates.
#'        If supplied, a TVTP specification is used.
#' @param K Integer. Number of observed series (columns of \code{y}).
#' @param N Integer. Number of regimes.
#'
#' @returns Numeric scalar: the \emph{negative} log-likelihood to be minimized by an optimizer.
#'
#' @details
#' \itemize{
#'   \item \strong{Transition dynamics:}
#'         \itemize{
#'           \item \emph{Fixed P (no \code{exog}):} \code{params} begins with transition
#'                 parameters. For \eqn{N=2}, the implementation maps them to
#'                 \eqn{P=\begin{pmatrix} p_{11} & 1-p_{11}\\ 1-p_{22} & p_{22}\end{pmatrix}}.
#'                 For \eqn{N \ge 3}, row \eqn{i} of \eqn{P} holds \eqn{N-1} free
#'                 probabilities (packed row-wise in \code{params}) filling columns
#'                 \eqn{1,\dots,N-1}; the last column is the row complement
#'                 \eqn{1-\sum_{j<N} P_{ij}}. Proposals with a negative entry
#'                 (row sum of free probabilities above 1) are penalized.
#'           \item \emph{TVTP:} with \code{exog}, for \eqn{N=2},
#'                 \eqn{p_{ii,t} = \mathrm{logit}^{-1}(X_t^\top \beta_i)} and
#'                 \eqn{p_{ij,t} = 1 - p_{ii,t}} (\eqn{j \ne i});
#'                 for \eqn{N \ge 3}, a softmax over \eqn{N-1} free logit vectors per row
#'                 with the \eqn{N}-th logit fixed at 0 (reference category).
#'         }
#'   \item \strong{Correlation build:} per regime, the lower-triangular vector is filled into
#'         a symmetric correlation matrix. Non-positive-definite proposals or \eqn{|\rho|\ge 1}
#'         are penalized via a large objective value.
#'   \item \strong{Evaluation:} delegates to \code{\link{rsdc_hamilton}}; if the filter returns
#'         \code{log_likelihood = -Inf}, a large penalty is returned.
#' }
#'
#' @examples
#' # Small toy example (N = 2, K = 3), fixed P (no exog)
#' set.seed(1)
#' T <- 40; K <- 3; N <- 2
#' y <- scale(matrix(rnorm(T * K), T, K), center = TRUE, scale = TRUE)
#'
#' # Pack parameters: trans (p11, p22), then rho by regime (lower-tri order)
#' p11 <- 0.9; p22 <- 0.8
#' rho1 <- c(0.10, 0.05, 0.00)  # (2,1), (3,1), (3,2)
#' rho2 <- c(0.60, 0.40, 0.30)
#' params <- c(p11, p22, rho1, rho2)
#'
#' nll <- rsdc_likelihood(params, y = y, exog = NULL, K = K, N = N)
#' nll
#'
#' # TVTP example: add X and beta (pack beta row-wise, then rho)
#' X <- cbind(1, scale(seq_len(T)))
#' beta <- rbind(c(1.2, 0.0),
#'               c(0.8, -0.1))
#' params_tvtp <- c(as.vector(t(beta)), rho1, rho2)
#' nll_tvtp <- rsdc_likelihood(params_tvtp, y = y, exog = X, K = K, N = N)
#' nll_tvtp
#'
#' @seealso \code{\link{rsdc_hamilton}} (filter),
#' \code{\link[stats]{optim}}, and \code{\link[DEoptim]{DEoptim}}
#'
#' @note The function is written for use inside optimizers; it performs inexpensive validation
#'       and returns large penalties for invalid parameterizations instead of stopping with errors.
#'
#' @export
rsdc_likelihood <- function(params, y, exog = NULL, K, N) {
  if (K < 2L) return(1e10)  # correlation model is undefined for a single series
  if (any(is.na(params)) || any(!is.finite(params))) return(1e10)

  # Parameter count calculation
  n_p <- ifelse(is.null(exog), N * (N - 1), 0)
  n_pairs <- K * (K - 1) / 2
  n_rho <- N * n_pairs
  n_beta <- if (!is.null(exog)) {
    if (N == 2) N * ncol(exog) else N * (N - 1L) * ncol(exog)
  } else 0

  # Reject a wrong-length parameter vector (avoids NA from out-of-range slicing).
  if (length(params) != (if (is.null(exog)) n_p else n_beta) + n_rho) return(1e10)

  # Parameter extraction
  if (is.null(exog)) {
    # Without exogenous variables, extract transition parameters first
    trans_params <- params[1:n_p]
    rho_matrix <- matrix(params[(n_p + 1):(n_p + n_rho)], nrow = N, byrow = TRUE)
    beta <- NULL

    # Build transition matrix. N=2 keeps the diagonal-stay parameterization;
    # for N >= 3 the free entries fill the first N-1 columns and the last column
    # is the row complement (generic for any N).
    P <- matrix(0, N, N)
    if (N == 2) {
      if (any(trans_params < 0 | trans_params > 1)) return(1e10)
      P <- matrix(c(trans_params[1], 1 - trans_params[1],
                    1 - trans_params[2], trans_params[2]),
                  nrow = N, byrow = TRUE)
    } else {
      for (i in 1:N) {
        free <- trans_params[((i - 1L) * (N - 1L) + 1L):(i * (N - 1L))]
        last <- 1 - sum(free)
        if (any(free < 0) || last < 0) return(1e10)
        P[i, ] <- c(free, last)
      }
    }
  } else {
    # With exogenous variables (N=2: logistic beta; N=3: softmax beta)
    beta <- matrix(params[1:n_beta], nrow = N, byrow = TRUE)
    rho_matrix <- matrix(params[(n_beta + 1):(n_beta + n_rho)], nrow = N, byrow = TRUE)
    P <- NULL
  }

  # Check correlation validity
  if (any(abs(rho_matrix) >= 1)) return(1e10)

  # Build regime correlation matrices (identical to rsdc_hamilton)
  sigma <- array(0, dim = c(K, K, N))
  for (m in 1:N) {
    Rm <- diag(K)
    Rm[lower.tri(Rm)] <- rho_matrix[m, ]
    Rm[upper.tri(Rm)] <- t(Rm)[upper.tri(Rm)]
    sigma[, , m] <- Rm
  }

  # Fast C++ log-likelihood (mirrors the rsdc_hamilton filter exactly).
  empty <- matrix(0, 0, 0)
  ll <- if (is.null(exog))
    rsdc_loglik_cpp(as.matrix(y), sigma, 0L, P, empty, empty, numeric(0))
  else
    rsdc_loglik_cpp(as.matrix(y), sigma, 1L, empty, as.matrix(exog),
                    as.matrix(beta), numeric(0))

  if (!is.finite(ll)) return(1e10)
  return(-ll)
}
