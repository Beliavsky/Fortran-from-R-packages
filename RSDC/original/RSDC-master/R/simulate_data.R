#' Simulate Multivariate Regime-Switching Data (TVTP)
#'
#' Simulates a multivariate time series from a regime-switching model with
#' \emph{time-varying transition probabilities} (TVTP) driven by covariates \code{X}.
#' Transition probabilities are generated via a multinomial logistic (softmax) link;
#' observations are drawn from regime-specific Gaussian distributions.
#'
#' @param n Integer. Number of time steps to simulate; must be at least 2.
#' @param X Numeric matrix \eqn{n \times p} of covariates used to form the transition
#'   probabilities. Row \code{X[1, ]} forms \eqn{P_1}, used only to draw the initial
#'   state from \eqn{t(P_1)\,\pi_0}; rows \code{X[t, ]} for \eqn{t \ge 2} form \eqn{P_t}.
#' @param beta Numeric matrix \eqn{N \times (N-1) \cdot p}. Transition coefficients
#'   packed row-wise by regime, matching the parameterization of
#'   \code{\link{rsdc_hamilton}} and \code{\link{rsdc_estimate}}. For \eqn{N=2},
#'   each row is a length-\eqn{p} logistic vector:
#'   \eqn{p_{ii,t} = \mathrm{logit}^{-1}(X_t^\top \beta_i)}. For \eqn{N \ge 3},
#'   each row packs \eqn{N-1} softmax vectors of length \eqn{p}; the \eqn{N}-th
#'   logit is the reference (zero).
#' @param mu Numeric matrix \eqn{N \times K}. Regime-specific mean vectors.
#' @param sigma Numeric array \eqn{K \times K \times N}. Regime-specific covariance
#'   (here, correlation/variance) matrices; each \eqn{K \times K} slice must be
#'   symmetric positive definite.
#' @param N Integer. Number of regimes.
#' @param seed Optional integer. If supplied, sets the RNG seed for reproducibility.
#'
#' @returns A list with:
#' \describe{
#'   \item{states}{Integer vector of length \code{n}; the simulated regime index at each time.}
#'   \item{observations}{Numeric matrix \eqn{n \times K}; the simulated observations.}
#'   \item{transition_matrices}{Array \eqn{N \times N \times n}; the transition matrix \eqn{P_t}
#'         used at each time step (with \eqn{P_1} undefined by construction; see Details).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Initial state and first draw:} The initial regime \eqn{S_1} is drawn from
#'         \eqn{t(P_1)\,\pi_0} with \eqn{\pi_0 = (1/N,\dots,1/N)} and \eqn{P_1} built from
#'         \eqn{X_1} via the same logistic/softmax link, matching the Hamilton-filter t=1
#'         prior in \code{\link{rsdc_hamilton}}; the first observation \eqn{y_1} is drawn
#'         from \eqn{\mathcal{N}(\mu_{S_1}, \Sigma_{S_1})}.
#'   \item \strong{TVTP transition:} For \eqn{t \ge 2} and \eqn{N=2},
#'         \eqn{p_{ii,t} = \mathrm{logit}^{-1}(X_t^\top \beta_i)}; for
#'         \eqn{N \ge 3}, a softmax is used with \eqn{N-1} free logit vectors per
#'         row and the \eqn{N}-th logit fixed at zero (reference category).
#'         Both use log-sum-exp stabilization.
#'   \item \strong{Sampling:} Given \eqn{S_{t-1}}, draw \eqn{S_t} from the categorical
#'         distribution with probabilities \eqn{P_t(S_{t-1}, \cdot)} and
#'         \eqn{y_t \sim \mathcal{N}(\mu_{S_t}, \Sigma_{S_t})}.
#'   \item \strong{First-slice NA:} \code{transition_matrices[,,1]} is always \code{NA}.
#'         At \eqn{t=1} there is no predecessor state, so no transition matrix is computed.
#'         The assignment loop begins at \eqn{t=2}. Users iterating over the returned array
#'         should skip the first slice or use \code{transition_matrices[,,2:n]}.
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 200; K <- 3; N <- 2; p <- 2
#' X <- cbind(1, scale(seq_len(n)))
#'
#' beta <- matrix(0, nrow = N, ncol = (N - 1) * p)
#' beta[1, ] <- c(1.2,  0.0)
#' beta[2, ] <- c(1.0, -0.1)
#'
#' mu <- rbind(c(0, 0, 0),
#'             c(0, 0, 0))
#' rho <- rbind(c(0.10, 0.05, 0.00),
#'              c(0.60, 0.40, 0.30))
#' Sig <- array(0, dim = c(K, K, N))
#' for (m in 1:N) {
#'   R <- diag(K); R[lower.tri(R)] <- rho[m, ]; R[upper.tri(R)] <- t(R)[upper.tri(R)]
#'   Sig[, , m] <- R
#' }
#' sim <- rsdc_simulate(n = n, X = X, beta = beta, mu = mu, sigma = Sig, N = N, seed = 99)
#'
#' @seealso \code{\link{rsdc_hamilton}} (filter/evaluation),
#'   \code{\link{rsdc_estimate}} (estimators),
#'   \code{\link{rsdc_forecast}} (forecasting)
#'
#' @note Requires \pkg{mvtnorm} for multivariate normal sampling (called as \code{mvtnorm::rmvnorm}).
#'
#'   X-timing: uses \code{X[t, ]} (contemporaneous) to form \eqn{P_t}, consistent with
#'   \code{\link{rsdc_hamilton}}.
#'
#' @export
rsdc_simulate <- function(n, X, beta, mu, sigma, N, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  K <- ncol(mu)
  p <- ncol(X)

  if (!is.matrix(beta) || nrow(beta) != N || ncol(beta) != (N - 1L) * p)
    stop(sprintf("beta must be a matrix of dim N x (N-1)*p = %d x %d.", N, (N - 1L) * p))
  if (!all(dim(mu) == c(N, K))) stop("mu must be a matrix of dim N x K")
  if (!all(dim(sigma) == c(K, K, N))) stop("sigma must be an array of dim K x K x N")
  for (m in seq_len(N)) {
    eig <- eigen(sigma[,,m], only.values = TRUE)$values
    if (any(eig < 1e-8))
      stop(sprintf("sigma[,,%d] is not positive definite (min eigenvalue = %g).", m, min(eig)))
  }
  if (n < 2L) stop("n must be at least 2.")
  if (nrow(X) != n) stop(sprintf("X must have n = %d rows (one per time step), got %d.", n, nrow(X)))

  states <- numeric(n)
  observations <- matrix(NA, nrow = n, ncol = K)
  transition_matrices <- array(NA, dim = c(N, N, n))

  # Initial state: align with the Hamilton filter, whose t=1 prior is
  # t(P_1) %*% pi_0 with pi_0 = (1/N, ..., 1/N) and P_1 built from X[1, ] using the
  # same logistic/softmax link. Draw S_1 from that same distribution so the simulated
  # S_1 marginal matches the filter (and X[1, ] enters the DGP exactly as it enters
  # the likelihood). transition_matrices[,,1] is left NA (P_1 is used only here).
  P_1 <- matrix(NA, nrow = N, ncol = N)
  for (i in 1:N) {
    if (N == 2) {
      p_ii       <- plogis(sum(X[1, ] * beta[i, ]))
      P_1[i,  i] <- p_ii
      P_1[i, -i] <- 1 - p_ii
    } else {
      logits <- numeric(N)
      for (j in seq_len(N - 1L))
        logits[j] <- sum(X[1, ] * beta[i, ((j - 1L) * p + 1L):(j * p)])
      logits[N] <- 0
      exp_l      <- exp(logits - max(logits))
      P_1[i, ]   <- exp_l / sum(exp_l)
    }
  }
  init_probs <- as.numeric(t(P_1) %*% rep(1 / N, N))
  states[1] <- sample(1:N, size = 1, prob = init_probs)
  observations[1, ] <- mvtnorm::rmvnorm(1, mean = mu[states[1], ], sigma = sigma[,,states[1]])

  # Iterate over time
  for (t in 2:n) {
    P_t <- matrix(NA, nrow = N, ncol = N)

    for (i in 1:N) {
      if (N == 2) {
        p_ii      <- plogis(sum(X[t, ] * beta[i, ]))
        P_t[i,  i] <- p_ii
        P_t[i, -i] <- 1 - p_ii
      } else {
        logits <- numeric(N)
        for (j in seq_len(N - 1L))
          logits[j] <- sum(X[t, ] * beta[i, ((j - 1L) * p + 1L):(j * p)])
        logits[N] <- 0
        exp_l    <- exp(logits - max(logits))
        P_t[i, ] <- exp_l / sum(exp_l)
      }
    }

    transition_matrices[,,t] <- P_t
    prev_state <- states[t - 1]
    states[t] <- sample(1:N, size = 1, prob = P_t[prev_state, ])
    observations[t, ] <- mvtnorm::rmvnorm(1, mean = mu[states[t], ], sigma = sigma[,,states[t]])
  }

  list(
    states = states,
    observations = observations,
    transition_matrices = transition_matrices
  )
}
