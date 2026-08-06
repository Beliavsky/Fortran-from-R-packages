#' Simulate Benchmark and Test-Asset Returns for Spanning Studies
#'
#' Generates a panel of benchmark returns \eqn{R_1} (\eqn{n \times K}) and
#' test-asset returns \eqn{R_2} (\eqn{n \times N}) from a factor model with a
#' controllable mean-variance spanning violation, matching the data-generating
#' processes used in the size/power study of Ardia and Sessinou (2025). It is a
#' fast, self-contained replacement for ad-hoc \code{fGarch::garchSim()} loops:
#' the innovations are drawn in one vectorised call and the AR / GARCH dynamics
#' are applied with a lightweight recursion.
#'
#' @details
#' The latent factors \eqn{z} are \eqn{K} innovation series, cross-sectionally
#' correlated through the Cholesky factor of a Toeplitz matrix with decay
#' \code{rho_factor}; the idiosyncratic terms are \eqn{N} innovation series with
#' Toeplitz decay \code{rho_error}. Test assets are
#' \deqn{R_2 = \alpha + z B + \varepsilon, \qquad
#'       B = (1+\text{ncp})\,\mathbf{1}_{K\times N},\quad
#'       \alpha_j = \text{ncp},}
#' and the benchmarks are \eqn{R_1 = [\,z_1,\; z_{-1}+z_1\,]}. Under
#' \code{ncp = 0} the benchmarks span the test assets (\eqn{\alpha = 0}); larger
#' \code{|ncp|} moves the intercepts and loadings away from the spanning null.
#' With \code{sparse = TRUE} the first \eqn{\lfloor N/2 \rfloor} intercepts are
#' set to zero (a sparse alternative); \eqn{\alpha} is added column-wise, i.e.
#' asset \eqn{j} receives \eqn{\alpha_j}.
#'
#' Each innovation series is defined by \code{innovation} (marginal law) and
#' \code{dynamics} (serial structure): \code{"iid"}, an AR(1) with coefficient
#' \code{ar}, a GARCH(1,1) with parameters \code{garch}, or their
#' \code{"ar-garch"} combination. Student-\eqn{t} innovations are standardised to
#' unit variance when \code{standardize = TRUE}; skew-\eqn{t} innovations
#' (standardised Fernandez-Steel skew-\eqn{t}) are always unit-variance, so
#' \code{standardize} has no effect on them.
#'
#' \code{dgp} is a convenience argument selecting one of the twelve processes of
#' Ardia and Sessinou (2025); when supplied it overrides \code{innovation},
#' \code{dynamics}, \code{df}, \code{xi} and \code{standardize}:
#' \tabular{rlrl}{
#'   1 \tab normal, iid                \tab 7  \tab normal, AR-GARCH \cr
#'   2 \tab Student-t(5), iid (raw)     \tab 8  \tab Student-t(4), AR-GARCH \cr
#'   3 \tab skew-t(4, 0.9), iid         \tab 9  \tab skew-t(4, 0.9), AR-GARCH \cr
#'   4 \tab normal, GARCH               \tab 10 \tab normal, AR \cr
#'   5 \tab Student-t(4), GARCH         \tab 11 \tab Student-t(5), AR (raw) \cr
#'   6 \tab skew-t(4, 0.9), GARCH       \tab 12 \tab skew-t(4, 0.9), AR
#' }
#' DGPs 2 and 11 use \emph{raw} (non-standardised) \eqn{t_5} innovations, as in
#' the original study; all other heavy-tailed processes are standardised.
#'
#' @param n Integer, number of time periods (sample size).
#' @param K Integer, number of benchmark assets.
#' @param N Integer, number of test assets.
#' @param ncp Non-centrality controlling the spanning violation (0 = null).
#' @param rho_factor,rho_error Toeplitz correlation decays for the factors and
#'   the idiosyncratic terms, in \eqn{[0, 1)}.
#' @param innovation Marginal law of the innovations: \code{"normal"},
#'   \code{"t"}, or \code{"skew-t"}. Ignored if \code{dgp} is supplied.
#' @param dynamics Serial structure: \code{"iid"}, \code{"ar"}, \code{"garch"},
#'   or \code{"ar-garch"}. Ignored if \code{dgp} is supplied.
#' @param sparse Logical; if \code{TRUE}, zero the first \eqn{\lfloor N/2\rfloor}
#'   intercepts (sparse alternative).
#' @param df Degrees of freedom for the Student-/skew-\eqn{t} innovations.
#' @param xi Skewness parameter for the skew-\eqn{t} innovations.
#' @param ar AR(1) coefficient (used by \code{"ar"} and \code{"ar-garch"}).
#' @param garch Named numeric vector \code{c(omega, alpha, beta)} for the
#'   GARCH(1,1) variance recursion.
#' @param standardize Logical; standardise \eqn{t} innovations to unit variance.
#' @param dgp Optional integer in \code{1:12}; selects a preset process (see
#'   Details) and overrides \code{innovation}/\code{dynamics}/\code{df}/\code{xi}/
#'   \code{standardize}.
#' @param burnin Integer, number of initial observations discarded to remove the
#'   AR/GARCH transient.
#'
#' @return A list with components \code{R1} (\eqn{n \times K} benchmark returns)
#'   and \code{R2} (\eqn{n \times N} test-asset returns), suitable as the
#'   \code{R1}/\code{R2} arguments of the spanning tests in this package.
#'
#' @note Set the RNG state with \code{\link{set.seed}} beforehand for
#'   reproducibility.
#'
#' @references
#' \insertRef{ArdiaSessinou2025}{spantest} \cr
#' \insertRef{GungorLuger2016}{spantest} \cr
#'
#' @examples
#' set.seed(1)
#' sim <- span_simulate(n = 250, K = 3, N = 10, ncp = 0)       # spanning null
#' span_grs(sim$R1, sim$R2)$pval
#'
#' sim2 <- span_simulate(n = 250, K = 3, N = 10, ncp = 0.2)    # alternative
#' span_grs(sim2$R1, sim2$R2)$pval
#'
#' @seealso [span_grs()], [span_gl_ad()], [span_as()]
#'
#' @importFrom stats rnorm rt runif filter toeplitz
#' @export
span_simulate <- function(n, K, N, ncp = 0,
                          rho_factor = 0.8, rho_error = 0.5,
                          innovation = c("normal", "t", "skew-t"),
                          dynamics = c("iid", "ar", "garch", "ar-garch"),
                          sparse = FALSE, df = 5, xi = 0.9, ar = 0.2,
                          garch = c(omega = 0.1, alpha = 0.1, beta = 0.8),
                          standardize = TRUE, dgp = NULL, burnin = 500L) {

  innovation <- match.arg(innovation)
  dynamics   <- match.arg(dynamics)

  if (!is.null(dgp)) {
    dgp <- as.integer(dgp)
    if (length(dgp) != 1L || is.na(dgp) || dgp < 1L || dgp > 12L) {
      stop("`dgp` must be a single integer in 1:12.")
    }
    presets <- list(
      c("normal",  "iid",      "5", "TRUE"),  c("t",      "iid",      "5", "FALSE"),
      c("skew-t",  "iid",      "4", "TRUE"),  c("normal", "garch",    "5", "TRUE"),
      c("t",       "garch",    "4", "TRUE"),  c("skew-t", "garch",    "4", "TRUE"),
      c("normal",  "ar-garch", "5", "TRUE"),  c("t",      "ar-garch", "4", "TRUE"),
      c("skew-t",  "ar-garch", "4", "TRUE"),  c("normal", "ar",       "5", "TRUE"),
      c("t",       "ar",       "5", "FALSE"), c("skew-t", "ar",       "4", "TRUE"))
    p <- presets[[dgp]]
    innovation  <- p[1]
    dynamics    <- p[2]
    df          <- as.numeric(p[3])
    standardize <- as.logical(p[4])
    xi          <- 0.9
  }

  stopifnot(
    "n, K and N must be positive" = n >= 1 && K >= 1 && N >= 1,
    "rho_factor and rho_error must be in [0, 1)" =
      rho_factor >= 0 && rho_factor < 1 && rho_error >= 0 && rho_error < 1,
    "df must be greater than 2" = df > 2,
    "xi must be positive" = xi > 0
  )
  if (dynamics %in% c("garch", "ar-garch")) {
    stopifnot(
      "garch must be a named vector c(omega, alpha, beta), omega > 0, alpha, beta >= 0" =
        all(c("omega", "alpha", "beta") %in% names(garch)) &&
        garch[["omega"]] > 0 && garch[["alpha"]] >= 0 && garch[["beta"]] >= 0,
      "GARCH must be stationary: alpha + beta < 1" =
        garch[["alpha"]] + garch[["beta"]] < 1
    )
  }
  if (dynamics %in% c("ar", "ar-garch")) {
    stopifnot("AR(1) must be stationary: abs(ar) < 1" = abs(ar) < 1)
  }

  om <- garch[["omega"]]; al <- garch[["alpha"]]; be <- garch[["beta"]]

  one_series <- function() {
    m <- n + burnin
    z <- switch(innovation,
      "normal" = stats::rnorm(m),
      "t"      = { r <- stats::rt(m, df); if (standardize) r / sqrt(df / (df - 2)) else r },
      "skew-t" = f_rsstd(m, df, xi))
    s <- switch(dynamics,
      "iid"      = z,
      "garch"    = garch_filter(z, om, al, be),
      "ar"       = as.numeric(stats::filter(z, ar, method = "recursive")),
      "ar-garch" = as.numeric(stats::filter(garch_filter(z, om, al, be), ar, method = "recursive")))
    s[(burnin + 1L):m]
  }

  Z <- matrix(0, n, K); for (j in seq_len(K)) Z[, j] <- one_series()
  E <- matrix(0, n, N); for (j in seq_len(N)) E[, j] <- one_series()

  z   <- Z %*% chol(stats::toeplitz(rho_factor ^ (seq_len(K) - 1L)))
  eps <- E %*% chol(stats::toeplitz(rho_error  ^ (seq_len(N) - 1L)))

  alpha <- rep(ncp, N)
  if (sparse) alpha[seq_len(floor(N / 2))] <- 0
  Beta <- matrix(1 + ncp, nrow = K, ncol = N)

  y <- sweep(z %*% Beta + eps, 2, alpha, "+")          # column-wise intercepts
  x <- cbind(z[, 1], z[, -1, drop = FALSE] + z[, 1])   # benchmark construction

  list(R1 = x, R2 = y)
}

#' Standardised skew-Student-t random draws (Fernandez-Steel)
#'
#' Zero-mean, unit-variance skewed Student-t innovations, a base-R re-
#' implementation of \code{fGarch::rsstd(n, nu, xi)} (identical draws given the
#' same RNG state), used so the package does not depend on \pkg{fGarch}.
#'
#' @param n Integer, number of draws.
#' @param nu Degrees of freedom (\code{> 2}).
#' @param xi Skewness parameter (\code{> 0}; 1 = symmetric).
#' @return Numeric vector of length \code{n}.
#' @keywords internal
#' @noRd
f_rsstd <- function(n, nu, xi) {
  weight <- xi / (xi + 1 / xi)
  z  <- stats::runif(n, -weight, 1 - weight)
  Xi <- xi ^ sign(z)
  rt_std <- stats::rt(n, df = nu) / sqrt(nu / (nu - 2))   # standardised Student-t
  rnd <- -abs(rt_std) / Xi * sign(z)
  m1 <- 2 * sqrt(nu - 2) / (nu - 1) / beta(1 / 2, nu / 2)
  mu <- m1 * (xi - 1 / xi)
  sigma <- sqrt((1 - m1^2) * (xi^2 + 1 / xi^2) + 2 * m1^2 - 1)
  (rnd - mu) / sigma
}
