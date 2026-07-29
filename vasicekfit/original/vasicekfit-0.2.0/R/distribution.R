#' The Vasicek Loss Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Vasicek credit loss distribution.
#'
#' @param x,q vector of quantiles (loss rates in (0, 1)).
#' @param prob vector of probabilities.
#' @param n number of observations to generate.
#' @param p probability of default, in (0, 1).
#' @param rho asset value correlation, in (0, 1).
#' @param kappa optional numeric vector of macro-factor sensitivities.
#' @param u optional numeric vector of macro-factor values, same length as
#'   \code{kappa}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities are
#'   given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @details
#' The Vasicek loss distribution arises from the single-factor Gaussian copula
#' model. When \code{kappa} and \code{u} are supplied, the effective PD is
#' \eqn{\Phi(\Phi^{-1}(p) + \sum \kappa_j u_j)}.
#'
#' @return \code{dvasicek} gives the density, \code{pvasicek} gives the
#'   distribution function, \code{qvasicek} gives the quantile function, and
#'   \code{rvasicek} generates random deviates.
#'
#' @examples
#' # Standard Vasicek density
#' curve(dvasicek(x, p = 0.03, rho = 0.10), from = 0.001, to = 0.20,
#'       ylab = "Density", main = "Vasicek loss distribution")
#'
#' # 99th percentile
#' qvasicek(0.99, p = 0.03, rho = 0.10)
#'
#' # Random sample
#' hist(rvasicek(10000, p = 0.03, rho = 0.10), breaks = 100)
#'
#' @name Vasicek
NULL

.vasicek_mu <- function(p, kappa = NULL, u = NULL) {
  mu <- stats::qnorm(p)
  if (!is.null(kappa) && !is.null(u)) {
    if (length(kappa) != length(u)) {
      stop("kappa and u must have the same length")
    }
    mu <- mu + sum(kappa * u)
  }
  mu
}

#' @rdname Vasicek
#' @export
dvasicek <- function(x, p, rho, kappa = NULL, u = NULL, log = FALSE) {
  if (p <= 0 || p >= 1) stop("p must be in (0, 1)")
  if (rho <= 0 || rho >= 1) stop("rho must be in (0, 1)")

  mu <- .vasicek_mu(p, kappa, u)
  z <- suppressWarnings(stats::qnorm(x)) # NaN outside (0, 1); overwritten below
  arg <- (sqrt(1 - rho) * z - mu) / sqrt(rho)

  log_d <- 0.5 * log((1 - rho) / rho) +
    0.5 * (z^2 - arg^2)

  # Handle boundary values
  log_d[x <= 0 | x >= 1] <- -Inf

  if (log) log_d else exp(log_d)
}

#' @rdname Vasicek
#' @export
pvasicek <- function(q, p, rho, kappa = NULL, u = NULL,
                     lower.tail = TRUE, log.p = FALSE) {
  if (p <= 0 || p >= 1) stop("p must be in (0, 1)")
  if (rho <= 0 || rho >= 1) stop("rho must be in (0, 1)")

  mu <- .vasicek_mu(p, kappa, u)

  z <- suppressWarnings(stats::qnorm(q)) # NaN outside (0, 1); overwritten below
  prob <- stats::pnorm(
    (sqrt(1 - rho) * z - mu) / sqrt(rho)
  )

  prob[q <= 0] <- 0
  prob[q >= 1] <- 1

  if (!lower.tail) prob <- 1 - prob
  if (log.p) prob <- log(prob)
  prob
}

#' @rdname Vasicek
#' @export
qvasicek <- function(prob, p, rho, kappa = NULL, u = NULL,
                     lower.tail = TRUE, log.p = FALSE) {
  if (p <= 0 || p >= 1) stop("p must be in (0, 1)")
  if (rho <= 0 || rho >= 1) stop("rho must be in (0, 1)")

  if (log.p) prob <- exp(prob)
  if (!lower.tail) prob <- 1 - prob

  mu <- .vasicek_mu(p, kappa, u)

  stats::pnorm(
    (mu + sqrt(rho) * stats::qnorm(prob)) / sqrt(1 - rho)
  )
}

#' @rdname Vasicek
#' @export
rvasicek <- function(n, p, rho, kappa = NULL, u = NULL) {
  if (p <= 0 || p >= 1) stop("p must be in (0, 1)")
  if (rho <= 0 || rho >= 1) stop("rho must be in (0, 1)")

  mu <- .vasicek_mu(p, kappa, u)
  z <- stats::rnorm(n)

  stats::pnorm((mu + sqrt(rho) * z) / sqrt(1 - rho))
}
