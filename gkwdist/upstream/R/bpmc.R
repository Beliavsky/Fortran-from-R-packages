# ============================================================================#
# MCDONALD (Mc) / BETA POWER DISTRIBUTION
# ============================================================================#
#
# Wrapper functions for the McDonald (Beta Power) distribution.
# C++ implementations are in src/mc.cpp
#
# Functions:
#   - dmc: Probability density function (PDF)
#   - pmc: Cumulative distribution function (CDF)
#   - qmc: Quantile function (inverse CDF)
#   - rmc: Random number generation
#   - llmc: Negative log-likelihood
#   - grmc: Gradient of negative log-likelihood
#   - hsmc: Hessian of negative log-likelihood
# ============================================================================#


# ----------------------------------------------------------------------------#
# 1. DENSITY FUNCTION (dmc)
# ----------------------------------------------------------------------------#

#' @title Density of the McDonald (Mc)/Beta Power Distribution Distribution
#' @author Lopes, J. E.
#' @keywords distribution density mcdonald
#'
#' @description
#' Computes the probability density function (PDF) for the McDonald (Mc)
#' distribution (also previously referred to as Beta Power) with parameters
#' \code{gamma} (\eqn{\gamma}), \code{delta} (\eqn{\delta}), and \code{lambda}
#' (\eqn{\lambda}). This distribution is defined on the interval (0, 1).
#'
#' @param x Vector of quantiles (values between 0 and 1).
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param lambda Shape parameter \code{lambda} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param log Logical; if \code{TRUE}, the logarithm of the density is
#'   returned (\eqn{\log(f(x))}). Default: \code{FALSE}.
#'
#' @return A vector of density values (\eqn{f(x)}) or log-density values
#'   (\eqn{\log(f(x))}). The length of the result is determined by the recycling
#'   rule applied to the arguments (\code{x}, \code{gamma}, \code{delta},
#'   \code{lambda}). Returns \code{0} (or \code{-Inf} if
#'   \code{log = TRUE}) for \code{x} outside the interval (0, 1), or
#'   \code{NaN} if parameters are invalid (e.g., \code{gamma <= 0},
#'   \code{delta < 0}, \code{lambda <= 0}).
#'
#' @details
#' The probability density function (PDF) of the McDonald (Mc) distribution
#' is given by:
#' \deqn{
#' f(x; \gamma, \delta, \lambda) = \frac{\lambda}{B(\gamma,\delta+1)} x^{\gamma \lambda - 1} (1 - x^\lambda)^\delta
#' }
#' for \eqn{0 < x < 1}, where \eqn{B(a,b)} is the Beta function
#' (\code{\link[base]{beta}}).
#'
#' The Mc distribution is a special case of the five-parameter
#' Generalized Kumaraswamy (GKw) distribution (\code{\link{dgkw}}) obtained
#' by setting the parameters \eqn{\alpha = 1} and \eqn{\beta = 1}.
#' It was introduced by McDonald (1984) and is related to the Generalized Beta
#' distribution of the first kind (GB1). When \eqn{\lambda=1}, it simplifies
#' to the standard Beta distribution with parameters \eqn{\gamma} and
#' \eqn{\delta+1}.
#'
#' @references
#' McDonald, J. B. (1984). Some generalized functions for the size distribution
#' of income. *Econometrica*, 52(3), 647-663.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' @seealso
#' \code{\link{dgkw}} (parent distribution density),
#' \code{\link{pmc}}, \code{\link{qmc}}, \code{\link{rmc}} (other Mc functions),
#' \code{\link[stats]{dbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' x_vals <- c(0.2, 0.5, 0.8)
#' gamma_par <- 2.0
#' delta_par <- 1.5
#' lambda_par <- 1.0 # Equivalent to Beta(gamma, delta+1)
#'
#' # Calculate density using dmc
#' densities <- dmc(x_vals, gamma_par, delta_par, lambda_par)
#' print(densities)
#' # Compare with Beta density
#' print(stats::dbeta(x_vals, shape1 = gamma_par, shape2 = delta_par + 1))
#'
#' # Calculate log-density
#' log_densities <- dmc(x_vals, gamma_par, delta_par, lambda_par, log = TRUE)
#' print(log_densities)
#'
#' # Compare with dgkw setting alpha = 1, beta = 1
#' densities_gkw <- dgkw(x_vals,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = lambda_par
#' )
#' print(paste("Max difference:", max(abs(densities - densities_gkw)))) # Should be near zero
#'
#' # Plot the density for different lambda values
#' curve_x <- seq(0.01, 0.99, length.out = 200)
#' curve_y1 <- dmc(curve_x, gamma = 2, delta = 3, lambda = 0.5)
#' curve_y2 <- dmc(curve_x, gamma = 2, delta = 3, lambda = 1.0) # Beta(2, 4)
#' curve_y3 <- dmc(curve_x, gamma = 2, delta = 3, lambda = 2.0)
#'
#' plot(curve_x, curve_y2,
#'   type = "l", main = "McDonald (Mc) Density (gamma=2, delta=3)",
#'   xlab = "x", ylab = "f(x)", col = "red", ylim = range(0, curve_y1, curve_y2, curve_y3)
#' )
#' lines(curve_x, curve_y1, col = "blue")
#' lines(curve_x, curve_y3, col = "green")
#' legend("topright",
#'   legend = c("lambda=0.5", "lambda=1.0 (Beta)", "lambda=2.0"),
#'   col = c("blue", "red", "green"), lty = 1, bty = "n"
#' )
#' }
#'
#' @export
dmc <- function(x, gamma = 1, delta = 0, lambda = 1, log = FALSE) {
  # Input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }
  if (!is.numeric(lambda) || any(lambda <= 0)) {
    stop("'lambda' must be positive")
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("'log' must be a single logical value")
  }

  # Call C++ implementation
  .Call("_gkwdist_dmc",
    as.numeric(x),
    as.numeric(gamma),
    as.numeric(delta),
    as.numeric(lambda),
    as.logical(log),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 2. DISTRIBUTION FUNCTION (pmc)
# ----------------------------------------------------------------------------#

#' @title CDF of the McDonald (Mc)/Beta Power Distribution
#' @author Lopes, J. E.
#' @keywords distribution cumulative mcdonald
#'
#' @description
#' Computes the cumulative distribution function (CDF), \eqn{F(q) = P(X \le q)},
#' for the McDonald (Mc) distribution (also known as Beta Power) with
#' parameters \code{gamma} (\eqn{\gamma}), \code{delta} (\eqn{\delta}), and
#' \code{lambda} (\eqn{\lambda}). This distribution is defined on the interval
#' (0, 1) and is a special case of the Generalized Kumaraswamy (GKw)
#' distribution where \eqn{\alpha = 1} and \eqn{\beta = 1}.
#'
#' @param q Vector of quantiles (values generally between 0 and 1).
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param lambda Shape parameter \code{lambda} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}, otherwise, \eqn{P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of probabilities, \eqn{F(q)}, or their logarithms/complements
#'   depending on \code{lower.tail} and \code{log.p}. The length of the result
#'   is determined by the recycling rule applied to the arguments (\code{q},
#'   \code{gamma}, \code{delta}, \code{lambda}). Returns \code{0} (or \code{-Inf}
#'   if \code{log.p = TRUE}) for \code{q <= 0} and \code{1} (or \code{0} if
#'   \code{log.p = TRUE}) for \code{q >= 1}. Returns \code{NaN} for invalid
#'   parameters.
#'
#' @details
#' The McDonald (Mc) distribution is a special case of the five-parameter
#' Generalized Kumaraswamy (GKw) distribution (\code{\link{pgkw}}) obtained
#' by setting parameters \eqn{\alpha = 1} and \eqn{\beta = 1}.
#'
#' The CDF of the GKw distribution is \eqn{F_{GKw}(q) = I_{y(q)}(\gamma, \delta+1)},
#' where \eqn{y(q) = [1-(1-q^{\alpha})^{\beta}]^{\lambda}} and \eqn{I_x(a,b)}
#' is the regularized incomplete beta function (\code{\link[stats]{pbeta}}).
#' Setting \eqn{\alpha=1} and \eqn{\beta=1} simplifies \eqn{y(q)} to \eqn{q^\lambda},
#' yielding the Mc CDF:
#' \deqn{
#' F(q; \gamma, \delta, \lambda) = I_{q^\lambda}(\gamma, \delta+1)
#' }
#' This is evaluated using the \code{\link[stats]{pbeta}} function as
#' \code{pbeta(q^lambda, shape1 = gamma, shape2 = delta + 1)}.
#'
#' @references
#' McDonald, J. B. (1984). Some generalized functions for the size distribution
#' of income. *Econometrica*, 52(3), 647-663.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' @seealso
#' \code{\link{pgkw}} (parent distribution CDF),
#' \code{\link{dmc}}, \code{\link{qmc}}, \code{\link{rmc}} (other Mc functions),
#' \code{\link[stats]{pbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' q_vals <- c(0.2, 0.5, 0.8)
#' gamma_par <- 2.0
#' delta_par <- 1.5
#' lambda_par <- 1.0 # Equivalent to Beta(gamma, delta+1)
#'
#' # Calculate CDF P(X <= q) using pmc
#' probs <- pmc(q_vals, gamma_par, delta_par, lambda_par)
#' print(probs)
#' # Compare with Beta CDF
#' print(stats::pbeta(q_vals, shape1 = gamma_par, shape2 = delta_par + 1))
#'
#' # Calculate upper tail P(X > q)
#' probs_upper <- pmc(q_vals, gamma_par, delta_par, lambda_par,
#'   lower.tail = FALSE
#' )
#' print(probs_upper)
#' # Check: probs + probs_upper should be 1
#' print(probs + probs_upper)
#'
#' # Calculate log CDF
#' logs <- pmc(q_vals, gamma_par, delta_par, lambda_par, log.p = TRUE)
#' print(logs)
#' # Check: should match log(probs)
#' print(log(probs))
#'
#' # Compare with pgkw setting alpha = 1, beta = 1
#' probs_gkw <- pgkw(q_vals,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = lambda_par
#' )
#' print(paste("Max difference:", max(abs(probs - probs_gkw)))) # Should be near zero
#'
#' # Plot the CDF for different lambda values
#' curve_q <- seq(0.01, 0.99, length.out = 200)
#' curve_p1 <- pmc(curve_q, gamma = 2, delta = 3, lambda = 0.5)
#' curve_p2 <- pmc(curve_q, gamma = 2, delta = 3, lambda = 1.0) # Beta(2, 4)
#' curve_p3 <- pmc(curve_q, gamma = 2, delta = 3, lambda = 2.0)
#'
#' plot(curve_q, curve_p2,
#'   type = "l", main = "Mc/Beta Power CDF (gamma=2, delta=3)",
#'   xlab = "q", ylab = "F(q)", col = "red", ylim = c(0, 1)
#' )
#' lines(curve_q, curve_p1, col = "blue")
#' lines(curve_q, curve_p3, col = "green")
#' legend("bottomright",
#'   legend = c("lambda=0.5", "lambda=1.0 (Beta)", "lambda=2.0"),
#'   col = c("blue", "red", "green"), lty = 1, bty = "n"
#' )
#' }
#'
#' @export
pmc <- function(q, gamma = 1, delta = 0, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(q)) stop("'q' must be numeric")
  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }
  if (!is.numeric(lambda) || any(lambda <= 0)) {
    stop("'lambda' must be positive")
  }
  if (!is.logical(lower.tail) || length(lower.tail) != 1) {
    stop("'lower.tail' must be a single logical value")
  }
  if (!is.logical(log.p) || length(log.p) != 1) {
    stop("'log.p' must be a single logical value")
  }

  # Call C++ implementation
  .Call("_gkwdist_pmc",
    as.numeric(q),
    as.numeric(gamma),
    as.numeric(delta),
    as.numeric(lambda),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 3. QUANTILE FUNCTION (qmc)
# ----------------------------------------------------------------------------#

#' @title Quantile Function of the McDonald (Mc)/Beta Power Distribution
#' @author Lopes, J. E.
#' @keywords distribution quantile mcdonald
#'
#' @description
#' Computes the quantile function (inverse CDF) for the McDonald (Mc) distribution
#' (also known as Beta Power) with parameters \code{gamma} (\eqn{\gamma}),
#' \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}). It finds the
#' value \code{q} such that \eqn{P(X \le q) = p}. This distribution is a special
#' case of the Generalized Kumaraswamy (GKw) distribution where \eqn{\alpha = 1}
#' and \eqn{\beta = 1}.
#'
#' @param p Vector of probabilities (values between 0 and 1).
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param lambda Shape parameter \code{lambda} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{p = P(X \le q)},
#'   otherwise, probabilities are \eqn{p = P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of quantiles corresponding to the given probabilities \code{p}.
#'   The length of the result is determined by the recycling rule applied to
#'   the arguments (\code{p}, \code{gamma}, \code{delta}, \code{lambda}).
#'   Returns:
#'   \itemize{
#'     \item \code{0} for \code{p = 0} (or \code{p = -Inf} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{1} for \code{p = 1} (or \code{p = 0} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{NaN} for \code{p < 0} or \code{p > 1} (or corresponding log scale).
#'     \item \code{NaN} for invalid parameters (e.g., \code{gamma <= 0},
#'           \code{delta < 0}, \code{lambda <= 0}).
#'   }
#'   Boundary return values are adjusted accordingly for \code{lower.tail = FALSE}.
#'
#' @details
#' The quantile function \eqn{Q(p)} is the inverse of the CDF \eqn{F(q)}. The CDF
#' for the Mc (\eqn{\alpha=1, \beta=1}) distribution is \eqn{F(q) = I_{q^\lambda}(\gamma, \delta+1)},
#' where \eqn{I_z(a,b)} is the regularized incomplete beta function (see \code{\link{pmc}}).
#'
#' To find the quantile \eqn{q}, we first invert the Beta function part: let
#' \eqn{y = I^{-1}_{p}(\gamma, \delta+1)}, where \eqn{I^{-1}_p(a,b)} is the
#' inverse computed via \code{\link[stats]{qbeta}}. We then solve \eqn{q^\lambda = y}
#' for \eqn{q}, yielding the quantile function:
#' \deqn{
#' Q(p) = \left[ I^{-1}_{p}(\gamma, \delta+1) \right]^{1/\lambda}
#' }
#' The function uses this formula, calculating \eqn{I^{-1}_{p}(\gamma, \delta+1)}
#' via \code{qbeta(p, gamma, delta + 1, ...)} while respecting the
#' \code{lower.tail} and \code{log.p} arguments. This is equivalent to the general
#' GKw quantile function (\code{\link{qgkw}}) evaluated with \eqn{\alpha=1, \beta=1}.
#'
#' @references
#' McDonald, J. B. (1984). Some generalized functions for the size distribution
#' of income. *Econometrica*, 52(3), 647-663.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' @seealso
#' \code{\link{qgkw}} (parent distribution quantile function),
#' \code{\link{dmc}}, \code{\link{pmc}}, \code{\link{rmc}} (other Mc functions),
#' \code{\link[stats]{qbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' p_vals <- c(0.1, 0.5, 0.9)
#' gamma_par <- 2.0
#' delta_par <- 1.5
#' lambda_par <- 1.0 # Equivalent to Beta(gamma, delta+1)
#'
#' # Calculate quantiles using qmc
#' quantiles <- qmc(p_vals, gamma_par, delta_par, lambda_par)
#' print(quantiles)
#' # Compare with Beta quantiles
#' print(stats::qbeta(p_vals, shape1 = gamma_par, shape2 = delta_par + 1))
#'
#' # Calculate quantiles for upper tail probabilities P(X > q) = p
#' quantiles_upper <- qmc(p_vals, gamma_par, delta_par, lambda_par,
#'   lower.tail = FALSE
#' )
#' print(quantiles_upper)
#' # Check: qmc(p, ..., lt=F) == qmc(1-p, ..., lt=T)
#' print(qmc(1 - p_vals, gamma_par, delta_par, lambda_par))
#'
#' # Calculate quantiles from log probabilities
#' log.p_vals <- log(p_vals)
#' quantiles_logp <- qmc(log.p_vals, gamma_par, delta_par, lambda_par, log.p = TRUE)
#' print(quantiles_logp)
#' # Check: should match original quantiles
#' print(quantiles)
#'
#' # Compare with qgkw setting alpha = 1, beta = 1
#' quantiles_gkw <- qgkw(p_vals,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = lambda_par
#' )
#' print(paste("Max difference:", max(abs(quantiles - quantiles_gkw)))) # Should be near zero
#'
#' # Verify inverse relationship with pmc
#' p_check <- 0.75
#' q_calc <- qmc(p_check, gamma_par, delta_par, lambda_par) # Use lambda != 1
#' p_recalc <- pmc(q_calc, gamma_par, delta_par, lambda_par)
#' print(paste("Original p:", p_check, " Recalculated p:", p_recalc))
#' # abs(p_check - p_recalc) < 1e-9 # Should be TRUE
#'
#' # Boundary conditions
#' print(qmc(c(0, 1), gamma_par, delta_par, lambda_par)) # Should be 0, 1
#' print(qmc(c(-Inf, 0), gamma_par, delta_par, lambda_par, log.p = TRUE)) # Should be 0, 1
#' }
#'
#' @export
qmc <- function(p, gamma = 1, delta = 0, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(p)) stop("'p' must be numeric")
  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }
  if (!is.numeric(lambda) || any(lambda <= 0)) {
    stop("'lambda' must be positive")
  }
  if (!is.logical(lower.tail) || length(lower.tail) != 1) {
    stop("'lower.tail' must be a single logical value")
  }
  if (!is.logical(log.p) || length(log.p) != 1) {
    stop("'log.p' must be a single logical value")
  }

  # Additional validation for probabilities
  if (!log.p && any(p < 0 | p > 1, na.rm = TRUE)) {
    warning("'p' values outside [0, 1] will produce NaN")
  }

  # Call C++ implementation
  .Call("_gkwdist_qmc",
    as.numeric(p),
    as.numeric(gamma),
    as.numeric(delta),
    as.numeric(lambda),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 4. RANDOM GENERATION (rmc)
# ----------------------------------------------------------------------------#

#' @title Random Number Generation for the McDonald (Mc)/Beta Power Distribution
#' @author Lopes, J. E.
#' @keywords distribution random mcdonald
#'
#' @description
#' Generates random deviates from the McDonald (Mc) distribution (also known as
#' Beta Power) with parameters \code{gamma} (\eqn{\gamma}), \code{delta}
#' (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}). This distribution is a
#' special case of the Generalized Kumaraswamy (GKw) distribution where
#' \eqn{\alpha = 1} and \eqn{\beta = 1}.
#'
#' @param n Number of observations. If \code{length(n) > 1}, the length is
#'   taken to be the number required. Must be a non-negative integer.
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param lambda Shape parameter \code{lambda} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#'
#' @return A vector of length \code{n} containing random deviates from the Mc
#'   distribution, with values in (0, 1). The length of the result is determined
#'   by \code{n} and the recycling rule applied to the parameters (\code{gamma},
#'   \code{delta}, \code{lambda}). Returns \code{NaN} if parameters
#'   are invalid (e.g., \code{gamma <= 0}, \code{delta < 0}, \code{lambda <= 0}).
#'
#' @details
#' The generation method uses the relationship between the GKw distribution and the
#' Beta distribution. The general procedure for GKw (\code{\link{rgkw}}) is:
#' If \eqn{W \sim \mathrm{Beta}(\gamma, \delta+1)}, then
#' \eqn{X = \{1 - [1 - W^{1/\lambda}]^{1/\beta}\}^{1/\alpha}} follows the
#' GKw(\eqn{\alpha, \beta, \gamma, \delta, \lambda}) distribution.
#'
#' For the Mc distribution, \eqn{\alpha=1} and \eqn{\beta=1}. Therefore, the
#' algorithm simplifies significantly:
#' \enumerate{
#'   \item Generate \eqn{U \sim \mathrm{Beta}(\gamma, \delta+1)} using
#'         \code{\link[stats]{rbeta}}.
#'   \item Compute the Mc variate \eqn{X = U^{1/\lambda}}.
#' }
#' This procedure is implemented efficiently, handling parameter recycling as needed.
#'
#' @references
#' McDonald, J. B. (1984). Some generalized functions for the size distribution
#' of income. *Econometrica*, 52(3), 647-663.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer-Verlag.
#' (General methods for random variate generation).
#'
#' @seealso
#' \code{\link{rgkw}} (parent distribution random generation),
#' \code{\link{dmc}}, \code{\link{pmc}}, \code{\link{qmc}} (other Mc functions),
#' \code{\link[stats]{rbeta}}
#'
#' @examples
#' \donttest{
#' set.seed(2028) # for reproducibility
#'
#' # Generate 1000 random values from a specific Mc distribution
#' gamma_par <- 2.0
#' delta_par <- 1.5
#' lambda_par <- 1.0 # Equivalent to Beta(gamma, delta+1)
#'
#' x_sample_mc <- rmc(1000,
#'   gamma = gamma_par, delta = delta_par,
#'   lambda = lambda_par
#' )
#' summary(x_sample_mc)
#'
#' # Histogram of generated values compared to theoretical density
#' hist(x_sample_mc,
#'   breaks = 30, freq = FALSE, # freq=FALSE for density
#'   main = "Histogram of Mc Sample (Beta Case)", xlab = "x"
#' )
#' curve(dmc(x, gamma = gamma_par, delta = delta_par, lambda = lambda_par),
#'   add = TRUE, col = "red", lwd = 2, n = 201
#' )
#' curve(stats::dbeta(x, gamma_par, delta_par + 1), add = TRUE, col = "blue", lty = 2)
#' legend("topright",
#'   legend = c("Theoretical Mc PDF", "Theoretical Beta PDF"),
#'   col = c("red", "blue"), lwd = c(2, 1), lty = c(1, 2), bty = "n"
#' )
#'
#' # Comparing empirical and theoretical quantiles (Q-Q plot)
#' lambda_par_qq <- 0.7 # Use lambda != 1 for non-Beta case
#' x_sample_mc_qq <- rmc(1000,
#'   gamma = gamma_par, delta = delta_par,
#'   lambda = lambda_par_qq
#' )
#' prob_points <- seq(0.01, 0.99, by = 0.01)
#' theo_quantiles <- qmc(prob_points,
#'   gamma = gamma_par, delta = delta_par,
#'   lambda = lambda_par_qq
#' )
#' emp_quantiles <- quantile(x_sample_mc_qq, prob_points, type = 7)
#'
#' plot(theo_quantiles, emp_quantiles,
#'   pch = 16, cex = 0.8,
#'   main = "Q-Q Plot for Mc Distribution",
#'   xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles (n=1000)"
#' )
#' abline(a = 0, b = 1, col = "blue", lty = 2)
#'
#' # Compare summary stats with rgkw(..., alpha=1, beta=1, ...)
#' # Note: individual values will differ due to randomness
#' x_sample_gkw <- rgkw(1000,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = lambda_par_qq
#' )
#' print("Summary stats for rmc sample:")
#' print(summary(x_sample_mc_qq))
#' print("Summary stats for rgkw(alpha=1, beta=1) sample:")
#' print(summary(x_sample_gkw)) # Should be similar
#' }
#'
#' @export
rmc <- function(n, gamma = 1, delta = 0, lambda = 1) {
  # Input validation
  if (length(n) > 1) n <- length(n)
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("'n' must be a positive integer")
  }
  n <- as.integer(n)

  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }
  if (!is.numeric(lambda) || any(lambda <= 0)) {
    stop("'lambda' must be positive")
  }

  # Call C++ implementation
  .Call("_gkwdist_rmc",
    as.integer(n),
    as.numeric(gamma),
    as.numeric(delta),
    as.numeric(lambda),
    PACKAGE = "gkwdist"
  )
}


# ============================================================================#
# MAXIMUM LIKELIHOOD ESTIMATION FUNCTIONS
# ============================================================================#

# ----------------------------------------------------------------------------#
# 5. NEGATIVE LOG-LIKELIHOOD (llmc)
# ----------------------------------------------------------------------------#

#' @title Negative Log-Likelihood for the McDonald (Mc)/Beta Power Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize mcdonald
#'
#' @description
#' Computes the negative log-likelihood function for the McDonald (Mc)
#' distribution (also known as Beta Power) with parameters \code{gamma}
#' (\eqn{\gamma}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}),
#' given a vector of observations. This distribution is the special case of the
#' Generalized Kumaraswamy (GKw) distribution where \eqn{\alpha = 1} and
#' \eqn{\beta = 1}. This function is suitable for maximum likelihood estimation.
#'
#' @param par A numeric vector of length 3 containing the distribution parameters
#'   in the order: \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}),
#'   \code{lambda} (\eqn{\lambda > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a single \code{double} value representing the negative
#'   log-likelihood (\eqn{-\ell(\theta|\mathbf{x})}). Returns \code{Inf}
#'   if any parameter values in \code{par} are invalid according to their
#'   constraints, or if any value in \code{data} is not in the interval (0, 1).
#'
#' @details
#' The McDonald (Mc) distribution is the GKw distribution (\code{\link{dmc}})
#' with \eqn{\alpha=1} and \eqn{\beta=1}. Its probability density function (PDF) is:
#' \deqn{
#' f(x | \theta) = \frac{\lambda}{B(\gamma,\delta+1)} x^{\gamma \lambda - 1} (1 - x^\lambda)^\delta
#' }
#' for \eqn{0 < x < 1}, \eqn{\theta = (\gamma, \delta, \lambda)}, and \eqn{B(a,b)}
#' is the Beta function (\code{\link[base]{beta}}).
#' The log-likelihood function \eqn{\ell(\theta | \mathbf{x})} for a sample
#' \eqn{\mathbf{x} = (x_1, \dots, x_n)} is \eqn{\sum_{i=1}^n \ln f(x_i | \theta)}:
#' \deqn{
#' \ell(\theta | \mathbf{x}) = n[\ln(\lambda) - \ln B(\gamma, \delta+1)]
#' + \sum_{i=1}^{n} [(\gamma\lambda - 1)\ln(x_i) + \delta\ln(1 - x_i^\lambda)]
#' }
#' This function computes and returns the *negative* log-likelihood, \eqn{-\ell(\theta|\mathbf{x})},
#' suitable for minimization using optimization routines like \code{\link[stats]{optim}}.
#' Numerical stability is maintained, including using the log-gamma function
#' (\code{\link[base]{lgamma}}) for the Beta function term.
#'
#' @references
#' McDonald, J. B. (1984). Some generalized functions for the size distribution
#' of income. *Econometrica*, 52(3), 647-663.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' @seealso
#' \code{\link{llgkw}} (parent distribution negative log-likelihood),
#' \code{\link{dmc}}, \code{\link{pmc}}, \code{\link{qmc}}, \code{\link{rmc}},
#' \code{grmc} (gradient, if available),
#' \code{hsmc} (Hessian, if available),
#' \code{\link[stats]{optim}}, \code{\link[base]{lbeta}}
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Log-Likelihood Evaluation
#'
#' # Generate sample data with more stable parameters
#' set.seed(123)
#' n <- 1000
#' true_params <- c(gamma = 2.0, delta = 2.5, lambda = 1.5)
#' data <- rmc(n,
#'   gamma = true_params[1], delta = true_params[2],
#'   lambda = true_params[3]
#' )
#'
#' # Evaluate negative log-likelihood at true parameters
#' nll_true <- llmc(par = true_params, data = data)
#' cat("Negative log-likelihood at true parameters:", nll_true, "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 2.0, 1.0),
#'   c(2.0, 2.5, 1.5),
#'   c(2.5, 3.0, 2.0)
#' )
#'
#' nll_values <- apply(test_params, 1, function(p) llmc(p, data))
#' results <- data.frame(
#'   Gamma = test_params[, 1],
#'   Delta = test_params[, 2],
#'   Lambda = test_params[, 3],
#'   NegLogLik = nll_values
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Maximum Likelihood Estimation
#'
#' # Optimization using BFGS with analytical gradient
#' fit <- optim(
#'   par = c(1.5, 2.0, 1.0),
#'   fn = llmc,
#'   gr = grmc,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("gamma", "delta", "lambda")
#' se <- sqrt(diag(solve(fit$hessian)))
#'
#' results <- data.frame(
#'   Parameter = c("gamma", "delta", "lambda"),
#'   True = true_params,
#'   MLE = mle,
#'   SE = se,
#'   CI_Lower = mle - 1.96 * se,
#'   CI_Upper = mle + 1.96 * se
#' )
#' print(results, digits = 4)
#'
#' cat("\nNegative log-likelihood at MLE:", fit$value, "\n")
#' cat("AIC:", 2 * fit$value + 2 * length(mle), "\n")
#' cat("BIC:", 2 * fit$value + length(mle) * log(n), "\n")
#'
#'
#' ## Example 3: Comparing Optimization Methods
#'
#' methods <- c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG")
#' start_params <- c(1.5, 2.0, 1.0)
#'
#' comparison <- data.frame(
#'   Method = character(),
#'   Gamma = numeric(),
#'   Delta = numeric(),
#'   Lambda = numeric(),
#'   NegLogLik = numeric(),
#'   Convergence = integer(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (method in methods) {
#'   if (method %in% c("BFGS", "CG")) {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llmc,
#'       gr = grmc,
#'       data = data,
#'       method = method
#'     )
#'   } else if (method == "L-BFGS-B") {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llmc,
#'       gr = grmc,
#'       data = data,
#'       method = method,
#'       lower = c(0.01, 0.01, 0.01),
#'       upper = c(100, 100, 100)
#'     )
#'   } else {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llmc,
#'       data = data,
#'       method = method
#'     )
#'   }
#'
#'   comparison <- rbind(comparison, data.frame(
#'     Method = method,
#'     Gamma = fit_temp$par[1],
#'     Delta = fit_temp$par[2],
#'     Lambda = fit_temp$par[3],
#'     NegLogLik = fit_temp$value,
#'     Convergence = fit_temp$convergence,
#'     stringsAsFactors = FALSE
#'   ))
#' }
#'
#' print(comparison, digits = 4, row.names = FALSE)
#'
#'
#' ## Example 4: Likelihood Ratio Test
#'
#' # Test H0: lambda = 1.5 vs H1: lambda free
#' loglik_full <- -fit$value
#'
#' restricted_ll <- function(params_restricted, data, lambda_fixed) {
#'   llmc(par = c(
#'     params_restricted[1], params_restricted[2],
#'     lambda_fixed
#'   ), data = data)
#' }
#'
#' fit_restricted <- optim(
#'   par = c(mle[1], mle[2]),
#'   fn = restricted_ll,
#'   data = data,
#'   lambda_fixed = 1.5,
#'   method = "BFGS"
#' )
#'
#' loglik_restricted <- -fit_restricted$value
#' lr_stat <- 2 * (loglik_full - loglik_restricted)
#' p_value <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
#'
#' cat("LR Statistic:", round(lr_stat, 4), "\n")
#' cat("P-value:", format.pval(p_value, digits = 4), "\n")
#'
#'
#' ## Example 5: Univariate Profile Likelihoods
#'
#' # Profile for gamma
#' gamma_grid <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = 50)
#' gamma_grid <- gamma_grid[gamma_grid > 0]
#' profile_ll_gamma <- numeric(length(gamma_grid))
#'
#' for (i in seq_along(gamma_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-1],
#'     fn = function(p) llmc(c(gamma_grid[i], p), data),
#'     method = "BFGS"
#'   )
#'   profile_ll_gamma[i] <- -profile_fit$value
#' }
#'
#' # Profile for delta
#' delta_grid <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = 50)
#' delta_grid <- delta_grid[delta_grid > 0]
#' profile_ll_delta <- numeric(length(delta_grid))
#'
#' for (i in seq_along(delta_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-2],
#'     fn = function(p) llmc(c(p[1], delta_grid[i], p[2]), data),
#'     method = "BFGS"
#'   )
#'   profile_ll_delta[i] <- -profile_fit$value
#' }
#'
#' # Profile for lambda
#' lambda_grid <- seq(mle[3] - 1.5, mle[3] + 1.5, length.out = 50)
#' lambda_grid <- lambda_grid[lambda_grid > 0]
#' profile_ll_lambda <- numeric(length(lambda_grid))
#'
#' for (i in seq_along(lambda_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-3],
#'     fn = function(p) llmc(c(p[1], p[2], lambda_grid[i]), data),
#'     method = "BFGS"
#'   )
#'   profile_ll_lambda[i] <- -profile_fit$value
#' }
#'
#' # 95% confidence threshold
#' chi_crit <- qchisq(0.95, df = 1)
#' threshold <- max(profile_ll_gamma) - chi_crit / 2
#'
#' # Plot all profiles
#'
#' plot(gamma_grid, profile_ll_gamma,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", gamma)), las = 1
#' )
#' abline(v = mle[1], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[1], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3, lwd = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#808080"),
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.8
#' )
#' grid(col = "gray90")
#'
#' plot(delta_grid, profile_ll_delta,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(delta), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", delta)), las = 1
#' )
#' abline(v = mle[2], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[2], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3, lwd = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#808080"),
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.8
#' )
#' grid(col = "gray90")
#'
#' plot(lambda_grid, profile_ll_lambda,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(lambda), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", lambda)), las = 1
#' )
#' abline(v = mle[3], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[3], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3, lwd = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#808080"),
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.8
#' )
#' grid(col = "gray90")
#'
#' ## Example 6: 2D Log-Likelihood Surfaces (All pairs side by side)
#'
#' # Create 2D grids with wider range (±1.5)
#' gamma_2d <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = round(n / 25))
#' delta_2d <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = round(n / 25))
#' lambda_2d <- seq(mle[3] - 1.5, mle[3] + 1.5, length.out = round(n / 25))
#'
#' gamma_2d <- gamma_2d[gamma_2d > 0]
#' delta_2d <- delta_2d[delta_2d > 0]
#' lambda_2d <- lambda_2d[lambda_2d > 0]
#'
#' # Compute all log-likelihood surfaces
#' ll_surface_gd <- matrix(NA, nrow = length(gamma_2d), ncol = length(delta_2d))
#' ll_surface_gl <- matrix(NA, nrow = length(gamma_2d), ncol = length(lambda_2d))
#' ll_surface_dl <- matrix(NA, nrow = length(delta_2d), ncol = length(lambda_2d))
#'
#' # Gamma vs Delta
#' for (i in seq_along(gamma_2d)) {
#'   for (j in seq_along(delta_2d)) {
#'     ll_surface_gd[i, j] <- -llmc(c(gamma_2d[i], delta_2d[j], mle[3]), data)
#'   }
#' }
#'
#' # Gamma vs Lambda
#' for (i in seq_along(gamma_2d)) {
#'   for (j in seq_along(lambda_2d)) {
#'     ll_surface_gl[i, j] <- -llmc(c(gamma_2d[i], mle[2], lambda_2d[j]), data)
#'   }
#' }
#'
#' # Delta vs Lambda
#' for (i in seq_along(delta_2d)) {
#'   for (j in seq_along(lambda_2d)) {
#'     ll_surface_dl[i, j] <- -llmc(c(mle[1], delta_2d[i], lambda_2d[j]), data)
#'   }
#' }
#'
#' # Confidence region levels
#' max_ll_gd <- max(ll_surface_gd, na.rm = TRUE)
#' max_ll_gl <- max(ll_surface_gl, na.rm = TRUE)
#' max_ll_dl <- max(ll_surface_dl, na.rm = TRUE)
#'
#' levels_95_gd <- max_ll_gd - qchisq(0.95, df = 2) / 2
#' levels_95_gl <- max_ll_gl - qchisq(0.95, df = 2) / 2
#' levels_95_dl <- max_ll_dl - qchisq(0.95, df = 2) / 2
#'
#' # Plot
#'
#' # Gamma vs Delta
#' contour(gamma_2d, delta_2d, ll_surface_gd,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Gamma vs Delta", las = 1,
#'   levels = seq(min(ll_surface_gd, na.rm = TRUE), max_ll_gd, length.out = 20),
#'   col = "#2E4057", lwd = 1
#' )
#' contour(gamma_2d, delta_2d, ll_surface_gd,
#'   levels = levels_95_gd, col = "#FF6347", lwd = 2.5, lty = 1, add = TRUE
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Gamma vs Lambda
#' contour(gamma_2d, lambda_2d, ll_surface_gl,
#'   xlab = expression(gamma), ylab = expression(lambda),
#'   main = "Gamma vs Lambda", las = 1,
#'   levels = seq(min(ll_surface_gl, na.rm = TRUE), max_ll_gl, length.out = 20),
#'   col = "#2E4057", lwd = 1
#' )
#' contour(gamma_2d, lambda_2d, ll_surface_gl,
#'   levels = levels_95_gl, col = "#FF6347", lwd = 2.5, lty = 1, add = TRUE
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Delta vs Lambda
#' contour(delta_2d, lambda_2d, ll_surface_dl,
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "Delta vs Lambda", las = 1,
#'   levels = seq(min(ll_surface_dl, na.rm = TRUE), max_ll_dl, length.out = 20),
#'   col = "#2E4057", lwd = 1
#' )
#' contour(delta_2d, lambda_2d, ll_surface_dl,
#'   levels = levels_95_dl, col = "#FF6347", lwd = 2.5, lty = 1, add = TRUE
#' )
#' points(mle[2], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' legend("topright",
#'   legend = c("MLE", "True", "95% CR"),
#'   col = c("#8B0000", "#006400", "#FF6347"),
#'   pch = c(19, 17, NA),
#'   lty = c(NA, NA, 1),
#'   lwd = c(NA, NA, 2.5),
#'   bty = "n", cex = 0.8
#' )
#' }
#'
#' @export
llmc <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 3) {
    stop("'par' must be a numeric vector of length 3")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_llmc",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 6. GRADIENT (grmc)
# ----------------------------------------------------------------------------#

#' @title Gradient of the Negative Log-Likelihood for the McDonald (Mc)/Beta Power Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize gradient mcdonald
#'
#' @description
#' Computes the gradient vector (vector of first partial derivatives) of the
#' negative log-likelihood function for the McDonald (Mc) distribution (also
#' known as Beta Power) with parameters \code{gamma} (\eqn{\gamma}), \code{delta}
#' (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}). This distribution is the
#' special case of the Generalized Kumaraswamy (GKw) distribution where
#' \eqn{\alpha = 1} and \eqn{\beta = 1}. The gradient is useful for optimization.
#'
#' @param par A numeric vector of length 3 containing the distribution parameters
#'   in the order: \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}),
#'   \code{lambda} (\eqn{\lambda > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a numeric vector of length 3 containing the partial derivatives
#'   of the negative log-likelihood function \eqn{-\ell(\theta | \mathbf{x})} with
#'   respect to each parameter:
#'   \eqn{(-\partial \ell/\partial \gamma, -\partial \ell/\partial \delta, -\partial \ell/\partial \lambda)}.
#'   Returns a vector of \code{NaN} if any parameter values are invalid according
#'   to their constraints, or if any value in \code{data} is not in the
#'   interval (0, 1).
#'
#' @details
#' The components of the gradient vector of the negative log-likelihood
#' (\eqn{-\nabla \ell(\theta | \mathbf{x})}) for the Mc (\eqn{\alpha=1, \beta=1})
#' model are:
#'
#' \deqn{
#' -\frac{\partial \ell}{\partial \gamma} = n[\psi(\gamma+\delta+1) - \psi(\gamma)] -
#' \lambda\sum_{i=1}^{n}\ln(x_i)
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \delta} = n[\psi(\gamma+\delta+1) - \psi(\delta+1)] -
#' \sum_{i=1}^{n}\ln(1-x_i^{\lambda})
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \lambda} = -\frac{n}{\lambda} - \gamma\sum_{i=1}^{n}\ln(x_i) +
#' \delta\sum_{i=1}^{n}\frac{x_i^{\lambda}\ln(x_i)}{1-x_i^{\lambda}}
#' }
#'
#' where \eqn{\psi(\cdot)} is the digamma function (\code{\link[base]{digamma}}).
#' These formulas represent the derivatives of \eqn{-\ell(\theta)}, consistent with
#' minimizing the negative log-likelihood. They correspond to the relevant components
#' of the general GKw gradient (\code{\link{grgkw}}) evaluated at \eqn{\alpha=1, \beta=1}.
#'
#' @references
#' McDonald, J. B. (1984). Some generalized functions for the size distribution
#' of income. *Econometrica*, 52(3), 647-663.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' (Note: Specific gradient formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{grgkw}} (parent distribution gradient),
#' \code{\link{llmc}} (negative log-likelihood for Mc),
#' \code{hsmc} (Hessian for Mc, if available),
#' \code{\link{dmc}} (density for Mc),
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{grad}} (for numerical gradient comparison),
#' \code{\link[base]{digamma}}.
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Examples
#'
#' # Generate sample data with more stable parameters
#' set.seed(123)
#' n <- 1000
#' true_params <- c(gamma = 2.0, delta = 2.5, lambda = 1.5)
#' data <- rmc(n,
#'   gamma = true_params[1], delta = true_params[2],
#'   lambda = true_params[3]
#' )
#'
#' # Evaluate Hessian at true parameters
#' hess_true <- hsmc(par = true_params, data = data)
#' cat("Hessian matrix at true parameters:\n")
#' print(hess_true, digits = 4)
#'
#' # Check symmetry
#' cat(
#'   "\nSymmetry check (max |H - H^T|):",
#'   max(abs(hess_true - t(hess_true))), "\n"
#' )
#'
#'
#' ## Example 2: Hessian Properties at MLE
#'
#' # Fit model
#' fit <- optim(
#'   par = c(1.5, 2.0, 1.0),
#'   fn = llmc,
#'   gr = grmc,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("gamma", "delta", "lambda")
#'
#' # Hessian at MLE
#' hessian_at_mle <- hsmc(par = mle, data = data)
#' cat("\nHessian at MLE:\n")
#' print(hessian_at_mle, digits = 4)
#'
#' # Compare with optim's numerical Hessian
#' cat("\nComparison with optim Hessian:\n")
#' cat(
#'   "Max absolute difference:",
#'   max(abs(hessian_at_mle - fit$hessian)), "\n"
#' )
#'
#' # Eigenvalue analysis
#' eigenvals <- eigen(hessian_at_mle, only.values = TRUE)$values
#' cat("\nEigenvalues:\n")
#' print(eigenvals)
#'
#' cat("\nPositive definite:", all(eigenvals > 0), "\n")
#' cat("Condition number:", max(eigenvals) / min(eigenvals), "\n")
#'
#'
#' ## Example 3: Standard Errors and Confidence Intervals
#'
#' # Observed information matrix
#' obs_info <- hessian_at_mle
#'
#' # Variance-covariance matrix
#' vcov_matrix <- solve(obs_info)
#' cat("\nVariance-Covariance Matrix:\n")
#' print(vcov_matrix, digits = 6)
#'
#' # Standard errors
#' se <- sqrt(diag(vcov_matrix))
#' names(se) <- c("gamma", "delta", "lambda")
#'
#' # Correlation matrix
#' corr_matrix <- cov2cor(vcov_matrix)
#' cat("\nCorrelation Matrix:\n")
#' print(corr_matrix, digits = 4)
#'
#' # Confidence intervals
#' z_crit <- qnorm(0.975)
#' results <- data.frame(
#'   Parameter = c("gamma", "delta", "lambda"),
#'   True = true_params,
#'   MLE = mle,
#'   SE = se,
#'   CI_Lower = mle - z_crit * se,
#'   CI_Upper = mle + z_crit * se
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 4: Determinant and Trace Analysis
#'
#' # Compute at different points
#' test_params <- rbind(
#'   c(1.5, 2.0, 1.0),
#'   c(2.0, 2.5, 1.5),
#'   mle,
#'   c(2.5, 3.0, 2.0)
#' )
#'
#' hess_properties <- data.frame(
#'   Gamma = numeric(),
#'   Delta = numeric(),
#'   Lambda = numeric(),
#'   Determinant = numeric(),
#'   Trace = numeric(),
#'   Min_Eigenval = numeric(),
#'   Max_Eigenval = numeric(),
#'   Cond_Number = numeric(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (i in 1:nrow(test_params)) {
#'   H <- hsmc(par = test_params[i, ], data = data)
#'   eigs <- eigen(H, only.values = TRUE)$values
#'
#'   hess_properties <- rbind(hess_properties, data.frame(
#'     Gamma = test_params[i, 1],
#'     Delta = test_params[i, 2],
#'     Lambda = test_params[i, 3],
#'     Determinant = det(H),
#'     Trace = sum(diag(H)),
#'     Min_Eigenval = min(eigs),
#'     Max_Eigenval = max(eigs),
#'     Cond_Number = max(eigs) / min(eigs)
#'   ))
#' }
#'
#' cat("\nHessian Properties at Different Points:\n")
#' print(hess_properties, digits = 4, row.names = FALSE)
#'
#'
#' ## Example 5: Curvature Visualization (All pairs side by side)
#'
#' # Create grids around MLE with wider range (±1.5)
#' gamma_grid <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = 25)
#' delta_grid <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = 25)
#' lambda_grid <- seq(mle[3] - 1.5, mle[3] + 1.5, length.out = 25)
#'
#' gamma_grid <- gamma_grid[gamma_grid > 0]
#' delta_grid <- delta_grid[delta_grid > 0]
#' lambda_grid <- lambda_grid[lambda_grid > 0]
#'
#' # Compute curvature measures for all pairs
#' determinant_surface_gd <- matrix(NA, nrow = length(gamma_grid), ncol = length(delta_grid))
#' trace_surface_gd <- matrix(NA, nrow = length(gamma_grid), ncol = length(delta_grid))
#'
#' determinant_surface_gl <- matrix(NA, nrow = length(gamma_grid), ncol = length(lambda_grid))
#' trace_surface_gl <- matrix(NA, nrow = length(gamma_grid), ncol = length(lambda_grid))
#'
#' determinant_surface_dl <- matrix(NA, nrow = length(delta_grid), ncol = length(lambda_grid))
#' trace_surface_dl <- matrix(NA, nrow = length(delta_grid), ncol = length(lambda_grid))
#'
#' # Gamma vs Delta
#' for (i in seq_along(gamma_grid)) {
#'   for (j in seq_along(delta_grid)) {
#'     H <- hsmc(c(gamma_grid[i], delta_grid[j], mle[3]), data)
#'     determinant_surface_gd[i, j] <- det(H)
#'     trace_surface_gd[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Gamma vs Lambda
#' for (i in seq_along(gamma_grid)) {
#'   for (j in seq_along(lambda_grid)) {
#'     H <- hsmc(c(gamma_grid[i], mle[2], lambda_grid[j]), data)
#'     determinant_surface_gl[i, j] <- det(H)
#'     trace_surface_gl[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Delta vs Lambda
#' for (i in seq_along(delta_grid)) {
#'   for (j in seq_along(lambda_grid)) {
#'     H <- hsmc(c(mle[1], delta_grid[i], lambda_grid[j]), data)
#'     determinant_surface_dl[i, j] <- det(H)
#'     trace_surface_dl[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Plot
#'
#' # Determinant plots
#' contour(gamma_grid, delta_grid, determinant_surface_gd,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Determinant: Gamma vs Delta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(gamma_grid, lambda_grid, determinant_surface_gl,
#'   xlab = expression(gamma), ylab = expression(lambda),
#'   main = "Determinant: Gamma vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(delta_grid, lambda_grid, determinant_surface_dl,
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "Determinant: Delta vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[2], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Trace plots
#' contour(gamma_grid, delta_grid, trace_surface_gd,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Trace: Gamma vs Delta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(gamma_grid, lambda_grid, trace_surface_gl,
#'   xlab = expression(gamma), ylab = expression(lambda),
#'   main = "Trace: Gamma vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(delta_grid, lambda_grid, trace_surface_dl,
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "Trace: Delta vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[2], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' legend("topright",
#'   legend = c("MLE", "True"),
#'   col = c("#8B0000", "#006400"),
#'   pch = c(19, 17),
#'   bty = "n", cex = 0.8
#' )
#'
#'
#' ## Example 6: Confidence Ellipses (All pairs side by side)
#'
#' # Extract all 2x2 submatrices
#' vcov_gd <- vcov_matrix[1:2, 1:2]
#' vcov_gl <- vcov_matrix[c(1, 3), c(1, 3)]
#' vcov_dl <- vcov_matrix[2:3, 2:3]
#'
#' # Create confidence ellipses
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' # Gamma vs Delta ellipse
#' eig_decomp_gd <- eigen(vcov_gd)
#' ellipse_gd <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_gd[i, ] <- mle[1:2] + sqrt(chi2_val) *
#'     (eig_decomp_gd$vectors %*% diag(sqrt(eig_decomp_gd$values)) %*% v)
#' }
#'
#' # Gamma vs Lambda ellipse
#' eig_decomp_gl <- eigen(vcov_gl)
#' ellipse_gl <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_gl[i, ] <- mle[c(1, 3)] + sqrt(chi2_val) *
#'     (eig_decomp_gl$vectors %*% diag(sqrt(eig_decomp_gl$values)) %*% v)
#' }
#'
#' # Delta vs Lambda ellipse
#' eig_decomp_dl <- eigen(vcov_dl)
#' ellipse_dl <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_dl[i, ] <- mle[2:3] + sqrt(chi2_val) *
#'     (eig_decomp_dl$vectors %*% diag(sqrt(eig_decomp_dl$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_gd <- sqrt(diag(vcov_gd))
#' ci_gamma_gd <- mle[1] + c(-1, 1) * 1.96 * se_gd[1]
#' ci_delta_gd <- mle[2] + c(-1, 1) * 1.96 * se_gd[2]
#'
#' se_gl <- sqrt(diag(vcov_gl))
#' ci_gamma_gl <- mle[1] + c(-1, 1) * 1.96 * se_gl[1]
#' ci_lambda_gl <- mle[3] + c(-1, 1) * 1.96 * se_gl[2]
#'
#' se_dl <- sqrt(diag(vcov_dl))
#' ci_delta_dl <- mle[2] + c(-1, 1) * 1.96 * se_dl[1]
#' ci_lambda_dl <- mle[3] + c(-1, 1) * 1.96 * se_dl[2]
#'
#' # Plot
#'
#' # Gamma vs Delta
#' plot(ellipse_gd[, 1], ellipse_gd[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Gamma vs Delta", las = 1, xlim = range(ellipse_gd[, 1], ci_gamma_gd),
#'   ylim = range(ellipse_gd[, 2], ci_delta_gd)
#' )
#' abline(v = ci_gamma_gd, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_delta_gd, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Gamma vs Lambda
#' plot(ellipse_gl[, 1], ellipse_gl[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = expression(lambda),
#'   main = "Gamma vs Lambda", las = 1, xlim = range(ellipse_gl[, 1], ci_gamma_gl),
#'   ylim = range(ellipse_gl[, 2], ci_lambda_gl)
#' )
#' abline(v = ci_gamma_gl, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_lambda_gl, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Delta vs Lambda
#' plot(ellipse_dl[, 1], ellipse_dl[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "Delta vs Lambda", las = 1, xlim = range(ellipse_dl[, 1], ci_delta_dl),
#'   ylim = range(ellipse_dl[, 2], ci_lambda_dl)
#' )
#' abline(v = ci_delta_dl, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_lambda_dl, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[2], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' legend("topright",
#'   legend = c("MLE", "True", "95% CR", "Marginal 95% CI"),
#'   col = c("#8B0000", "#006400", "#2E4057", "#808080"),
#'   pch = c(19, 17, NA, NA),
#'   lty = c(NA, NA, 1, 3),
#'   lwd = c(NA, NA, 2, 1.5),
#'   bty = "n", cex = 0.8
#' )
#' }
#'
#' @export
grmc <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 3) {
    stop("'par' must be a numeric vector of length 3")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_grmc",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 7. HESSIAN (hsmc)
# ----------------------------------------------------------------------------#

#' @title Hessian Matrix of the Negative Log-Likelihood for the McDonald (Mc)/Beta Power Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize hessian mcdonald
#'
#' @description
#' Computes the analytic 3x3 Hessian matrix (matrix of second partial derivatives)
#' of the negative log-likelihood function for the McDonald (Mc) distribution
#' (also known as Beta Power) with parameters \code{gamma} (\eqn{\gamma}),
#' \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}). This distribution
#' is the special case of the Generalized Kumaraswamy (GKw) distribution where
#' \eqn{\alpha = 1} and \eqn{\beta = 1}. The Hessian is useful for estimating
#' standard errors and in optimization algorithms.
#'
#' @param par A numeric vector of length 3 containing the distribution parameters
#'   in the order: \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}),
#'   \code{lambda} (\eqn{\lambda > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a 3x3 numeric matrix representing the Hessian matrix of the
#'   negative log-likelihood function, \eqn{-\partial^2 \ell / (\partial \theta_i \partial \theta_j)},
#'   where \eqn{\theta = (\gamma, \delta, \lambda)}.
#'   Returns a 3x3 matrix populated with \code{NaN} if any parameter values are
#'   invalid according to their constraints, or if any value in \code{data} is
#'   not in the interval (0, 1).
#'
#' @details
#' This function calculates the analytic second partial derivatives of the
#' negative log-likelihood function (\eqn{-\ell(\theta|\mathbf{x})}).
#' The components are based on the second derivatives of the log-likelihood \eqn{\ell}
#' (derived from the PDF in \code{\link{dmc}}).
#'
#' **Note:** The formulas below represent the second derivatives of the positive
#' log-likelihood (\eqn{\ell}). The function returns the **negative** of these values.
#' Users should verify these formulas independently if using for critical applications.
#'
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \gamma^2} = -n[\psi'(\gamma) - \psi'(\gamma+\delta+1)]
#' }
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \gamma \partial \delta} = -n\psi'(\gamma+\delta+1)
#' }
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \gamma \partial \lambda} = \sum_{i=1}^{n}\ln(x_i)
#' }
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \delta^2} = -n[\psi'(\delta+1) - \psi'(\gamma+\delta+1)]
#' }
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \delta \partial \lambda} = -\sum_{i=1}^{n}\frac{x_i^{\lambda}\ln(x_i)}{1-x_i^{\lambda}}
#' }
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \lambda^2} = -\frac{n}{\lambda^2} -
#' \delta\sum_{i=1}^{n}\frac{x_i^{\lambda}[\ln(x_i)]^2}{(1-x_i^{\lambda})^2}
#' }
#'
#' where \eqn{\psi'(\cdot)} is the trigamma function (\code{\link[base]{trigamma}}).
#' (*Note: The formula for \eqn{\partial^2 \ell / \partial \lambda^2} provided in the source
#' comment was different and potentially related to the expected information matrix;
#' the formula shown here is derived from the gradient provided earlier. Verification
#' is recommended.*)
#'
#' The returned matrix is symmetric, with rows/columns corresponding to
#' \eqn{\gamma, \delta, \lambda}.
#'
#' @references
#' McDonald, J. B. (1984). Some generalized functions for the size distribution
#' of income. *Econometrica*, 52(3), 647-663.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' (Note: Specific Hessian formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{hsgkw}} (parent distribution Hessian),
#' \code{\link{llmc}} (negative log-likelihood for Mc),
#' \code{\link{grmc}} (gradient for Mc, if available),
#' \code{\link{dmc}} (density for Mc),
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{hessian}} (for numerical Hessian comparison),
#' \code{\link[base]{trigamma}}.
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Hessian Evaluation
#'
#' # Generate sample data with more stable parameters
#' set.seed(123)
#' n <- 1000
#' true_params <- c(gamma = 2.0, delta = 2.5, lambda = 1.5)
#' data <- rmc(n,
#'   gamma = true_params[1], delta = true_params[2],
#'   lambda = true_params[3]
#' )
#'
#' # Evaluate Hessian at true parameters
#' hess_true <- hsmc(par = true_params, data = data)
#' cat("Hessian matrix at true parameters:\n")
#' print(hess_true, digits = 4)
#'
#' # Check symmetry
#' cat(
#'   "\nSymmetry check (max |H - H^T|):",
#'   max(abs(hess_true - t(hess_true))), "\n"
#' )
#'
#'
#' ## Example 2: Hessian Properties at MLE
#'
#' # Fit model
#' fit <- optim(
#'   par = c(1.5, 2.0, 1.0),
#'   fn = llmc,
#'   gr = grmc,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("gamma", "delta", "lambda")
#'
#' # Hessian at MLE
#' hessian_at_mle <- hsmc(par = mle, data = data)
#' cat("\nHessian at MLE:\n")
#' print(hessian_at_mle, digits = 4)
#'
#' # Compare with optim's numerical Hessian
#' cat("\nComparison with optim Hessian:\n")
#' cat(
#'   "Max absolute difference:",
#'   max(abs(hessian_at_mle - fit$hessian)), "\n"
#' )
#'
#' # Eigenvalue analysis
#' eigenvals <- eigen(hessian_at_mle, only.values = TRUE)$values
#' cat("\nEigenvalues:\n")
#' print(eigenvals)
#'
#' cat("\nPositive definite:", all(eigenvals > 0), "\n")
#' cat("Condition number:", max(eigenvals) / min(eigenvals), "\n")
#'
#'
#' ## Example 3: Standard Errors and Confidence Intervals
#'
#' # Observed information matrix
#' obs_info <- hessian_at_mle
#'
#' # Variance-covariance matrix
#' vcov_matrix <- solve(obs_info)
#' cat("\nVariance-Covariance Matrix:\n")
#' print(vcov_matrix, digits = 6)
#'
#' # Standard errors
#' se <- sqrt(diag(vcov_matrix))
#' names(se) <- c("gamma", "delta", "lambda")
#'
#' # Correlation matrix
#' corr_matrix <- cov2cor(vcov_matrix)
#' cat("\nCorrelation Matrix:\n")
#' print(corr_matrix, digits = 4)
#'
#' # Confidence intervals
#' z_crit <- qnorm(0.975)
#' results <- data.frame(
#'   Parameter = c("gamma", "delta", "lambda"),
#'   True = true_params,
#'   MLE = mle,
#'   SE = se,
#'   CI_Lower = mle - z_crit * se,
#'   CI_Upper = mle + z_crit * se
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 4: Determinant and Trace Analysis
#'
#' # Compute at different points
#' test_params <- rbind(
#'   c(1.5, 2.0, 1.0),
#'   c(2.0, 2.5, 1.5),
#'   mle,
#'   c(2.5, 3.0, 2.0)
#' )
#'
#' hess_properties <- data.frame(
#'   Gamma = numeric(),
#'   Delta = numeric(),
#'   Lambda = numeric(),
#'   Determinant = numeric(),
#'   Trace = numeric(),
#'   Min_Eigenval = numeric(),
#'   Max_Eigenval = numeric(),
#'   Cond_Number = numeric(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (i in 1:nrow(test_params)) {
#'   H <- hsmc(par = test_params[i, ], data = data)
#'   eigs <- eigen(H, only.values = TRUE)$values
#'
#'   hess_properties <- rbind(hess_properties, data.frame(
#'     Gamma = test_params[i, 1],
#'     Delta = test_params[i, 2],
#'     Lambda = test_params[i, 3],
#'     Determinant = det(H),
#'     Trace = sum(diag(H)),
#'     Min_Eigenval = min(eigs),
#'     Max_Eigenval = max(eigs),
#'     Cond_Number = max(eigs) / min(eigs)
#'   ))
#' }
#'
#' cat("\nHessian Properties at Different Points:\n")
#' print(hess_properties, digits = 4, row.names = FALSE)
#'
#'
#' ## Example 5: Curvature Visualization (All pairs side by side)
#'
#' # Create grids around MLE with wider range (±1.5)
#' gamma_grid <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = 25)
#' delta_grid <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = 25)
#' lambda_grid <- seq(mle[3] - 1.5, mle[3] + 1.5, length.out = 25)
#'
#' gamma_grid <- gamma_grid[gamma_grid > 0]
#' delta_grid <- delta_grid[delta_grid > 0]
#' lambda_grid <- lambda_grid[lambda_grid > 0]
#'
#' # Compute curvature measures for all pairs
#' determinant_surface_gd <- matrix(NA, nrow = length(gamma_grid), ncol = length(delta_grid))
#' trace_surface_gd <- matrix(NA, nrow = length(gamma_grid), ncol = length(delta_grid))
#'
#' determinant_surface_gl <- matrix(NA, nrow = length(gamma_grid), ncol = length(lambda_grid))
#' trace_surface_gl <- matrix(NA, nrow = length(gamma_grid), ncol = length(lambda_grid))
#'
#' determinant_surface_dl <- matrix(NA, nrow = length(delta_grid), ncol = length(lambda_grid))
#' trace_surface_dl <- matrix(NA, nrow = length(delta_grid), ncol = length(lambda_grid))
#'
#' # Gamma vs Delta
#' for (i in seq_along(gamma_grid)) {
#'   for (j in seq_along(delta_grid)) {
#'     H <- hsmc(c(gamma_grid[i], delta_grid[j], mle[3]), data)
#'     determinant_surface_gd[i, j] <- det(H)
#'     trace_surface_gd[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Gamma vs Lambda
#' for (i in seq_along(gamma_grid)) {
#'   for (j in seq_along(lambda_grid)) {
#'     H <- hsmc(c(gamma_grid[i], mle[2], lambda_grid[j]), data)
#'     determinant_surface_gl[i, j] <- det(H)
#'     trace_surface_gl[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Delta vs Lambda
#' for (i in seq_along(delta_grid)) {
#'   for (j in seq_along(lambda_grid)) {
#'     H <- hsmc(c(mle[1], delta_grid[i], lambda_grid[j]), data)
#'     determinant_surface_dl[i, j] <- det(H)
#'     trace_surface_dl[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Plot
#'
#' # Determinant plots
#' contour(gamma_grid, delta_grid, determinant_surface_gd,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Determinant: Gamma vs Delta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(gamma_grid, lambda_grid, determinant_surface_gl,
#'   xlab = expression(gamma), ylab = expression(lambda),
#'   main = "Determinant: Gamma vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(delta_grid, lambda_grid, determinant_surface_dl,
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "Determinant: Delta vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[2], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Trace plots
#' contour(gamma_grid, delta_grid, trace_surface_gd,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Trace: Gamma vs Delta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(gamma_grid, lambda_grid, trace_surface_gl,
#'   xlab = expression(gamma), ylab = expression(lambda),
#'   main = "Trace: Gamma vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(delta_grid, lambda_grid, trace_surface_dl,
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "Trace: Delta vs Lambda", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[2], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' legend("topright",
#'   legend = c("MLE", "True"),
#'   col = c("#8B0000", "#006400"),
#'   pch = c(19, 17),
#'   bty = "n", cex = 0.8
#' )
#'
#' ## Example 6: Confidence Ellipses (All pairs side by side)
#'
#' # Extract all 2x2 submatrices
#' vcov_gd <- vcov_matrix[1:2, 1:2]
#' vcov_gl <- vcov_matrix[c(1, 3), c(1, 3)]
#' vcov_dl <- vcov_matrix[2:3, 2:3]
#'
#' # Create confidence ellipses
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' # Gamma vs Delta ellipse
#' eig_decomp_gd <- eigen(vcov_gd)
#' ellipse_gd <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_gd[i, ] <- mle[1:2] + sqrt(chi2_val) *
#'     (eig_decomp_gd$vectors %*% diag(sqrt(eig_decomp_gd$values)) %*% v)
#' }
#'
#' # Gamma vs Lambda ellipse
#' eig_decomp_gl <- eigen(vcov_gl)
#' ellipse_gl <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_gl[i, ] <- mle[c(1, 3)] + sqrt(chi2_val) *
#'     (eig_decomp_gl$vectors %*% diag(sqrt(eig_decomp_gl$values)) %*% v)
#' }
#'
#' # Delta vs Lambda ellipse
#' eig_decomp_dl <- eigen(vcov_dl)
#' ellipse_dl <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_dl[i, ] <- mle[2:3] + sqrt(chi2_val) *
#'     (eig_decomp_dl$vectors %*% diag(sqrt(eig_decomp_dl$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_gd <- sqrt(diag(vcov_gd))
#' ci_gamma_gd <- mle[1] + c(-1, 1) * 1.96 * se_gd[1]
#' ci_delta_gd <- mle[2] + c(-1, 1) * 1.96 * se_gd[2]
#'
#' se_gl <- sqrt(diag(vcov_gl))
#' ci_gamma_gl <- mle[1] + c(-1, 1) * 1.96 * se_gl[1]
#' ci_lambda_gl <- mle[3] + c(-1, 1) * 1.96 * se_gl[2]
#'
#' se_dl <- sqrt(diag(vcov_dl))
#' ci_delta_dl <- mle[2] + c(-1, 1) * 1.96 * se_dl[1]
#' ci_lambda_dl <- mle[3] + c(-1, 1) * 1.96 * se_dl[2]
#'
#' # Plot
#'
#' # Gamma vs Delta
#' plot(ellipse_gd[, 1], ellipse_gd[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Gamma vs Delta", las = 1, xlim = range(ellipse_gd[, 1], ci_gamma_gd),
#'   ylim = range(ellipse_gd[, 2], ci_delta_gd)
#' )
#' abline(v = ci_gamma_gd, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_delta_gd, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Gamma vs Lambda
#' plot(ellipse_gl[, 1], ellipse_gl[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = expression(lambda),
#'   main = "Gamma vs Lambda", las = 1, xlim = range(ellipse_gl[, 1], ci_gamma_gl),
#'   ylim = range(ellipse_gl[, 2], ci_lambda_gl)
#' )
#' abline(v = ci_gamma_gl, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_lambda_gl, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Delta vs Lambda
#' plot(ellipse_dl[, 1], ellipse_dl[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "Delta vs Lambda", las = 1, xlim = range(ellipse_dl[, 1], ci_delta_dl),
#'   ylim = range(ellipse_dl[, 2], ci_lambda_dl)
#' )
#' abline(v = ci_delta_dl, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_lambda_dl, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[2], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' legend("topright",
#'   legend = c("MLE", "True", "95% CR", "Marginal 95% CI"),
#'   col = c("#8B0000", "#006400", "#2E4057", "#808080"),
#'   pch = c(19, 17, NA, NA),
#'   lty = c(NA, NA, 1, 3),
#'   lwd = c(NA, NA, 2, 1.5),
#'   bty = "n", cex = 0.8
#' )
#' }
#'
#' @export
hsmc <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 3) {
    stop("'par' must be a numeric vector of length 3")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_hsmc",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}
