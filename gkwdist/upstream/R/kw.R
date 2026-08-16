# ============================================================================#
# KUMARASWAMY (Kw) DISTRIBUTION
# ============================================================================#
#
# Wrapper functions for the two-parameter Kumaraswamy distribution.
# C++ implementations are in src/kw.cpp
#
# Functions:
#   - dkw: Probability density function (PDF)
#   - pkw: Cumulative distribution function (CDF)
#   - qkw: Quantile function (inverse CDF)
#   - rkw: Random number generation
#   - llkw: Negative log-likelihood
#   - grkw: Gradient of negative log-likelihood
#   - hskw: Hessian of negative log-likelihood
# ============================================================================#


# ----------------------------------------------------------------------------#
# 1. DENSITY FUNCTION (dkw)
# ----------------------------------------------------------------------------#


#' @title Density of the Kumaraswamy (Kw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution density kumaraswamy
#'
#' @description
#' Computes the probability density function (PDF) for the two-parameter
#' Kumaraswamy (Kw) distribution with shape parameters \code{alpha} (\eqn{\alpha})
#' and \code{beta} (\eqn{\beta}). This distribution is defined on the interval (0, 1).
#'
#' @param x Vector of quantiles (values between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param log Logical; if \code{TRUE}, the logarithm of the density is
#'   returned (\eqn{\log(f(x))}). Default: \code{FALSE}.
#'
#' @return A vector of density values (\eqn{f(x)}) or log-density values
#'   (\eqn{\log(f(x))}). The length of the result is determined by the recycling
#'   rule applied to the arguments (\code{x}, \code{alpha}, \code{beta}).
#'   Returns \code{0} (or \code{-Inf} if \code{log = TRUE}) for \code{x}
#'   outside the interval (0, 1), or \code{NaN} if parameters are invalid
#'   (e.g., \code{alpha <= 0}, \code{beta <= 0}).
#'
#' @details
#' The probability density function (PDF) of the Kumaraswamy (Kw) distribution
#' is given by:
#' \deqn{
#' f(x; \alpha, \beta) = \alpha \beta x^{\alpha-1} (1 - x^\alpha)^{\beta-1}
#' }
#' for \eqn{0 < x < 1}, \eqn{\alpha > 0}, and \eqn{\beta > 0}.
#'
#' The Kumaraswamy distribution is identical to the Generalized Kumaraswamy (GKw)
#' distribution (\code{\link{dgkw}}) with parameters \eqn{\gamma = 1},
#' \eqn{\delta = 0}, and \eqn{\lambda = 1}. It is also a special case of the
#' Exponentiated Kumaraswamy (\code{\link{dekw}}) with \eqn{\lambda = 1}, and
#' the Kumaraswamy-Kumaraswamy (\code{\link{dkkw}}) with \eqn{\delta = 0}
#' and \eqn{\lambda = 1}.
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
#' with some tractability advantages. *Statistical Methodology*, *6*(1), 70-81.
#'
#'
#' @seealso
#' \code{\link{dgkw}} (parent distribution density),
#' \code{\link{dekw}}, \code{\link{dkkw}},
#' \code{\link{pkw}}, \code{\link{qkw}}, \code{\link{rkw}} (other Kw functions),
#' \code{\link[stats]{dbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' x_vals <- c(0.2, 0.5, 0.8)
#' alpha_par <- 2.0
#' beta_par <- 3.0
#'
#' # Calculate density using dkw
#' densities <- dkw(x_vals, alpha_par, beta_par)
#' print(densities)
#'
#' # Calculate log-density
#' log_densities <- dkw(x_vals, alpha_par, beta_par, log = TRUE)
#' print(log_densities)
#' # Check: should match log(densities)
#' print(log(densities))
#'
#' # Compare with dgkw setting gamma = 1, delta = 0, lambda = 1
#' densities_gkw <- dgkw(x_vals, alpha_par, beta_par,
#'   gamma = 1.0, delta = 0.0,
#'   lambda = 1.0
#' )
#' print(paste("Max difference:", max(abs(densities - densities_gkw)))) # Should be near zero
#'
#' # Plot the density for different shape parameter combinations
#' curve_x <- seq(0.001, 0.999, length.out = 200)
#' plot(curve_x, dkw(curve_x, alpha = 2, beta = 3),
#'   type = "l",
#'   main = "Kumaraswamy Density Examples", xlab = "x", ylab = "f(x)",
#'   col = "blue", ylim = c(0, 4)
#' )
#' lines(curve_x, dkw(curve_x, alpha = 3, beta = 2), col = "red")
#' lines(curve_x, dkw(curve_x, alpha = 0.5, beta = 0.5), col = "green") # U-shaped
#' lines(curve_x, dkw(curve_x, alpha = 5, beta = 1), col = "purple") # J-shaped
#' lines(curve_x, dkw(curve_x, alpha = 1, beta = 3), col = "orange") # J-shaped (reversed)
#' legend("top",
#'   legend = c("a=2, b=3", "a=3, b=2", "a=0.5, b=0.5", "a=5, b=1", "a=1, b=3"),
#'   col = c("blue", "red", "green", "purple", "orange"), lty = 1, bty = "n", ncol = 2
#' )
#' }
#' @export
dkw <- function(x, alpha = 1, beta = 1, log = FALSE) {
  # Input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("'log' must be a single logical value")
  }

  # Call C++ implementation
  .Call("_gkwdist_dkw",
    as.numeric(x),
    as.numeric(alpha),
    as.numeric(beta),
    as.logical(log),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 2. DISTRIBUTION FUNCTION (pkw)
# ----------------------------------------------------------------------------#

#' @title Cumulative Distribution Function (CDF) of the Kumaraswamy (Kw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution cumulative kumaraswamy
#'
#' @description
#' Computes the cumulative distribution function (CDF), \eqn{P(X \le q)}, for the
#' two-parameter Kumaraswamy (Kw) distribution with shape parameters \code{alpha}
#' (\eqn{\alpha}) and \code{beta} (\eqn{\beta}). This distribution is defined
#' on the interval (0, 1).
#'
#' @param q Vector of quantiles (values generally between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}, otherwise, \eqn{P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of probabilities, \eqn{F(q)}, or their logarithms/complements
#'   depending on \code{lower.tail} and \code{log.p}. The length of the result
#'   is determined by the recycling rule applied to the arguments (\code{q},
#'   \code{alpha}, \code{beta}). Returns \code{0} (or \code{-Inf} if
#'   \code{log.p = TRUE}) for \code{q <= 0} and \code{1} (or \code{0} if
#'   \code{log.p = TRUE}) for \code{q >= 1}. Returns \code{NaN} for invalid
#'   parameters.
#'
#' @details
#' The cumulative distribution function (CDF) of the Kumaraswamy (Kw)
#' distribution is given by:
#' \deqn{
#' F(x; \alpha, \beta) = 1 - (1 - x^\alpha)^\beta
#' }
#' for \eqn{0 < x < 1}, \eqn{\alpha > 0}, and \eqn{\beta > 0}.
#'
#' The Kw distribution is a special case of several generalized distributions:
#' \itemize{
#'  \item Generalized Kumaraswamy (\code{\link{pgkw}}) with \eqn{\gamma=1, \delta=0, \lambda=1}.
#'  \item Exponentiated Kumaraswamy (\code{\link{pekw}}) with \eqn{\lambda=1}.
#'  \item Kumaraswamy-Kumaraswamy (\code{\link{pkkw}}) with \eqn{\delta=0, \lambda=1}.
#' }
#' The implementation uses the closed-form expression for efficiency.
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
#' with some tractability advantages. *Statistical Methodology*, *6*(1), 70-81.
#'
#'
#' @seealso
#' \code{\link{pgkw}}, \code{\link{pekw}}, \code{\link{pkkw}} (related generalized CDFs),
#' \code{\link{dkw}}, \code{\link{qkw}}, \code{\link{rkw}} (other Kw functions),
#' \code{\link[stats]{pbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' q_vals <- c(0.2, 0.5, 0.8)
#' alpha_par <- 2.0
#' beta_par <- 3.0
#'
#' # Calculate CDF P(X <= q) using pkw
#' probs <- pkw(q_vals, alpha_par, beta_par)
#' print(probs)
#'
#' # Calculate upper tail P(X > q)
#' probs_upper <- pkw(q_vals, alpha_par, beta_par, lower.tail = FALSE)
#' print(probs_upper)
#' # Check: probs + probs_upper should be 1
#' print(probs + probs_upper)
#'
#' # Calculate log CDF
#' logs <- pkw(q_vals, alpha_par, beta_par, log.p = TRUE)
#' print(logs)
#' # Check: should match log(probs)
#' print(log(probs))
#'
#' # Compare with pgkw setting gamma = 1, delta = 0, lambda = 1
#' probs_gkw <- pgkw(q_vals, alpha_par, beta_par,
#'   gamma = 1.0, delta = 0.0,
#'   lambda = 1.0
#' )
#' print(paste("Max difference:", max(abs(probs - probs_gkw)))) # Should be near zero
#'
#' # Plot the CDF for different shape parameter combinations
#' curve_q <- seq(0.001, 0.999, length.out = 200)
#' plot(curve_q, pkw(curve_q, alpha = 2, beta = 3),
#'   type = "l",
#'   main = "Kumaraswamy CDF Examples", xlab = "q", ylab = "F(q)",
#'   col = "blue", ylim = c(0, 1)
#' )
#' lines(curve_q, pkw(curve_q, alpha = 3, beta = 2), col = "red")
#' lines(curve_q, pkw(curve_q, alpha = 0.5, beta = 0.5), col = "green")
#' lines(curve_q, pkw(curve_q, alpha = 5, beta = 1), col = "purple")
#' lines(curve_q, pkw(curve_q, alpha = 1, beta = 3), col = "orange")
#' legend("bottomright",
#'   legend = c("a=2, b=3", "a=3, b=2", "a=0.5, b=0.5", "a=5, b=1", "a=1, b=3"),
#'   col = c("blue", "red", "green", "purple", "orange"), lty = 1, bty = "n", ncol = 2
#' )
#' }
#'
#' @export
pkw <- function(q, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(q)) stop("'q' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
  }
  if (!is.logical(lower.tail) || length(lower.tail) != 1) {
    stop("'lower.tail' must be a single logical value")
  }
  if (!is.logical(log.p) || length(log.p) != 1) {
    stop("'log.p' must be a single logical value")
  }

  # Call C++ implementation
  .Call("_gkwdist_pkw",
    as.numeric(q),
    as.numeric(alpha),
    as.numeric(beta),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 3. QUANTILE FUNCTION (qkw)
# ----------------------------------------------------------------------------#

#' @title Quantile Function of the Kumaraswamy (Kw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution quantile kumaraswamy
#'
#' @description
#' Computes the quantile function (inverse CDF) for the two-parameter
#' Kumaraswamy (Kw) distribution with shape parameters \code{alpha} (\eqn{\alpha})
#' and \code{beta} (\eqn{\beta}). It finds the value \code{q} such that
#' \eqn{P(X \le q) = p}.
#'
#' @param p Vector of probabilities (values between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{p = P(X \le q)},
#'   otherwise, probabilities are \eqn{p = P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of quantiles corresponding to the given probabilities \code{p}.
#'   The length of the result is determined by the recycling rule applied to
#'   the arguments (\code{p}, \code{alpha}, \code{beta}).
#'   Returns:
#'   \itemize{
#'     \item \code{0} for \code{p = 0} (or \code{p = -Inf} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{1} for \code{p = 1} (or \code{p = 0} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{NaN} for \code{p < 0} or \code{p > 1} (or corresponding log scale).
#'     \item \code{NaN} for invalid parameters (e.g., \code{alpha <= 0},
#'           \code{beta <= 0}).
#'   }
#'   Boundary return values are adjusted accordingly for \code{lower.tail = FALSE}.
#'
#' @details
#' The quantile function \eqn{Q(p)} is the inverse of the CDF \eqn{F(q)}. The CDF
#' for the Kumaraswamy distribution is \eqn{F(q) = 1 - (1 - q^\alpha)^\beta}
#' (see \code{\link{pkw}}). Inverting this equation for \eqn{q} yields the
#' quantile function:
#' \deqn{
#' Q(p) = \left\{ 1 - (1 - p)^{1/\beta} \right\}^{1/\alpha}
#' }
#' The function uses this closed-form expression and correctly handles the
#' \code{lower.tail} and \code{log.p} arguments by transforming \code{p}
#' appropriately before applying the formula. This is equivalent to the general
#' GKw quantile function (\code{\link{qgkw}}) evaluated with \eqn{\gamma=1, \delta=0, \lambda=1}.
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
#' with some tractability advantages. *Statistical Methodology*, *6*(1), 70-81.
#'
#'
#' @seealso
#' \code{\link{qgkw}} (parent distribution quantile function),
#' \code{\link{dkw}}, \code{\link{pkw}}, \code{\link{rkw}} (other Kw functions),
#' \code{\link[stats]{qbeta}}, \code{\link[stats]{qunif}}
#'
#' @examples
#' \donttest{
#' # Example values
#' p_vals <- c(0.1, 0.5, 0.9)
#' alpha_par <- 2.0
#' beta_par <- 3.0
#'
#' # Calculate quantiles using qkw
#' quantiles <- qkw(p_vals, alpha_par, beta_par)
#' print(quantiles)
#'
#' # Calculate quantiles for upper tail probabilities P(X > q) = p
#' quantiles_upper <- qkw(p_vals, alpha_par, beta_par, lower.tail = FALSE)
#' print(quantiles_upper)
#'
#' # Calculate quantiles from log probabilities
#' log.p_vals <- log(p_vals)
#' quantiles_logp <- qkw(log.p_vals, alpha_par, beta_par, log.p = TRUE)
#' print(quantiles_logp)
#' # Check: should match original quantiles
#' print(quantiles)
#'
#' # Compare with qgkw setting gamma = 1, delta = 0, lambda = 1
#' quantiles_gkw <- qgkw(p_vals,
#'   alpha = alpha_par, beta = beta_par,
#'   gamma = 1.0, delta = 0.0, lambda = 1.0
#' )
#' print(paste("Max difference:", max(abs(quantiles - quantiles_gkw)))) # Should be near zero
#'
#' # Verify inverse relationship with pkw
#' p_check <- 0.75
#' q_calc <- qkw(p_check, alpha_par, beta_par)
#' p_recalc <- pkw(q_calc, alpha_par, beta_par)
#' print(paste("Original p:", p_check, " Recalculated p:", p_recalc))
#' # abs(p_check - p_recalc) < 1e-9 # Should be TRUE
#'
#' # Boundary conditions
#' print(qkw(c(0, 1), alpha_par, beta_par)) # Should be 0, 1
#' print(qkw(c(-Inf, 0), alpha_par, beta_par, log.p = TRUE)) # Should be 0, 1
#' }
#'
#' @export
qkw <- function(p, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(p)) stop("'p' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
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
  .Call("_gkwdist_qkw",
    as.numeric(p),
    as.numeric(alpha),
    as.numeric(beta),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 4. RANDOM GENERATION (rkw)
# ----------------------------------------------------------------------------#

#' @title Random Number Generation for the Kumaraswamy (Kw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution random kumaraswamy
#'
#' @description
#' Generates random deviates from the two-parameter Kumaraswamy (Kw)
#' distribution with shape parameters \code{alpha} (\eqn{\alpha}) and
#' \code{beta} (\eqn{\beta}).
#'
#' @param n Number of observations. If \code{length(n) > 1}, the length is
#'   taken to be the number required. Must be a non-negative integer.
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#'
#' @return A vector of length \code{n} containing random deviates from the Kw
#'   distribution, with values in (0, 1). The length of the result is determined
#'   by \code{n} and the recycling rule applied to the parameters (\code{alpha},
#'   \code{beta}). Returns \code{NaN} if parameters are invalid (e.g.,
#'   \code{alpha <= 0}, \code{beta <= 0}).
#'
#' @details
#' The generation method uses the inverse transform (quantile) method.
#' That is, if \eqn{U} is a random variable following a standard Uniform
#' distribution on (0, 1), then \eqn{X = Q(U)} follows the Kw distribution,
#' where \eqn{Q(p)} is the Kw quantile function (\code{\link{qkw}}):
#' \deqn{
#' Q(p) = \left\{ 1 - (1 - p)^{1/\beta} \right\}^{1/\alpha}
#' }
#' The implementation generates \eqn{U} using \code{\link[stats]{runif}}
#' and applies this transformation. This is equivalent to the general GKw
#' generation method (\code{\link{rgkw}}) evaluated at \eqn{\gamma=1, \delta=0, \lambda=1}.
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
#' with some tractability advantages. *Statistical Methodology*, *6*(1), 70-81.
#'
#'
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer-Verlag.
#' (General methods for random variate generation).
#'
#' @seealso
#' \code{\link{rgkw}} (parent distribution random generation),
#' \code{\link{dkw}}, \code{\link{pkw}}, \code{\link{qkw}} (other Kw functions),
#' \code{\link[stats]{runif}}
#'
#' @examples
#' \donttest{
#' set.seed(2029) # for reproducibility
#'
#' # Generate 1000 random values from a specific Kw distribution
#' alpha_par <- 2.0
#' beta_par <- 3.0
#'
#' x_sample_kw <- rkw(1000, alpha = alpha_par, beta = beta_par)
#' summary(x_sample_kw)
#'
#' # Histogram of generated values compared to theoretical density
#' hist(x_sample_kw,
#'   breaks = 30, freq = FALSE, # freq=FALSE for density
#'   main = "Histogram of Kw Sample", xlab = "x", ylim = c(0, 2.5)
#' )
#' curve(dkw(x, alpha = alpha_par, beta = beta_par),
#'   add = TRUE, col = "red", lwd = 2, n = 201
#' )
#' legend("top", legend = "Theoretical PDF", col = "red", lwd = 2, bty = "n")
#'
#' # Comparing empirical and theoretical quantiles (Q-Q plot)
#' prob_points <- seq(0.01, 0.99, by = 0.01)
#' theo_quantiles <- qkw(prob_points, alpha = alpha_par, beta = beta_par)
#' emp_quantiles <- quantile(x_sample_kw, prob_points, type = 7)
#'
#' plot(theo_quantiles, emp_quantiles,
#'   pch = 16, cex = 0.8,
#'   main = "Q-Q Plot for Kw Distribution",
#'   xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles (n=1000)"
#' )
#' abline(a = 0, b = 1, col = "blue", lty = 2)
#'
#' # Compare summary stats with rgkw(..., gamma=1, delta=0, lambda=1)
#' # Note: individual values will differ due to randomness
#' x_sample_gkw <- rgkw(1000,
#'   alpha = alpha_par, beta = beta_par, gamma = 1.0,
#'   delta = 0.0, lambda = 1.0
#' )
#' print("Summary stats for rkw sample:")
#' print(summary(x_sample_kw))
#' print("Summary stats for rgkw(gamma=1, delta=0, lambda=1) sample:")
#' print(summary(x_sample_gkw)) # Should be similar
#' }
#'
#' @export
rkw <- function(n, alpha = 1, beta = 1) {
  # Input validation
  if (length(n) > 1) n <- length(n)
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("'n' must be a positive integer")
  }
  n <- as.integer(n)

  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
  }

  # Call C++ implementation
  .Call("_gkwdist_rkw",
    as.integer(n),
    as.numeric(alpha),
    as.numeric(beta),
    PACKAGE = "gkwdist"
  )
}


# ============================================================================#
# MAXIMUM LIKELIHOOD ESTIMATION FUNCTIONS
# ============================================================================#

# ----------------------------------------------------------------------------#
# 5. NEGATIVE LOG-LIKELIHOOD (llkw)
# ----------------------------------------------------------------------------#

#' @title Negative Log-Likelihood of the Kumaraswamy (Kw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize kumaraswamy
#'
#' @description
#' Computes the negative log-likelihood function for the two-parameter
#' Kumaraswamy (Kw) distribution with parameters \code{alpha} (\eqn{\alpha})
#' and \code{beta} (\eqn{\beta}), given a vector of observations. This function
#' is suitable for maximum likelihood estimation.
#'
#' @param par A numeric vector of length 2 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a single \code{double} value representing the negative
#'   log-likelihood (\eqn{-\ell(\theta|\mathbf{x})}). Returns \code{Inf}
#'   if any parameter values in \code{par} are invalid according to their
#'   constraints, or if any value in \code{data} is not in the interval (0, 1).
#'
#' @details
#' The Kumaraswamy (Kw) distribution's probability density function (PDF) is
#' (see \code{\link{dkw}}):
#' \deqn{
#' f(x | \theta) = \alpha \beta x^{\alpha-1} (1 - x^\alpha)^{\beta-1}
#' }
#' for \eqn{0 < x < 1} and \eqn{\theta = (\alpha, \beta)}.
#' The log-likelihood function \eqn{\ell(\theta | \mathbf{x})} for a sample
#' \eqn{\mathbf{x} = (x_1, \dots, x_n)} is \eqn{\sum_{i=1}^n \ln f(x_i | \theta)}:
#' \deqn{
#' \ell(\theta | \mathbf{x}) = n[\ln(\alpha) + \ln(\beta)]
#' + \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta-1)\ln(v_i)]
#' }
#' where \eqn{v_i = 1 - x_i^{\alpha}}.
#' This function computes and returns the *negative* log-likelihood, \eqn{-\ell(\theta|\mathbf{x})},
#' suitable for minimization using optimization routines like \code{\link[stats]{optim}}.
#' It is equivalent to the negative log-likelihood of the GKw distribution
#' (\code{\link{llgkw}}) evaluated at \eqn{\gamma=1, \delta=0, \lambda=1}.
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
#' with some tractability advantages. *Statistical Methodology*, *6*(1), 70-81.
#'
#'
#' @seealso
#' \code{\link{llgkw}} (parent distribution negative log-likelihood),
#' \code{\link{dkw}}, \code{\link{pkw}}, \code{\link{qkw}}, \code{\link{rkw}},
#' \code{grkw} (gradient, if available),
#' \code{hskw} (Hessian, if available),
#' \code{\link[stats]{optim}}
#'
#' @examples
#' \donttest{
#' ## Example 1: Maximum Likelihood Estimation with Analytical Gradient
#'
#' # Generate sample data
#' set.seed(123)
#' n <- 1000
#' true_params <- c(alpha = 2.5, beta = 3.5)
#' data <- rkw(n, alpha = true_params[1], beta = true_params[2])
#'
#' # Optimization using BFGS with analytical gradient
#' fit <- optim(
#'   par = c(2, 2),
#'   fn = llkw,
#'   gr = grkw,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE
#' )
#'
#' # Extract results
#' mle <- fit$par
#' names(mle) <- c("alpha", "beta")
#' se <- sqrt(diag(solve(fit$hessian)))
#' ci_lower <- mle - 1.96 * se
#' ci_upper <- mle + 1.96 * se
#'
#' # Summary table
#' results <- data.frame(
#'   Parameter = c("alpha", "beta"),
#'   True = true_params,
#'   MLE = mle,
#'   SE = se,
#'   CI_Lower = ci_lower,
#'   CI_Upper = ci_upper
#' )
#' print(results, digits = 4)
#'
#' ## Example 2: Verifying Gradient at MLE
#'
#' # At MLE, gradient should be approximately zero
#' gradient_at_mle <- grkw(par = mle, data = data)
#' print(gradient_at_mle)
#' cat("Max absolute score:", max(abs(gradient_at_mle)), "\n")
#'
#' ## Example 3: Checking Hessian Properties
#'
#' # Hessian at MLE
#' hessian_at_mle <- hskw(par = mle, data = data)
#' print(hessian_at_mle, digits = 4)
#'
#' # Check positive definiteness via eigenvalues
#' eigenvals <- eigen(hessian_at_mle, only.values = TRUE)$values
#' print(eigenvals)
#' all(eigenvals > 0)
#'
#' # Condition number
#' cond_number <- max(eigenvals) / min(eigenvals)
#' cat("Condition number:", format(cond_number, scientific = TRUE), "\n")
#'
#' ## Example 4: Comparing Optimization Methods
#'
#' methods <- c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG")
#' start_params <- c(2, 2)
#'
#' comparison <- data.frame(
#'   Method = character(),
#'   Alpha_Est = numeric(),
#'   Beta_Est = numeric(),
#'   NegLogLik = numeric(),
#'   Convergence = integer(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (method in methods) {
#'   if (method %in% c("BFGS", "CG")) {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llkw,
#'       gr = grkw,
#'       data = data,
#'       method = method
#'     )
#'   } else if (method == "L-BFGS-B") {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llkw,
#'       gr = grkw,
#'       data = data,
#'       method = method,
#'       lower = c(0.01, 0.01),
#'       upper = c(100, 100)
#'     )
#'   } else {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llkw,
#'       data = data,
#'       method = method
#'     )
#'   }
#'
#'   comparison <- rbind(comparison, data.frame(
#'     Method = method,
#'     Alpha_Est = fit_temp$par[1],
#'     Beta_Est = fit_temp$par[2],
#'     NegLogLik = fit_temp$value,
#'     Convergence = fit_temp$convergence,
#'     stringsAsFactors = FALSE
#'   ))
#' }
#'
#' print(comparison, digits = 4, row.names = FALSE)
#'
#' ## Example 5: Likelihood Ratio Test
#'
#' # Test H0: beta = 3 vs H1: beta free
#' loglik_full <- -fit$value
#'
#' # Restricted model: fix beta = 3
#' restricted_ll <- function(alpha, data, beta_fixed) {
#'   llkw(par = c(alpha, beta_fixed), data = data)
#' }
#'
#' fit_restricted <- optimize(
#'   f = restricted_ll,
#'   interval = c(0.1, 10),
#'   data = data,
#'   beta_fixed = 3,
#'   maximum = FALSE
#' )
#'
#' loglik_restricted <- -fit_restricted$objective
#' lr_stat <- 2 * (loglik_full - loglik_restricted)
#' p_value <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
#'
#' cat("LR Statistic:", round(lr_stat, 4), "\n")
#' cat("P-value:", format.pval(p_value, digits = 4), "\n")
#'
#' ## Example 6: Univariate Profile Likelihoods
#'
#' # Grid for alpha
#' alpha_grid <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = 50)
#' alpha_grid <- alpha_grid[alpha_grid > 0]
#' profile_ll_alpha <- numeric(length(alpha_grid))
#'
#' for (i in seq_along(alpha_grid)) {
#'   profile_fit <- optimize(
#'     f = function(beta) llkw(c(alpha_grid[i], beta), data),
#'     interval = c(0.1, 10),
#'     maximum = FALSE
#'   )
#'   profile_ll_alpha[i] <- -profile_fit$objective
#' }
#'
#' # Grid for beta
#' beta_grid <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = 50)
#' beta_grid <- beta_grid[beta_grid > 0]
#' profile_ll_beta <- numeric(length(beta_grid))
#'
#' for (i in seq_along(beta_grid)) {
#'   profile_fit <- optimize(
#'     f = function(alpha) llkw(c(alpha, beta_grid[i]), data),
#'     interval = c(0.1, 10),
#'     maximum = FALSE
#'   )
#'   profile_ll_beta[i] <- -profile_fit$objective
#' }
#'
#' # 95% confidence threshold
#' chi_crit <- qchisq(0.95, df = 1)
#' threshold <- max(profile_ll_alpha) - chi_crit / 2
#'
#' # Plot
#'
#' # Profile for alpha
#' plot(alpha_grid, profile_ll_alpha,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile Likelihood: ", alpha)), las = 1
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
#' # Profile for beta
#' plot(beta_grid, profile_ll_beta,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(beta), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile Likelihood: ", beta)), las = 1
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
#' ## Example 7: 2D Profile Likelihood Surface
#'
#' # Create 2D grid
#' alpha_2d <- seq(mle[1] - 1, mle[1] + 1, length.out = round(n / 4))
#' beta_2d <- seq(mle[2] - 1, mle[2] + 1, length.out = round(n / 4))
#' alpha_2d <- alpha_2d[alpha_2d > 0]
#' beta_2d <- beta_2d[beta_2d > 0]
#'
#' # Compute log-likelihood surface
#' ll_surface <- matrix(NA, nrow = length(alpha_2d), ncol = length(beta_2d))
#'
#' for (i in seq_along(alpha_2d)) {
#'   for (j in seq_along(beta_2d)) {
#'     ll_surface[i, j] <- -llkw(c(alpha_2d[i], beta_2d[j]), data)
#'   }
#' }
#'
#' # Confidence region levels
#' max_ll <- max(ll_surface, na.rm = TRUE)
#' levels_90 <- max_ll - qchisq(0.90, df = 2) / 2
#' levels_95 <- max_ll - qchisq(0.95, df = 2) / 2
#' levels_99 <- max_ll - qchisq(0.99, df = 2) / 2
#'
#' # Plot contour
#' contour(alpha_2d, beta_2d, ll_surface,
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "2D Profile Log-Likelihood",
#'   levels = seq(min(ll_surface, na.rm = TRUE), max_ll, length.out = round(n / 4)),
#'   col = "#2E4057", las = 1, lwd = 1
#' )
#'
#' # Add confidence region contours
#' contour(alpha_2d, beta_2d, ll_surface,
#'   levels = c(levels_90, levels_95, levels_99),
#'   col = c("#FFA07A", "#FF6347", "#8B0000"),
#'   lwd = c(2, 2.5, 3), lty = c(3, 2, 1),
#'   add = TRUE, labcex = 0.8
#' )
#'
#' # Mark points
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#'
#' legend("topright",
#'   legend = c("MLE", "True", "90% CR", "95% CR", "99% CR"),
#'   col = c("#8B0000", "#006400", "#FFA07A", "#FF6347", "#8B0000"),
#'   pch = c(19, 17, NA, NA, NA),
#'   lty = c(NA, NA, 3, 2, 1),
#'   lwd = c(NA, NA, 2, 2.5, 3),
#'   bty = "n", cex = 0.8
#' )
#' grid(col = "gray90")
#'
#' ## Example 8: Combined View - Profiles with 2D Surface
#'
#' # Top left: Profile for alpha
#' plot(alpha_grid, profile_ll_alpha,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", alpha)), las = 1
#' )
#' abline(v = mle[1], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[1], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3)
#' grid(col = "gray90")
#'
#' # Top right: Profile for beta
#' plot(beta_grid, profile_ll_beta,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(beta), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", beta)), las = 1
#' )
#' abline(v = mle[2], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[2], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3)
#' grid(col = "gray90")
#'
#' # Bottom left: 2D contour
#' contour(alpha_2d, beta_2d, ll_surface,
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "2D Log-Likelihood Surface",
#'   levels = seq(min(ll_surface, na.rm = TRUE), max_ll, length.out = 15),
#'   col = "#2E4057", las = 1, lwd = 1
#' )
#' contour(alpha_2d, beta_2d, ll_surface,
#'   levels = c(levels_95),
#'   col = "#8B0000", lwd = 2.5, add = TRUE
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' ## Example 9: Numerical Gradient Verification
#'
#' # Manual finite difference gradient
#' numerical_gradient <- function(f, x, data, h = 1e-7) {
#'   grad <- numeric(length(x))
#'   for (i in seq_along(x)) {
#'     x_plus <- x_minus <- x
#'     x_plus[i] <- x[i] + h
#'     x_minus[i] <- x[i] - h
#'     grad[i] <- (f(x_plus, data) - f(x_minus, data)) / (2 * h)
#'   }
#'   return(grad)
#' }
#'
#' # Compare
#' grad_analytical <- grkw(par = mle, data = data)
#' grad_numerical <- numerical_gradient(llkw, mle, data)
#'
#' comparison_grad <- data.frame(
#'   Parameter = c("alpha", "beta"),
#'   Analytical = grad_analytical,
#'   Numerical = grad_numerical,
#'   Difference = abs(grad_analytical - grad_numerical)
#' )
#' print(comparison_grad, digits = 8)
#'
#' ## Example 10: Bootstrap Confidence Intervals
#'
#' n_boot <- round(n / 4)
#' boot_estimates <- matrix(NA, nrow = n_boot, ncol = 2)
#'
#' set.seed(456)
#' for (b in 1:n_boot) {
#'   boot_data <- rkw(n, alpha = mle[1], beta = mle[2])
#'   boot_fit <- optim(
#'     par = mle,
#'     fn = llkw,
#'     gr = grkw,
#'     data = boot_data,
#'     method = "BFGS",
#'     control = list(maxit = 500)
#'   )
#'   if (boot_fit$convergence == 0) {
#'     boot_estimates[b, ] <- boot_fit$par
#'   }
#' }
#'
#' boot_estimates <- boot_estimates[complete.cases(boot_estimates), ]
#' boot_ci <- apply(boot_estimates, 2, quantile, probs = c(0.025, 0.975))
#' colnames(boot_ci) <- c("alpha", "beta")
#'
#' print(t(boot_ci), digits = 4)
#'
#' # Plot bootstrap distributions
#'
#' hist(boot_estimates[, 1],
#'   breaks = 20, col = "#87CEEB", border = "white",
#'   main = expression(paste("Bootstrap: ", hat(alpha))),
#'   xlab = expression(hat(alpha)), las = 1
#' )
#' abline(v = mle[1], col = "#8B0000", lwd = 2)
#' abline(v = true_params[1], col = "#006400", lwd = 2, lty = 2)
#' abline(v = boot_ci[, 1], col = "#2E4057", lwd = 2, lty = 3)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#2E4057"),
#'   lwd = 2, lty = c(1, 2, 3), bty = "n"
#' )
#'
#' hist(boot_estimates[, 2],
#'   breaks = 20, col = "#FFA07A", border = "white",
#'   main = expression(paste("Bootstrap: ", hat(beta))),
#'   xlab = expression(hat(beta)), las = 1
#' )
#' abline(v = mle[2], col = "#8B0000", lwd = 2)
#' abline(v = true_params[2], col = "#006400", lwd = 2, lty = 2)
#' abline(v = boot_ci[, 2], col = "#2E4057", lwd = 2, lty = 3)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#2E4057"),
#'   lwd = 2, lty = c(1, 2, 3), bty = "n"
#' )
#' }
#'
#' @export
llkw <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 2) {
    stop("'par' must be a numeric vector of length 2")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_llkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 6. GRADIENT (grkw)
# ----------------------------------------------------------------------------#

#' @title Gradient of the Negative Log-Likelihood for the Kumaraswamy (Kw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize gradient kumaraswamy
#'
#' @description
#' Computes the gradient vector (vector of first partial derivatives) of the
#' negative log-likelihood function for the two-parameter Kumaraswamy (Kw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}) and \code{beta}
#' (\eqn{\beta}). This provides the analytical gradient often used for efficient
#' optimization via maximum likelihood estimation.
#'
#' @param par A numeric vector of length 2 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a numeric vector of length 2 containing the partial derivatives
#'   of the negative log-likelihood function \eqn{-\ell(\theta | \mathbf{x})} with
#'   respect to each parameter: \eqn{(-\partial \ell/\partial \alpha, -\partial \ell/\partial \beta)}.
#'   Returns a vector of \code{NaN} if any parameter values are invalid according
#'   to their constraints, or if any value in \code{data} is not in the
#'   interval (0, 1).
#'
#' @details
#' The components of the gradient vector of the negative log-likelihood
#' (\eqn{-\nabla \ell(\theta | \mathbf{x})}) for the Kw model are:
#'
#' \deqn{
#' -\frac{\partial \ell}{\partial \alpha} = -\frac{n}{\alpha} - \sum_{i=1}^{n}\ln(x_i)
#' + (\beta-1)\sum_{i=1}^{n}\frac{x_i^{\alpha}\ln(x_i)}{v_i}
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \beta} = -\frac{n}{\beta} - \sum_{i=1}^{n}\ln(v_i)
#' }
#'
#' where \eqn{v_i = 1 - x_i^{\alpha}}.
#' These formulas represent the derivatives of \eqn{-\ell(\theta)}, consistent with
#' minimizing the negative log-likelihood. They correspond to the relevant components
#' of the general GKw gradient (\code{\link{grgkw}}) evaluated at \eqn{\gamma=1, \delta=0, \lambda=1}.
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
#' with some tractability advantages. *Statistical Methodology*, *6*(1), 70-81.
#'
#' (Note: Specific gradient formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{grgkw}} (parent distribution gradient),
#' \code{\link{llkw}} (negative log-likelihood for Kw),
#' \code{hskw} (Hessian for Kw, if available),
#' \code{\link{dkw}} (density for Kw),
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{grad}} (for numerical gradient comparison).
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Gradient Evaluation
#'
#' # Generate sample data
#' set.seed(123)
#' n <- 1000
#' true_params <- c(alpha = 2.5, beta = 3.5)
#' data <- rkw(n, alpha = true_params[1], beta = true_params[2])
#'
#' # Evaluate gradient at true parameters
#' grad_true <- grkw(par = true_params, data = data)
#' cat("Gradient at true parameters:\n")
#' print(grad_true)
#' cat("Norm:", sqrt(sum(grad_true^2)), "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 2.5),
#'   c(2.0, 3.0),
#'   c(2.5, 3.5),
#'   c(3.0, 4.0)
#' )
#'
#' grad_norms <- apply(test_params, 1, function(p) {
#'   g <- grkw(p, data)
#'   sqrt(sum(g^2))
#' })
#'
#' results <- data.frame(
#'   Alpha = test_params[, 1],
#'   Beta = test_params[, 2],
#'   Grad_Norm = grad_norms
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Gradient in Optimization
#'
#' # Optimization with analytical gradient
#' fit_with_grad <- optim(
#'   par = c(2, 2),
#'   fn = llkw,
#'   gr = grkw,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE,
#'   control = list(trace = 0)
#' )
#'
#' # Optimization without gradient
#' fit_no_grad <- optim(
#'   par = c(2, 2),
#'   fn = llkw,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE,
#'   control = list(trace = 0)
#' )
#'
#' comparison <- data.frame(
#'   Method = c("With Gradient", "Without Gradient"),
#'   Alpha = c(fit_with_grad$par[1], fit_no_grad$par[1]),
#'   Beta = c(fit_with_grad$par[2], fit_no_grad$par[2]),
#'   NegLogLik = c(fit_with_grad$value, fit_no_grad$value),
#'   Iterations = c(fit_with_grad$counts[1], fit_no_grad$counts[1])
#' )
#' print(comparison, digits = 4, row.names = FALSE)
#'
#'
#' ## Example 3: Verifying Gradient at MLE
#'
#' mle <- fit_with_grad$par
#' names(mle) <- c("alpha", "beta")
#'
#' # At MLE, gradient should be approximately zero
#' gradient_at_mle <- grkw(par = mle, data = data)
#' cat("\nGradient at MLE:\n")
#' print(gradient_at_mle)
#' cat("Max absolute component:", max(abs(gradient_at_mle)), "\n")
#' cat("Gradient norm:", sqrt(sum(gradient_at_mle^2)), "\n")
#'
#'
#' ## Example 4: Numerical vs Analytical Gradient
#'
#' # Manual finite difference gradient
#' numerical_gradient <- function(f, x, data, h = 1e-7) {
#'   grad <- numeric(length(x))
#'   for (i in seq_along(x)) {
#'     x_plus <- x_minus <- x
#'     x_plus[i] <- x[i] + h
#'     x_minus[i] <- x[i] - h
#'     grad[i] <- (f(x_plus, data) - f(x_minus, data)) / (2 * h)
#'   }
#'   return(grad)
#' }
#'
#' # Compare at several points
#' test_points <- rbind(
#'   c(1.5, 2.5),
#'   c(2.0, 3.0),
#'   mle,
#'   c(3.0, 4.0)
#' )
#'
#' cat("\nNumerical vs Analytical Gradient Comparison:\n")
#' for (i in 1:nrow(test_points)) {
#'   grad_analytical <- grkw(par = test_points[i, ], data = data)
#'   grad_numerical <- numerical_gradient(llkw, test_points[i, ], data)
#'
#'   cat(
#'     "\nPoint", i, ": alpha =", test_points[i, 1],
#'     ", beta =", test_points[i, 2], "\n"
#'   )
#'
#'   comparison <- data.frame(
#'     Parameter = c("alpha", "beta"),
#'     Analytical = grad_analytical,
#'     Numerical = grad_numerical,
#'     Abs_Diff = abs(grad_analytical - grad_numerical),
#'     Rel_Error = abs(grad_analytical - grad_numerical) /
#'       (abs(grad_analytical) + 1e-10)
#'   )
#'   print(comparison, digits = 8)
#' }
#'
#'
#' ## Example 5: Gradient Path Visualization
#'
#' # Create grid
#' alpha_grid <- seq(mle[1] - 1, mle[1] + 1, length.out = 20)
#' beta_grid <- seq(mle[2] - 1, mle[2] + 1, length.out = 20)
#' alpha_grid <- alpha_grid[alpha_grid > 0]
#' beta_grid <- beta_grid[beta_grid > 0]
#'
#' # Compute gradient vectors
#' grad_alpha <- matrix(NA, nrow = length(alpha_grid), ncol = length(beta_grid))
#' grad_beta <- matrix(NA, nrow = length(alpha_grid), ncol = length(beta_grid))
#'
#' for (i in seq_along(alpha_grid)) {
#'   for (j in seq_along(beta_grid)) {
#'     g <- grkw(c(alpha_grid[i], beta_grid[j]), data)
#'     grad_alpha[i, j] <- -g[1] # Negative for gradient ascent
#'     grad_beta[i, j] <- -g[2]
#'   }
#' }
#'
#' # Plot gradient field
#'
#' plot(mle[1], mle[2],
#'   pch = 19, col = "#8B0000", cex = 1.5,
#'   xlim = range(alpha_grid), ylim = range(beta_grid),
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Gradient Vector Field", las = 1
#' )
#'
#' # Subsample for clearer visualization
#' step <- 2
#' for (i in seq(1, length(alpha_grid), by = step)) {
#'   for (j in seq(1, length(beta_grid), by = step)) {
#'     arrows(alpha_grid[i], beta_grid[j],
#'       alpha_grid[i] + 0.05 * grad_alpha[i, j],
#'       beta_grid[j] + 0.05 * grad_beta[i, j],
#'       length = 0.05, col = "#2E4057", lwd = 1
#'     )
#'   }
#' }
#'
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True"),
#'   col = c("#8B0000", "#006400"),
#'   pch = c(19, 17), bty = "n"
#' )
#' grid(col = "gray90")
#'
#'
#' ## Example 6: Score Test Statistic
#'
#' # Score test for H0: theta = theta0
#' theta0 <- c(2, 3)
#' score_theta0 <- -grkw(par = theta0, data = data) # Score is negative gradient
#'
#' # Fisher information at theta0 (using Hessian)
#' fisher_info <- hskw(par = theta0, data = data)
#'
#' # Score test statistic
#' score_stat <- t(score_theta0) %*% solve(fisher_info) %*% score_theta0
#' p_value <- pchisq(score_stat, df = 2, lower.tail = FALSE)
#'
#' cat("\nScore Test:\n")
#' cat("H0: alpha = 2, beta = 3\n")
#' cat("Score vector:", score_theta0, "\n")
#' cat("Test statistic:", score_stat, "\n")
#' cat("P-value:", format.pval(p_value, digits = 4), "\n")
#' }
#' @export
grkw <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 2) {
    stop("'par' must be a numeric vector of length 2")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_grkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 7. HESSIAN (hskw)
# ----------------------------------------------------------------------------#

#' @title Hessian Matrix of the Negative Log-Likelihood for the Kw Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize hessian kumaraswamy
#'
#' @description
#' Computes the analytic 2x2 Hessian matrix (matrix of second partial derivatives)
#' of the negative log-likelihood function for the two-parameter Kumaraswamy (Kw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}) and \code{beta}
#' (\eqn{\beta}). The Hessian is useful for estimating standard errors and in
#' optimization algorithms.
#'
#' @param par A numeric vector of length 2 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a 2x2 numeric matrix representing the Hessian matrix of the
#'   negative log-likelihood function, \eqn{-\partial^2 \ell / (\partial \theta_i \partial \theta_j)},
#'   where \eqn{\theta = (\alpha, \beta)}.
#'   Returns a 2x2 matrix populated with \code{NaN} if any parameter values are
#'   invalid according to their constraints, or if any value in \code{data} is
#'   not in the interval (0, 1).
#'
#' @details
#' This function calculates the analytic second partial derivatives of the
#' negative log-likelihood function (\eqn{-\ell(\theta|\mathbf{x})}). The components
#' are the negative of the second derivatives of the log-likelihood \eqn{\ell}
#' (derived from the PDF in \code{\link{dkw}}).
#'
#' Let \eqn{v_i = 1 - x_i^{\alpha}}. The second derivatives of the positive log-likelihood (\eqn{\ell}) are:
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \alpha^2} = -\frac{n}{\alpha^2} -
#' (\beta-1)\sum_{i=1}^{n}\frac{x_i^{\alpha}(\ln(x_i))^2}{v_i^2}
#' }
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \alpha \partial \beta} = -
#' \sum_{i=1}^{n}\frac{x_i^{\alpha}\ln(x_i)}{v_i}
#' }
#' \deqn{
#' \frac{\partial^2 \ell}{\partial \beta^2} = -\frac{n}{\beta^2}
#' }
#' The function returns the Hessian matrix containing the negative of these values.
#'
#' Key properties of the returned matrix:
#' \itemize{
#'   \item Dimensions: 2x2.
#'   \item Symmetry: The matrix is symmetric.
#'   \item Ordering: Rows and columns correspond to the parameters in the order
#'     \eqn{\alpha, \beta}.
#'   \item Content: Analytic second derivatives of the *negative* log-likelihood.
#' }
#' This corresponds to the relevant 2x2 submatrix of the 5x5 GKw Hessian (\code{\link{hsgkw}})
#' evaluated at \eqn{\gamma=1, \delta=0, \lambda=1}.
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
#' with some tractability advantages. *Statistical Methodology*, *6*(1), 70-81.
#'
#' (Note: Specific Hessian formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{hsgkw}} (parent distribution Hessian),
#' \code{\link{llkw}} (negative log-likelihood for Kw),
#' \code{grkw} (gradient for Kw, if available),
#' \code{\link{dkw}} (density for Kw),
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{hessian}} (for numerical Hessian comparison).
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Hessian Evaluation
#'
#' # Generate sample data
#' set.seed(123)
#' n <- 1000
#' true_params <- c(alpha = 2.5, beta = 3.5)
#' data <- rkw(n, alpha = true_params[1], beta = true_params[2])
#'
#' # Evaluate Hessian at true parameters
#' hess_true <- hskw(par = true_params, data = data)
#' cat("Hessian matrix at true parameters:\n")
#' print(hess_true, digits = 4)
#'
#' # Check symmetry
#' cat(
#'   "\nSymmetry check (max |H - H^T|):",
#'   max(abs(hess_true - t(hess_true))), "\n"
#' )
#'
#' ## Example 2: Hessian Properties at MLE
#'
#' # Fit model
#' fit <- optim(
#'   par = c(2, 2),
#'   fn = llkw,
#'   gr = grkw,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("alpha", "beta")
#'
#' # Hessian at MLE
#' hessian_at_mle <- hskw(par = mle, data = data)
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
#' # Observed information matrix (negative Hessian for neg-loglik)
#' obs_info <- hessian_at_mle
#'
#' # Variance-covariance matrix
#' vcov_matrix <- solve(obs_info)
#' cat("\nVariance-Covariance Matrix:\n")
#' print(vcov_matrix, digits = 6)
#'
#' # Standard errors
#' se <- sqrt(diag(vcov_matrix))
#' names(se) <- c("alpha", "beta")
#'
#' # Correlation matrix
#' corr_matrix <- cov2cor(vcov_matrix)
#' cat("\nCorrelation Matrix:\n")
#' print(corr_matrix, digits = 4)
#'
#' # Confidence intervals
#' z_crit <- qnorm(0.975)
#' results <- data.frame(
#'   Parameter = c("alpha", "beta"),
#'   True = true_params,
#'   MLE = mle,
#'   SE = se,
#'   CI_Lower = mle - z_crit * se,
#'   CI_Upper = mle + z_crit * se
#' )
#' print(results, digits = 4)
#'
#' ## Example 4: Determinant and Trace Analysis
#'
#' # Compute at different points
#' test_params <- rbind(
#'   c(1.5, 2.5),
#'   c(2.0, 3.0),
#'   mle,
#'   c(3.0, 4.0)
#' )
#'
#' hess_properties <- data.frame(
#'   Alpha = numeric(),
#'   Beta = numeric(),
#'   Determinant = numeric(),
#'   Trace = numeric(),
#'   Min_Eigenval = numeric(),
#'   Max_Eigenval = numeric(),
#'   Cond_Number = numeric(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (i in 1:nrow(test_params)) {
#'   H <- hskw(par = test_params[i, ], data = data)
#'   eigs <- eigen(H, only.values = TRUE)$values
#'
#'   hess_properties <- rbind(hess_properties, data.frame(
#'     Alpha = test_params[i, 1],
#'     Beta = test_params[i, 2],
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
#' ## Example 5: Curvature Visualization
#'
#' # Create grid around MLE
#' alpha_grid <- seq(mle[1] - 0.5, mle[1] + 0.5, length.out = 30)
#' beta_grid <- seq(mle[2] - 0.5, mle[2] + 0.5, length.out = 30)
#' alpha_grid <- alpha_grid[alpha_grid > 0]
#' beta_grid <- beta_grid[beta_grid > 0]
#'
#' # Compute curvature measures
#' determinant_surface <- matrix(NA,
#'   nrow = length(alpha_grid),
#'   ncol = length(beta_grid)
#' )
#' trace_surface <- matrix(NA,
#'   nrow = length(alpha_grid),
#'   ncol = length(beta_grid)
#' )
#'
#' for (i in seq_along(alpha_grid)) {
#'   for (j in seq_along(beta_grid)) {
#'     H <- hskw(c(alpha_grid[i], beta_grid[j]), data)
#'     determinant_surface[i, j] <- det(H)
#'     trace_surface[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Plot
#'
#' contour(alpha_grid, beta_grid, determinant_surface,
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Hessian Determinant", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(alpha_grid, beta_grid, trace_surface,
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Hessian Trace", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' ## Example 6: Fisher Information and Asymptotic Efficiency
#'
#' # Observed information (at MLE)
#' obs_fisher <- hessian_at_mle
#'
#' # Asymptotic covariance matrix
#' asymp_cov <- solve(obs_fisher)
#'
#' cat("\nAsymptotic Standard Errors:\n")
#' cat("SE(alpha):", sqrt(asymp_cov[1, 1]), "\n")
#' cat("SE(beta):", sqrt(asymp_cov[2, 2]), "\n")
#'
#' # Cramér-Rao Lower Bound
#' cat("\nCramér-Rao Lower Bounds:\n")
#' cat("CRLB(alpha):", sqrt(asymp_cov[1, 1]), "\n")
#' cat("CRLB(beta):", sqrt(asymp_cov[2, 2]), "\n")
#'
#' # Efficiency ellipse (95% confidence region)
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' # Eigendecomposition
#' eig_decomp <- eigen(asymp_cov)
#'
#' # Ellipse points
#' ellipse <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse[i, ] <- mle + sqrt(chi2_val) *
#'     (eig_decomp$vectors %*% diag(sqrt(eig_decomp$values)) %*% v)
#' }
#'
#' # Plot confidence ellipse
#'
#' plot(ellipse[, 1], ellipse[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "95% Confidence Ellipse", las = 1
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CR"),
#'   col = c("#8B0000", "#006400", "#2E4057"),
#'   pch = c(19, 17, NA), lty = c(NA, NA, 1),
#'   lwd = c(NA, NA, 2), bty = "n"
#' )
#' grid(col = "gray90")
#' }
#'
#' @export
hskw <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 2) {
    stop("'par' must be a numeric vector of length 2")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_hskw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}

# ============================================================================#
# END
# ============================================================================#
