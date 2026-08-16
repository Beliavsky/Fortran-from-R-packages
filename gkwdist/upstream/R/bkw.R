# ============================================================================#
# BETA-KUMARASWAMY (BKw) DISTRIBUTION
# ============================================================================#
#
# Wrapper functions for the Beta-Kumaraswamy distribution.
# C++ implementations are in src/bkw.cpp
#
# Functions:
#   - dbkw: Probability density function (PDF)
#   - pbkw: Cumulative distribution function (CDF)
#   - qbkw: Quantile function (inverse CDF)
#   - rbkw: Random number generation
#   - llbkw: Negative log-likelihood
#   - grbkw: Gradient of negative log-likelihood
#   - hsbkw: Hessian of negative log-likelihood
# ============================================================================#


# ----------------------------------------------------------------------------#
# 1. DENSITY FUNCTION (dbkw)
# ----------------------------------------------------------------------------#

#' @title Density of the Beta-Kumaraswamy (BKw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution density
#'
#' @description
#' Computes the probability density function (PDF) for the Beta-Kumaraswamy
#' (BKw) distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{gamma} (\eqn{\gamma}), and \code{delta} (\eqn{\delta}).
#' This distribution is defined on the interval (0, 1).
#'
#' @param x Vector of quantiles (values between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param log Logical; if \code{TRUE}, the logarithm of the density is
#'   returned (\eqn{\log(f(x))}). Default: \code{FALSE}.
#'
#' @return A vector of density values (\eqn{f(x)}) or log-density values
#'   (\eqn{\log(f(x))}). The length of the result is determined by the recycling
#'   rule applied to the arguments (\code{x}, \code{alpha}, \code{beta},
#'   \code{gamma}, \code{delta}). Returns \code{0} (or \code{-Inf} if
#'   \code{log = TRUE}) for \code{x} outside the interval (0, 1), or
#'   \code{NaN} if parameters are invalid (e.g., \code{alpha <= 0}, \code{beta <= 0},
#'   \code{gamma <= 0}, \code{delta < 0}).
#'
#' @details
#' The probability density function (PDF) of the Beta-Kumaraswamy (BKw)
#' distribution is given by:
#' \deqn{
#' f(x; \alpha, \beta, \gamma, \delta) = \frac{\alpha \beta}{B(\gamma, \delta+1)} x^{\alpha - 1} \bigl(1 - x^\alpha\bigr)^{\beta(\delta+1) - 1} \bigl[1 - \bigl(1 - x^\alpha\bigr)^\beta\bigr]^{\gamma - 1}
#' }
#' for \eqn{0 < x < 1}, where \eqn{B(a,b)} is the Beta function
#' (\code{\link[base]{beta}}).
#'
#' The BKw distribution is a special case of the five-parameter
#' Generalized Kumaraswamy (GKw) distribution (\code{\link{dgkw}}) obtained
#' by setting the parameter \eqn{\lambda = 1}.
#' Numerical evaluation is performed using algorithms similar to those for `dgkw`,
#' ensuring stability.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' @seealso
#' \code{\link{dgkw}} (parent distribution density),
#' \code{\link{pbkw}}, \code{\link{qbkw}}, \code{\link{rbkw}} (other BKw functions),
#'
#' @examples
#' \donttest{
#' # Example values
#' x_vals <- c(0.2, 0.5, 0.8)
#' alpha_par <- 2.0
#' beta_par <- 1.5
#' gamma_par <- 1.0 # Equivalent to Kw when gamma=1
#' delta_par <- 0.5
#'
#' # Calculate density
#' densities <- dbkw(x_vals, alpha_par, beta_par, gamma_par, delta_par)
#' print(densities)
#'
#' # Calculate log-density
#' log_densities <- dbkw(x_vals, alpha_par, beta_par, gamma_par, delta_par,
#'   log = TRUE
#' )
#' print(log_densities)
#' # Check: should match log(densities)
#' print(log(densities))
#'
#' # Compare with dgkw setting lambda = 1
#' densities_gkw <- dgkw(x_vals, alpha_par, beta_par,
#'   gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print(paste("Max difference:", max(abs(densities - densities_gkw)))) # Should be near zero
#'
#' # Plot the density for different gamma values
#' curve_x <- seq(0.01, 0.99, length.out = 200)
#' curve_y1 <- dbkw(curve_x, alpha = 2, beta = 3, gamma = 0.5, delta = 1)
#' curve_y2 <- dbkw(curve_x, alpha = 2, beta = 3, gamma = 1.0, delta = 1)
#' curve_y3 <- dbkw(curve_x, alpha = 2, beta = 3, gamma = 2.0, delta = 1)
#'
#' plot(curve_x, curve_y1,
#'   type = "l", main = "BKw Density Examples (alpha=2, beta=3, delta=1)",
#'   xlab = "x", ylab = "f(x)", col = "blue", ylim = range(0, curve_y1, curve_y2, curve_y3)
#' )
#' lines(curve_x, curve_y2, col = "red")
#' lines(curve_x, curve_y3, col = "green")
#' legend("topright",
#'   legend = c("gamma=0.5", "gamma=1.0", "gamma=2.0"),
#'   col = c("blue", "red", "green"), lty = 1, bty = "n"
#' )
#' }
#'
#' @export
dbkw <- function(x, alpha = 1, beta = 1, gamma = 1, delta = 0, log = FALSE) {
  # Input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
  }
  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("'log' must be a single logical value")
  }

  # Call C++ implementation
  .Call("_gkwdist_dbkw",
    as.numeric(x),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(gamma),
    as.numeric(delta),
    as.logical(log),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 2. DISTRIBUTION FUNCTION (pbkw)
# ----------------------------------------------------------------------------#

#' @title Cumulative Distribution Function (CDF) of the Beta-Kumaraswamy (BKw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution cumulative
#'
#' @description
#' Computes the cumulative distribution function (CDF), \eqn{P(X \le q)}, for the
#' Beta-Kumaraswamy (BKw) distribution with parameters \code{alpha} (\eqn{\alpha}),
#' \code{beta} (\eqn{\beta}), \code{gamma} (\eqn{\gamma}), and \code{delta}
#' (\eqn{\delta}). This distribution is defined on the interval (0, 1) and is
#' a special case of the Generalized Kumaraswamy (GKw) distribution where
#' \eqn{\lambda = 1}.
#'
#' @param q Vector of quantiles (values generally between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}, otherwise, \eqn{P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of probabilities, \eqn{F(q)}, or their logarithms/complements
#'   depending on \code{lower.tail} and \code{log.p}. The length of the result
#'   is determined by the recycling rule applied to the arguments (\code{q},
#'   \code{alpha}, \code{beta}, \code{gamma}, \code{delta}). Returns \code{0}
#'   (or \code{-Inf} if \code{log.p = TRUE}) for \code{q <= 0} and \code{1}
#'   (or \code{0} if \code{log.p = TRUE}) for \code{q >= 1}. Returns \code{NaN}
#'   for invalid parameters.
#'
#' @details
#' The Beta-Kumaraswamy (BKw) distribution is a special case of the
#' five-parameter Generalized Kumaraswamy distribution (\code{\link{pgkw}})
#' obtained by setting the shape parameter \eqn{\lambda = 1}.
#'
#' The CDF of the GKw distribution is \eqn{F_{GKw}(q) = I_{y(q)}(\gamma, \delta+1)},
#' where \eqn{y(q) = [1-(1-q^{\alpha})^{\beta}]^{\lambda}} and \eqn{I_x(a,b)}
#' is the regularized incomplete beta function (\code{\link[stats]{pbeta}}).
#' Setting \eqn{\lambda=1} simplifies \eqn{y(q)} to \eqn{1 - (1 - q^\alpha)^\beta},
#' yielding the BKw CDF:
#' \deqn{
#' F(q; \alpha, \beta, \gamma, \delta) = I_{1 - (1 - q^\alpha)^\beta}(\gamma, \delta+1)
#' }
#' This is evaluated using the \code{\link[stats]{pbeta}} function.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' @seealso
#' \code{\link{pgkw}} (parent distribution CDF),
#' \code{\link{dbkw}}, \code{\link{qbkw}}, \code{\link{rbkw}} (other BKw functions),
#' \code{\link[stats]{pbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' q_vals <- c(0.2, 0.5, 0.8)
#' alpha_par <- 2.0
#' beta_par <- 1.5
#' gamma_par <- 1.0
#' delta_par <- 0.5
#'
#' # Calculate CDF P(X <= q)
#' probs <- pbkw(q_vals, alpha_par, beta_par, gamma_par, delta_par)
#' print(probs)
#'
#' # Calculate upper tail P(X > q)
#' probs_upper <- pbkw(q_vals, alpha_par, beta_par, gamma_par, delta_par,
#'   lower.tail = FALSE
#' )
#' print(probs_upper)
#' # Check: probs + probs_upper should be 1
#' print(probs + probs_upper)
#'
#' # Calculate log CDF
#' logs <- pbkw(q_vals, alpha_par, beta_par, gamma_par, delta_par,
#'   log.p = TRUE
#' )
#' print(logs)
#' # Check: should match log(probs)
#' print(log(probs))
#'
#' # Compare with pgkw setting lambda = 1
#' probs_gkw <- pgkw(q_vals, alpha_par, beta_par,
#'   gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print(paste("Max difference:", max(abs(probs - probs_gkw)))) # Should be near zero
#'
#' # Plot the CDF
#' curve_q <- seq(0.01, 0.99, length.out = 200)
#' curve_p <- pbkw(curve_q, alpha = 2, beta = 3, gamma = 0.5, delta = 1)
#' plot(curve_q, curve_p,
#'   type = "l", main = "BKw CDF Example",
#'   xlab = "q", ylab = "F(q)", col = "blue", ylim = c(0, 1)
#' )
#' }
#'
#' @export
pbkw <- function(q, alpha = 1, beta = 1, gamma = 1, delta = 0, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(q)) stop("'q' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
  }
  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }
  if (!is.logical(lower.tail) || length(lower.tail) != 1) {
    stop("'lower.tail' must be a single logical value")
  }
  if (!is.logical(log.p) || length(log.p) != 1) {
    stop("'log.p' must be a single logical value")
  }

  # Call C++ implementation
  .Call("_gkwdist_pbkw",
    as.numeric(q),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(gamma),
    as.numeric(delta),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 3. QUANTILE FUNCTION (qbkw)
# ----------------------------------------------------------------------------#


#' @title Quantile Function of the Beta-Kumaraswamy (BKw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution quantile
#'
#' @description
#' Computes the quantile function (inverse CDF) for the Beta-Kumaraswamy (BKw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{gamma} (\eqn{\gamma}), and \code{delta} (\eqn{\delta}).
#' It finds the value \code{q} such that \eqn{P(X \le q) = p}. This distribution
#' is a special case of the Generalized Kumaraswamy (GKw) distribution where
#' the parameter \eqn{\lambda = 1}.
#'
#' @param p Vector of probabilities (values between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{p = P(X \le q)},
#'   otherwise, probabilities are \eqn{p = P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of quantiles corresponding to the given probabilities \code{p}.
#'   The length of the result is determined by the recycling rule applied to
#'   the arguments (\code{p}, \code{alpha}, \code{beta}, \code{gamma}, \code{delta}).
#'   Returns:
#'   \itemize{
#'     \item \code{0} for \code{p = 0} (or \code{p = -Inf} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{1} for \code{p = 1} (or \code{p = 0} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{NaN} for \code{p < 0} or \code{p > 1} (or corresponding log scale).
#'     \item \code{NaN} for invalid parameters (e.g., \code{alpha <= 0},
#'           \code{beta <= 0}, \code{gamma <= 0}, \code{delta < 0}).
#'   }
#'   Boundary return values are adjusted accordingly for \code{lower.tail = FALSE}.
#'
#' @details
#' The quantile function \eqn{Q(p)} is the inverse of the CDF \eqn{F(q)}. The CDF
#' for the BKw (\eqn{\lambda=1}) distribution is \eqn{F(q) = I_{y(q)}(\gamma, \delta+1)},
#' where \eqn{y(q) = 1 - (1 - q^\alpha)^\beta} and \eqn{I_z(a,b)} is the
#' regularized incomplete beta function (see \code{\link{pbkw}}).
#'
#' To find the quantile \eqn{q}, we first invert the outer Beta part: let
#' \eqn{y = I^{-1}_{p}(\gamma, \delta+1)}, where \eqn{I^{-1}_p(a,b)} is the
#' inverse of the regularized incomplete beta function, computed via
#' \code{\link[stats]{qbeta}}. Then, we invert the inner Kumaraswamy part:
#' \eqn{y = 1 - (1 - q^\alpha)^\beta}, which leads to \eqn{q = \{1 - (1-y)^{1/\beta}\}^{1/\alpha}}.
#' Substituting \eqn{y} gives the quantile function:
#' \deqn{
#' Q(p) = \left\{ 1 - \left[ 1 - I^{-1}_{p}(\gamma, \delta+1) \right]^{1/\beta} \right\}^{1/\alpha}
#' }
#' The function uses this formula, calculating \eqn{I^{-1}_{p}(\gamma, \delta+1)}
#' via \code{qbeta(p, gamma, delta + 1, ...)} while respecting the
#' \code{lower.tail} and \code{log.p} arguments.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' @seealso
#' \code{\link{qgkw}} (parent distribution quantile function),
#' \code{\link{dbkw}}, \code{\link{pbkw}}, \code{\link{rbkw}} (other BKw functions),
#' \code{\link[stats]{qbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' p_vals <- c(0.1, 0.5, 0.9)
#' alpha_par <- 2.0
#' beta_par <- 1.5
#' gamma_par <- 1.0
#' delta_par <- 0.5
#'
#' # Calculate quantiles
#' quantiles <- qbkw(p_vals, alpha_par, beta_par, gamma_par, delta_par)
#' print(quantiles)
#'
#' # Calculate quantiles for upper tail probabilities P(X > q) = p
#' quantiles_upper <- qbkw(p_vals, alpha_par, beta_par, gamma_par, delta_par,
#'   lower.tail = FALSE
#' )
#' print(quantiles_upper)
#' # Check: qbkw(p, ..., lt=F) == qbkw(1-p, ..., lt=T)
#' print(qbkw(1 - p_vals, alpha_par, beta_par, gamma_par, delta_par))
#'
#' # Calculate quantiles from log probabilities
#' log.p_vals <- log(p_vals)
#' quantiles_logp <- qbkw(log.p_vals, alpha_par, beta_par, gamma_par, delta_par,
#'   log.p = TRUE
#' )
#' print(quantiles_logp)
#' # Check: should match original quantiles
#' print(quantiles)
#'
#' # Compare with qgkw setting lambda = 1
#' quantiles_gkw <- qgkw(p_vals, alpha_par, beta_par,
#'   gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print(paste("Max difference:", max(abs(quantiles - quantiles_gkw)))) # Should be near zero
#'
#' # Verify inverse relationship with pbkw
#' p_check <- 0.75
#' q_calc <- qbkw(p_check, alpha_par, beta_par, gamma_par, delta_par)
#' p_recalc <- pbkw(q_calc, alpha_par, beta_par, gamma_par, delta_par)
#' print(paste("Original p:", p_check, " Recalculated p:", p_recalc))
#' # abs(p_check - p_recalc) < 1e-9 # Should be TRUE
#'
#' # Boundary conditions
#' print(qbkw(c(0, 1), alpha_par, beta_par, gamma_par, delta_par)) # Should be 0, 1
#' print(qbkw(c(-Inf, 0), alpha_par, beta_par, gamma_par, delta_par, log.p = TRUE)) # Should be 0, 1
#' }
#'
#' @export
qbkw <- function(p, alpha = 1, beta = 1, gamma = 1, delta = 0, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(p)) stop("'p' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
  }
  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
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
  .Call("_gkwdist_qbkw",
    as.numeric(p),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(gamma),
    as.numeric(delta),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 4. RANDOM GENERATION (rbkw)
# ----------------------------------------------------------------------------#

#' @title Random Number Generation for the Beta-Kumaraswamy (BKw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution random
#'
#' @description
#' Generates random deviates from the Beta-Kumaraswamy (BKw) distribution
#' with parameters \code{alpha} (\eqn{\alpha}), \code{beta} (\eqn{\beta}),
#' \code{gamma} (\eqn{\gamma}), and \code{delta} (\eqn{\delta}). This distribution
#' is a special case of the Generalized Kumaraswamy (GKw) distribution where
#' the parameter \eqn{\lambda = 1}.
#'
#' @param n Number of observations. If \code{length(n) > 1}, the length is
#'   taken to be the number required. Must be a non-negative integer.
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param gamma Shape parameter \code{gamma} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#'
#' @return A vector of length \code{n} containing random deviates from the BKw
#'   distribution. The length of the result is determined by \code{n} and the
#'   recycling rule applied to the parameters (\code{alpha}, \code{beta},
#'   \code{gamma}, \code{delta}). Returns \code{NaN} if parameters
#'   are invalid (e.g., \code{alpha <= 0}, \code{beta <= 0}, \code{gamma <= 0},
#'   \code{delta < 0}).
#'
#' @details
#' The generation method uses the relationship between the GKw distribution and the
#' Beta distribution. The general procedure for GKw (\code{\link{rgkw}}) is:
#' If \eqn{W \sim \mathrm{Beta}(\gamma, \delta+1)}, then
#' \eqn{X = \{1 - [1 - W^{1/\lambda}]^{1/\beta}\}^{1/\alpha}} follows the
#' GKw(\eqn{\alpha, \beta, \gamma, \delta, \lambda}) distribution.
#'
#' For the BKw distribution, \eqn{\lambda=1}. Therefore, the algorithm simplifies to:
#' \enumerate{
#'   \item Generate \eqn{V \sim \mathrm{Beta}(\gamma, \delta+1)} using
#'         \code{\link[stats]{rbeta}}.
#'   \item Compute the BKw variate \eqn{X = \{1 - (1 - V)^{1/\beta}\}^{1/\alpha}}.
#' }
#' This procedure is implemented efficiently, handling parameter recycling as needed.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
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
#' \code{\link{dbkw}}, \code{\link{pbkw}}, \code{\link{qbkw}} (other BKw functions),
#' \code{\link[stats]{rbeta}}
#'
#' @examples
#' \donttest{
#' set.seed(2026) # for reproducibility
#'
#' # Generate 1000 random values from a specific BKw distribution
#' alpha_par <- 2.0
#' beta_par <- 1.5
#' gamma_par <- 1.0
#' delta_par <- 0.5
#'
#' x_sample_bkw <- rbkw(1000,
#'   alpha = alpha_par, beta = beta_par,
#'   gamma = gamma_par, delta = delta_par
#' )
#' summary(x_sample_bkw)
#'
#' # Histogram of generated values compared to theoretical density
#' hist(x_sample_bkw,
#'   breaks = 30, freq = FALSE, # freq=FALSE for density
#'   main = "Histogram of BKw Sample", xlab = "x", ylim = c(0, 2.5)
#' )
#' curve(
#'   dbkw(x,
#'     alpha = alpha_par, beta = beta_par, gamma = gamma_par,
#'     delta = delta_par
#'   ),
#'   add = TRUE, col = "red", lwd = 2, n = 201
#' )
#' legend("topright", legend = "Theoretical PDF", col = "red", lwd = 2, bty = "n")
#'
#' # Comparing empirical and theoretical quantiles (Q-Q plot)
#' prob_points <- seq(0.01, 0.99, by = 0.01)
#' theo_quantiles <- qbkw(prob_points,
#'   alpha = alpha_par, beta = beta_par,
#'   gamma = gamma_par, delta = delta_par
#' )
#' emp_quantiles <- quantile(x_sample_bkw, prob_points, type = 7)
#'
#' plot(theo_quantiles, emp_quantiles,
#'   pch = 16, cex = 0.8,
#'   main = "Q-Q Plot for BKw Distribution",
#'   xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles (n=1000)"
#' )
#' abline(a = 0, b = 1, col = "blue", lty = 2)
#'
#' # Compare summary stats with rgkw(..., lambda=1, ...)
#' # Note: individual values will differ due to randomness
#' x_sample_gkw <- rgkw(1000,
#'   alpha = alpha_par, beta = beta_par, gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print("Summary stats for rbkw sample:")
#' print(summary(x_sample_bkw))
#' print("Summary stats for rgkw(lambda=1) sample:")
#' print(summary(x_sample_gkw)) # Should be similar
#' }
#'
#' @export
rbkw <- function(n, alpha = 1, beta = 1, gamma = 1, delta = 0) {
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
  if (!is.numeric(gamma) || any(gamma <= 0)) {
    stop("'gamma' must be positive")
  }
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }

  # Call C++ implementation
  .Call("_gkwdist_rbkw",
    as.integer(n),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(gamma),
    as.numeric(delta),
    PACKAGE = "gkwdist"
  )
}


# ============================================================================#
# MAXIMUM LIKELIHOOD ESTIMATION FUNCTIONS
# ============================================================================#

# ----------------------------------------------------------------------------#
# 5. NEGATIVE LOG-LIKELIHOOD (llbkw)
# ----------------------------------------------------------------------------#

#' @title Negative Log-Likelihood for Beta-Kumaraswamy (BKw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize
#'
#' @description
#' Computes the negative log-likelihood function for the Beta-Kumaraswamy (BKw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{gamma} (\eqn{\gamma}), and \code{delta} (\eqn{\delta}),
#' given a vector of observations. This distribution is the special case of the
#' Generalized Kumaraswamy (GKw) distribution where \eqn{\lambda = 1}. This function
#' is typically used for maximum likelihood estimation via numerical optimization.
#'
#' @param par A numeric vector of length 4 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
#'   \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a single \code{double} value representing the negative
#'   log-likelihood (\eqn{-\ell(\theta|\mathbf{x})}). Returns \code{Inf}
#'   if any parameter values in \code{par} are invalid according to their
#'   constraints, or if any value in \code{data} is not in the interval (0, 1).
#'
#' @details
#' The Beta-Kumaraswamy (BKw) distribution is the GKw distribution (\code{\link{dgkw}})
#' with \eqn{\lambda=1}. Its probability density function (PDF) is:
#' \deqn{
#' f(x | \theta) = \frac{\alpha \beta}{B(\gamma, \delta+1)} x^{\alpha - 1} \bigl(1 - x^\alpha\bigr)^{\beta(\delta+1) - 1} \bigl[1 - \bigl(1 - x^\alpha\bigr)^\beta\bigr]^{\gamma - 1}
#' }
#' for \eqn{0 < x < 1}, \eqn{\theta = (\alpha, \beta, \gamma, \delta)}, and \eqn{B(a,b)}
#' is the Beta function (\code{\link[base]{beta}}).
#' The log-likelihood function \eqn{\ell(\theta | \mathbf{x})} for a sample
#' \eqn{\mathbf{x} = (x_1, \dots, x_n)} is \eqn{\sum_{i=1}^n \ln f(x_i | \theta)}:
#' \deqn{
#' \ell(\theta | \mathbf{x}) = n[\ln(\alpha) + \ln(\beta) - \ln B(\gamma, \delta+1)]
#' + \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta(\delta+1)-1)\ln(v_i) + (\gamma-1)\ln(w_i)]
#' }
#' where:
#' \itemize{
#'   \item \eqn{v_i = 1 - x_i^{\alpha}}
#'   \item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
#' }
#' This function computes and returns the *negative* log-likelihood, \eqn{-\ell(\theta|\mathbf{x})},
#' suitable for minimization using optimization routines like \code{\link[stats]{optim}}.
#' Numerical stability is maintained similarly to \code{\link{llgkw}}.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#'
#' @seealso
#' \code{\link{llgkw}} (parent distribution negative log-likelihood),
#' \code{\link{dbkw}}, \code{\link{pbkw}}, \code{\link{qbkw}}, \code{\link{rbkw}},
#' \code{grbkw} (gradient, if available),
#' \code{hsbkw} (Hessian, if available),
#' \code{\link[stats]{optim}}, \code{\link[base]{lbeta}}
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Log-Likelihood Evaluation
#' # Generate sample data
#' set.seed(2203)
#' n <- 1000
#' true_params <- c(alpha = 2.0, beta = 1.5, gamma = 1.5, delta = 0.5)
#' data <- rbkw(n,
#'   alpha = true_params[1], beta = true_params[2],
#'   gamma = true_params[3], delta = true_params[4]
#' )
#'
#' # Evaluate negative log-likelihood at true parameters
#' nll_true <- llbkw(par = true_params, data = data)
#' cat("Negative log-likelihood at true parameters:", nll_true, "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 1.0, 1.0, 0.3),
#'   c(2.0, 1.5, 1.5, 0.5),
#'   c(2.5, 2.0, 2.0, 0.7)
#' )
#'
#' nll_values <- apply(test_params, 1, function(p) llbkw(p, data))
#' results <- data.frame(
#'   Alpha = test_params[, 1],
#'   Beta = test_params[, 2],
#'   Gamma = test_params[, 3],
#'   Delta = test_params[, 4],
#'   NegLogLik = nll_values
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Maximum Likelihood Estimation
#'
#' # Optimization using BFGS with no analytical gradient
#' fit <- optim(
#'   par = c(0.5, 1, 1.1, 0.3),
#'   fn = llbkw,
#'   # gr = grbkw,
#'   data = data,
#'   method = "BFGS",
#'   control = list(maxit = 2000),
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("alpha", "beta", "gamma", "delta")
#' se <- sqrt(diag(solve(fit$hessian)))
#'
#' results <- data.frame(
#'   Parameter = c("alpha", "beta", "gamma", "delta"),
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
#' start_params <- c(1.8, 1.2, 1.1, 0.3)
#'
#' comparison <- data.frame(
#'   Method = character(),
#'   Alpha = numeric(),
#'   Beta = numeric(),
#'   Gamma = numeric(),
#'   Delta = numeric(),
#'   NegLogLik = numeric(),
#'   Convergence = integer(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (method in methods) {
#'   if (method %in% c("BFGS", "CG")) {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llbkw,
#'       gr = grbkw,
#'       data = data,
#'       method = method
#'     )
#'   } else if (method == "L-BFGS-B") {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llbkw,
#'       gr = grbkw,
#'       data = data,
#'       method = method,
#'       lower = c(0.01, 0.01, 0.01, 0.01),
#'       upper = c(100, 100, 100, 100)
#'     )
#'   } else {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llbkw,
#'       data = data,
#'       method = method
#'     )
#'   }
#'
#'   comparison <- rbind(comparison, data.frame(
#'     Method = method,
#'     Alpha = fit_temp$par[1],
#'     Beta = fit_temp$par[2],
#'     Gamma = fit_temp$par[3],
#'     Delta = fit_temp$par[4],
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
#' # Test H0: delta = 0.5 vs H1: delta free
#' loglik_full <- -fit$value
#'
#' restricted_ll <- function(params_restricted, data, delta_fixed) {
#'   llbkw(par = c(
#'     params_restricted[1], params_restricted[2],
#'     params_restricted[3], delta_fixed
#'   ), data = data)
#' }
#'
#' fit_restricted <- optim(
#'   par = mle[1:3],
#'   fn = restricted_ll,
#'   data = data,
#'   delta_fixed = 0.5,
#'   method = "Nelder-Mead"
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
#' # Profile for alpha
#' alpha_grid <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = 50)
#' alpha_grid <- alpha_grid[alpha_grid > 0]
#' profile_ll_alpha <- numeric(length(alpha_grid))
#'
#' for (i in seq_along(alpha_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-1],
#'     fn = function(p) llbkw(c(alpha_grid[i], p), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_alpha[i] <- -profile_fit$value
#' }
#'
#' # Profile for beta
#' beta_grid <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = 50)
#' beta_grid <- beta_grid[beta_grid > 0]
#' profile_ll_beta <- numeric(length(beta_grid))
#'
#' for (i in seq_along(beta_grid)) {
#'   profile_fit <- optim(
#'     par = c(mle[1], mle[3], mle[4]),
#'     fn = function(p) llbkw(c(mle[1], beta_grid[i], p[1], p[2]), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_beta[i] <- -profile_fit$value
#' }
#'
#' # Profile for gamma
#' gamma_grid <- seq(mle[3] - 1.5, mle[3] + 1.5, length.out = 50)
#' gamma_grid <- gamma_grid[gamma_grid > 0]
#' profile_ll_gamma <- numeric(length(gamma_grid))
#'
#' for (i in seq_along(gamma_grid)) {
#'   profile_fit <- optim(
#'     par = c(mle[1], mle[2], mle[4]),
#'     fn = function(p) llbkw(c(p[1], mle[2], gamma_grid[i], p[2]), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_gamma[i] <- -profile_fit$value
#' }
#'
#' # Profile for delta
#' delta_grid <- seq(mle[4] - 1.5, mle[4] + 1.5, length.out = 50)
#' delta_grid <- delta_grid[delta_grid > 0]
#' profile_ll_delta <- numeric(length(delta_grid))
#'
#' for (i in seq_along(delta_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-4],
#'     fn = function(p) llbkw(c(p[1], p[2], p[3], delta_grid[i]), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_delta[i] <- -profile_fit$value
#' }
#'
#' # 95% confidence threshold
#' chi_crit <- qchisq(0.95, df = 1)
#' threshold <- max(profile_ll_alpha) - chi_crit / 2
#'
#' # Plot all profiles
#'
#' plot(alpha_grid, profile_ll_alpha,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", alpha)), las = 1
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
#' plot(beta_grid, profile_ll_beta,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(beta), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", beta)), las = 1
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
#' plot(gamma_grid, profile_ll_gamma,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", gamma)), las = 1
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
#' plot(delta_grid, profile_ll_delta,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(delta), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", delta)), las = 1
#' )
#' abline(v = mle[4], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[4], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3, lwd = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#808080"),
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.8
#' )
#' grid(col = "gray90")
#'
#'
#' ## Example 6: 2D Log-Likelihood Surfaces (Selected pairs)
#'
#' # Create 2D grids with wider range (±1.5)
#' alpha_2d <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = round(n / 25))
#' beta_2d <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = round(n / 25))
#' gamma_2d <- seq(mle[3] - 1.5, mle[3] + 1.5, length.out = round(n / 25))
#' delta_2d <- seq(mle[4] - 1.5, mle[4] + 1.5, length.out = round(n / 25))
#'
#' alpha_2d <- alpha_2d[alpha_2d > 0]
#' beta_2d <- beta_2d[beta_2d > 0]
#' gamma_2d <- gamma_2d[gamma_2d > 0]
#' delta_2d <- delta_2d[delta_2d > 0]
#'
#' # Compute selected log-likelihood surfaces
#' ll_surface_ab <- matrix(NA, nrow = length(alpha_2d), ncol = length(beta_2d))
#' ll_surface_ag <- matrix(NA, nrow = length(alpha_2d), ncol = length(gamma_2d))
#' ll_surface_bd <- matrix(NA, nrow = length(beta_2d), ncol = length(delta_2d))
#'
#' # Alpha vs Beta
#' for (i in seq_along(alpha_2d)) {
#'   for (j in seq_along(beta_2d)) {
#'     ll_surface_ab[i, j] <- -llbkw(c(alpha_2d[i], beta_2d[j], mle[3], mle[4]), data)
#'   }
#' }
#'
#' # Alpha vs Gamma
#' for (i in seq_along(alpha_2d)) {
#'   for (j in seq_along(gamma_2d)) {
#'     ll_surface_ag[i, j] <- -llbkw(c(alpha_2d[i], mle[2], gamma_2d[j], mle[4]), data)
#'   }
#' }
#'
#' # Beta vs Delta
#' for (i in seq_along(beta_2d)) {
#'   for (j in seq_along(delta_2d)) {
#'     ll_surface_bd[i, j] <- -llbkw(c(mle[1], beta_2d[i], mle[3], delta_2d[j]), data)
#'   }
#' }
#'
#' # Confidence region levels
#' max_ll_ab <- max(ll_surface_ab, na.rm = TRUE)
#' max_ll_ag <- max(ll_surface_ag, na.rm = TRUE)
#' max_ll_bd <- max(ll_surface_bd, na.rm = TRUE)
#'
#' levels_95_ab <- max_ll_ab - qchisq(0.95, df = 2) / 2
#' levels_95_ag <- max_ll_ag - qchisq(0.95, df = 2) / 2
#' levels_95_bd <- max_ll_bd - qchisq(0.95, df = 2) / 2
#'
#' # Plot selected surfaces
#'
#' # Alpha vs Beta
#' contour(alpha_2d, beta_2d, ll_surface_ab,
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Alpha vs Beta", las = 1,
#'   levels = seq(min(ll_surface_ab, na.rm = TRUE), max_ll_ab, length.out = 20),
#'   col = "#2E4057", lwd = 1
#' )
#' contour(alpha_2d, beta_2d, ll_surface_ab,
#'   levels = levels_95_ab, col = "#FF6347", lwd = 2.5, lty = 1, add = TRUE
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Alpha vs Gamma
#' contour(alpha_2d, gamma_2d, ll_surface_ag,
#'   xlab = expression(alpha), ylab = expression(gamma),
#'   main = "Alpha vs Gamma", las = 1,
#'   levels = seq(min(ll_surface_ag, na.rm = TRUE), max_ll_ag, length.out = 20),
#'   col = "#2E4057", lwd = 1
#' )
#' contour(alpha_2d, gamma_2d, ll_surface_ag,
#'   levels = levels_95_ag, col = "#FF6347", lwd = 2.5, lty = 1, add = TRUE
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Beta vs Delta
#' contour(beta_2d, delta_2d, ll_surface_bd,
#'   xlab = expression(beta), ylab = expression(delta),
#'   main = "Beta vs Delta", las = 1,
#'   levels = seq(min(ll_surface_bd, na.rm = TRUE), max_ll_bd, length.out = 20),
#'   col = "#2E4057", lwd = 1
#' )
#' contour(beta_2d, delta_2d, ll_surface_bd,
#'   levels = levels_95_bd, col = "#FF6347", lwd = 2.5, lty = 1, add = TRUE
#' )
#' points(mle[2], mle[4], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[4], pch = 17, col = "#006400", cex = 1.5)
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
llbkw <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 4) {
    stop("'par' must be a numeric vector of length 4")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_llbkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 6. GRADIENT (grbkw)
# ----------------------------------------------------------------------------#

#' @title Gradient of the Negative Log-Likelihood for the BKw Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize gradient
#'
#' @description
#' Computes the gradient vector (vector of first partial derivatives) of the
#' negative log-likelihood function for the Beta-Kumaraswamy (BKw) distribution
#' with parameters \code{alpha} (\eqn{\alpha}), \code{beta} (\eqn{\beta}),
#' \code{gamma} (\eqn{\gamma}), and \code{delta} (\eqn{\delta}). This distribution
#' is the special case of the Generalized Kumaraswamy (GKw) distribution where
#' \eqn{\lambda = 1}. The gradient is typically used in optimization algorithms
#' for maximum likelihood estimation.
#'
#' @param par A numeric vector of length 4 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
#'   \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a numeric vector of length 4 containing the partial derivatives
#'   of the negative log-likelihood function \eqn{-\ell(\theta | \mathbf{x})} with
#'   respect to each parameter:
#'   \eqn{(-\partial \ell/\partial \alpha, -\partial \ell/\partial \beta, -\partial \ell/\partial \gamma, -\partial \ell/\partial \delta)}.
#'   Returns a vector of \code{NaN} if any parameter values are invalid according
#'   to their constraints, or if any value in \code{data} is not in the
#'   interval (0, 1).
#'
#' @details
#' The components of the gradient vector of the negative log-likelihood
#' (\eqn{-\nabla \ell(\theta | \mathbf{x})}) for the BKw (\eqn{\lambda=1}) model are:
#'
#' \deqn{
#' -\frac{\partial \ell}{\partial \alpha} = -\frac{n}{\alpha} - \sum_{i=1}^{n}\ln(x_i)
#' + \sum_{i=1}^{n}\left[x_i^{\alpha} \ln(x_i) \left(\frac{\beta(\delta+1)-1}{v_i} -
#' \frac{(\gamma-1) \beta v_i^{\beta-1}}{w_i}\right)\right]
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \beta} = -\frac{n}{\beta} - (\delta+1)\sum_{i=1}^{n}\ln(v_i)
#' + \sum_{i=1}^{n}\left[\frac{(\gamma-1) v_i^{\beta} \ln(v_i)}{w_i}\right]
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \gamma} = n[\psi(\gamma) - \psi(\gamma+\delta+1)] -
#' \sum_{i=1}^{n}\ln(w_i)
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \delta} = n[\psi(\delta+1) - \psi(\gamma+\delta+1)] -
#' \beta\sum_{i=1}^{n}\ln(v_i)
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{v_i = 1 - x_i^{\alpha}}
#'   \item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
#'   \item \eqn{\psi(\cdot)} is the digamma function (\code{\link[base]{digamma}}).
#' }
#' These formulas represent the derivatives of \eqn{-\ell(\theta)}, consistent with
#' minimizing the negative log-likelihood. They correspond to the general GKw
#' gradient (\code{\link{grgkw}}) components for \eqn{\alpha, \beta, \gamma, \delta}
#' evaluated at \eqn{\lambda=1}. Note that the component for \eqn{\lambda} is omitted.
#' Numerical stability is maintained through careful implementation.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' (Note: Specific gradient formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{grgkw}} (parent distribution gradient),
#' \code{\link{llbkw}} (negative log-likelihood for BKw),
#' \code{\link{hsbkw}} (Hessian for BKw, if available),
#' \code{\link{dbkw}} (density for BKw),
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{grad}} (for numerical gradient comparison),
#' \code{\link[base]{digamma}}.
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Gradient Evaluation
#' # Generate sample data
#' set.seed(2203)
#' n <- 1000
#' true_params <- c(alpha = 2.0, beta = 1.5, gamma = 1.5, delta = 0.5)
#' data <- rbkw(n,
#'   alpha = true_params[1], beta = true_params[2],
#'   gamma = true_params[3], delta = true_params[4]
#' )
#'
#' # Evaluate gradient at true parameters
#' grad_true <- grbkw(par = true_params, data = data)
#' cat("Gradient at true parameters:\n")
#' print(grad_true)
#' cat("Norm:", sqrt(sum(grad_true^2)), "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 1.0, 1.0, 0.3),
#'   c(2.0, 1.5, 1.5, 0.5),
#'   c(2.5, 2.0, 2.0, 0.7)
#' )
#'
#' grad_norms <- apply(test_params, 1, function(p) {
#'   g <- grbkw(p, data)
#'   sqrt(sum(g^2))
#' })
#'
#' results <- data.frame(
#'   Alpha = test_params[, 1],
#'   Beta = test_params[, 2],
#'   Gamma = test_params[, 3],
#'   Delta = test_params[, 4],
#'   Grad_Norm = grad_norms
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Gradient in Optimization
#'
#' # Optimization with analytical gradient
#' fit_with_grad <- optim(
#'   par = c(1.8, 1.2, 1.1, 0.3),
#'   fn = llbkw,
#'   gr = grbkw,
#'   data = data,
#'   method = "Nelder-Mead",
#'   hessian = TRUE,
#'   control = list(trace = 0)
#' )
#'
#' # Optimization without gradient
#' fit_no_grad <- optim(
#'   par = c(1.8, 1.2, 1.1, 0.3),
#'   fn = llbkw,
#'   data = data,
#'   method = "Nelder-Mead",
#'   hessian = TRUE,
#'   control = list(trace = 0)
#' )
#'
#' comparison <- data.frame(
#'   Method = c("With Gradient", "Without Gradient"),
#'   Alpha = c(fit_with_grad$par[1], fit_no_grad$par[1]),
#'   Beta = c(fit_with_grad$par[2], fit_no_grad$par[2]),
#'   Gamma = c(fit_with_grad$par[3], fit_no_grad$par[3]),
#'   Delta = c(fit_with_grad$par[4], fit_no_grad$par[4]),
#'   NegLogLik = c(fit_with_grad$value, fit_no_grad$value),
#'   Iterations = c(fit_with_grad$counts[1], fit_no_grad$counts[1])
#' )
#' print(comparison, digits = 4, row.names = FALSE)
#'
#'
#' ## Example 3: Verifying Gradient at MLE
#'
#' mle <- fit_with_grad$par
#' names(mle) <- c("alpha", "beta", "gamma", "delta")
#'
#' # At MLE, gradient should be approximately zero
#' gradient_at_mle <- grbkw(par = mle, data = data)
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
#' # Compare at MLE
#' grad_analytical <- grbkw(par = mle, data = data)
#' grad_numerical <- numerical_gradient(llbkw, mle, data)
#'
#' comparison_grad <- data.frame(
#'   Parameter = c("alpha", "beta", "gamma", "delta"),
#'   Analytical = grad_analytical,
#'   Numerical = grad_numerical,
#'   Abs_Diff = abs(grad_analytical - grad_numerical),
#'   Rel_Error = abs(grad_analytical - grad_numerical) /
#'     (abs(grad_analytical) + 1e-10)
#' )
#' print(comparison_grad, digits = 8)
#'
#'
#' ## Example 5: Score Test Statistic
#'
#' # Score test for H0: theta = theta0
#' theta0 <- c(1.8, 1.3, 1.2, 0.4)
#' score_theta0 <- -grbkw(par = theta0, data = data)
#'
#' # Fisher information at theta0
#' fisher_info <- hsbkw(par = theta0, data = data)
#'
#' # Score test statistic
#' score_stat <- t(score_theta0) %*% solve(fisher_info) %*% score_theta0
#' p_value <- pchisq(score_stat, df = 4, lower.tail = FALSE)
#'
#' cat("\nScore Test:\n")
#' cat("H0: alpha=1.8, beta=1.3, gamma=1.2, delta=0.4\n")
#' cat("Test statistic:", score_stat, "\n")
#' cat("P-value:", format.pval(p_value, digits = 4), "\n")
#'
#'
#' ## Example 6: Confidence Ellipses (Selected pairs)
#'
#' # Observed information
#' obs_info <- hsbkw(par = mle, data = data)
#' vcov_full <- solve(obs_info)
#'
#' # Create confidence ellipses
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' # Alpha vs Beta ellipse
#' vcov_ab <- vcov_full[1:2, 1:2]
#' eig_decomp_ab <- eigen(vcov_ab)
#' ellipse_ab <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_ab[i, ] <- mle[1:2] + sqrt(chi2_val) *
#'     (eig_decomp_ab$vectors %*% diag(sqrt(eig_decomp_ab$values)) %*% v)
#' }
#'
#' # Alpha vs Gamma ellipse
#' vcov_ag <- vcov_full[c(1, 3), c(1, 3)]
#' eig_decomp_ag <- eigen(vcov_ag)
#' ellipse_ag <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_ag[i, ] <- mle[c(1, 3)] + sqrt(chi2_val) *
#'     (eig_decomp_ag$vectors %*% diag(sqrt(eig_decomp_ag$values)) %*% v)
#' }
#'
#' # Beta vs Delta ellipse
#' vcov_bd <- vcov_full[c(2, 4), c(2, 4)]
#' eig_decomp_bd <- eigen(vcov_bd)
#' ellipse_bd <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_bd[i, ] <- mle[c(2, 4)] + sqrt(chi2_val) *
#'     (eig_decomp_bd$vectors %*% diag(sqrt(eig_decomp_bd$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_ab <- sqrt(diag(vcov_ab))
#' ci_alpha_ab <- mle[1] + c(-1, 1) * 1.96 * se_ab[1]
#' ci_beta_ab <- mle[2] + c(-1, 1) * 1.96 * se_ab[2]
#'
#' se_ag <- sqrt(diag(vcov_ag))
#' ci_alpha_ag <- mle[1] + c(-1, 1) * 1.96 * se_ag[1]
#' ci_gamma_ag <- mle[3] + c(-1, 1) * 1.96 * se_ag[2]
#'
#' se_bd <- sqrt(diag(vcov_bd))
#' ci_beta_bd <- mle[2] + c(-1, 1) * 1.96 * se_bd[1]
#' ci_delta_bd <- mle[4] + c(-1, 1) * 1.96 * se_bd[2]
#'
#' # Plot selected ellipses
#'
#' # Alpha vs Beta
#' plot(ellipse_ab[, 1], ellipse_ab[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Alpha vs Beta", las = 1, xlim = range(ellipse_ab[, 1], ci_alpha_ab),
#'   ylim = range(ellipse_ab[, 2], ci_beta_ab)
#' )
#' abline(v = ci_alpha_ab, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_beta_ab, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Alpha vs Gamma
#' plot(ellipse_ag[, 1], ellipse_ag[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = expression(gamma),
#'   main = "Alpha vs Gamma", las = 1, xlim = range(ellipse_ag[, 1], ci_alpha_ag),
#'   ylim = range(ellipse_ag[, 2], ci_gamma_ag)
#' )
#' abline(v = ci_alpha_ag, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_gamma_ag, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Beta vs Delta
#' plot(ellipse_bd[, 1], ellipse_bd[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(beta), ylab = expression(delta),
#'   main = "Beta vs Delta", las = 1, xlim = range(ellipse_bd[, 1], ci_beta_bd),
#'   ylim = range(ellipse_bd[, 2], ci_delta_bd)
#' )
#' abline(v = ci_beta_bd, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_delta_bd, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[2], mle[4], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[4], pch = 17, col = "#006400", cex = 1.5)
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
grbkw <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 4) {
    stop("'par' must be a numeric vector of length 4")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_grbkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 7. HESSIAN (hsbkw)
# ----------------------------------------------------------------------------#

#' @title Hessian Matrix of the Negative Log-Likelihood for the BKw Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize hessian
#'
#' @description
#' Computes the analytic 4x4 Hessian matrix (matrix of second partial derivatives)
#' of the negative log-likelihood function for the Beta-Kumaraswamy (BKw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{gamma} (\eqn{\gamma}), and \code{delta} (\eqn{\delta}).
#' This distribution is the special case of the Generalized Kumaraswamy (GKw)
#' distribution where \eqn{\lambda = 1}. The Hessian is useful for estimating
#' standard errors and in optimization algorithms.
#'
#' @param par A numeric vector of length 4 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
#'   \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a 4x4 numeric matrix representing the Hessian matrix of the
#'   negative log-likelihood function, \eqn{-\partial^2 \ell / (\partial \theta_i \partial \theta_j)},
#'   where \eqn{\theta = (\alpha, \beta, \gamma, \delta)}.
#'   Returns a 4x4 matrix populated with \code{NaN} if any parameter values are
#'   invalid according to their constraints, or if any value in \code{data} is
#'   not in the interval (0, 1).
#'
#' @details
#' This function calculates the analytic second partial derivatives of the
#' negative log-likelihood function based on the BKw log-likelihood
#' (\eqn{\lambda=1} case of GKw, see \code{\link{llbkw}}):
#' \deqn{
#' \ell(\theta | \mathbf{x}) = n[\ln(\alpha) + \ln(\beta) - \ln B(\gamma, \delta+1)]
#' + \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta(\delta+1)-1)\ln(v_i) + (\gamma-1)\ln(w_i)]
#' }
#' where \eqn{\theta = (\alpha, \beta, \gamma, \delta)}, \eqn{B(a,b)}
#' is the Beta function (\code{\link[base]{beta}}), and intermediate terms are:
#' \itemize{
#'   \item \eqn{v_i = 1 - x_i^{\alpha}}
#'   \item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
#' }
#' The Hessian matrix returned contains the elements \eqn{- \frac{\partial^2 \ell(\theta | \mathbf{x})}{\partial \theta_i \partial \theta_j}}
#' for \eqn{\theta_i, \theta_j \in \{\alpha, \beta, \gamma, \delta\}}.
#'
#' Key properties of the returned matrix:
#' \itemize{
#'   \item Dimensions: 4x4.
#'   \item Symmetry: The matrix is symmetric.
#'   \item Ordering: Rows and columns correspond to the parameters in the order
#'     \eqn{\alpha, \beta, \gamma, \delta}.
#'   \item Content: Analytic second derivatives of the *negative* log-likelihood.
#' }
#' This corresponds to the relevant 4x4 submatrix of the 5x5 GKw Hessian (\code{\link{hsgkw}})
#' evaluated at \eqn{\lambda=1}. The exact analytical formulas are implemented directly.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' (Note: Specific Hessian formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{hsgkw}} (parent distribution Hessian),
#' \code{\link{llbkw}} (negative log-likelihood for BKw),
#' \code{\link{grbkw}} (gradient for BKw, if available),
#' \code{\link{dbkw}} (density for BKw),
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{hessian}} (for numerical Hessian comparison).
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Hessian Evaluation
#' # Generate sample data
#' set.seed(2203)
#' n <- 1000
#' true_params <- c(alpha = 2.0, beta = 1.5, gamma = 1.5, delta = 0.5)
#' data <- rbkw(n,
#'   alpha = true_params[1], beta = true_params[2],
#'   gamma = true_params[3], delta = true_params[4]
#' )
#'
#' # Evaluate Hessian at true parameters
#' hess_true <- hsbkw(par = true_params, data = data)
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
#'   par = c(1.8, 1.2, 1.1, 0.3),
#'   fn = llbkw,
#'   gr = grbkw,
#'   data = data,
#'   method = "Nelder-Mead",
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("alpha", "beta", "gamma", "delta")
#'
#' # Hessian at MLE
#' hessian_at_mle <- hsbkw(par = mle, data = data)
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
#' names(se) <- c("alpha", "beta", "gamma", "delta")
#'
#' # Correlation matrix
#' corr_matrix <- cov2cor(vcov_matrix)
#' cat("\nCorrelation Matrix:\n")
#' print(corr_matrix, digits = 4)
#'
#' # Confidence intervals
#' z_crit <- qnorm(0.975)
#' results <- data.frame(
#'   Parameter = c("alpha", "beta", "gamma", "delta"),
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
#'   c(1.5, 1.0, 1.0, 0.3),
#'   c(2.0, 1.5, 1.5, 0.5),
#'   mle,
#'   c(2.5, 2.0, 2.0, 0.7)
#' )
#'
#' hess_properties <- data.frame(
#'   Alpha = numeric(),
#'   Beta = numeric(),
#'   Gamma = numeric(),
#'   Delta = numeric(),
#'   Determinant = numeric(),
#'   Trace = numeric(),
#'   Min_Eigenval = numeric(),
#'   Max_Eigenval = numeric(),
#'   Cond_Number = numeric(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (i in 1:nrow(test_params)) {
#'   H <- hsbkw(par = test_params[i, ], data = data)
#'   eigs <- eigen(H, only.values = TRUE)$values
#'
#'   hess_properties <- rbind(hess_properties, data.frame(
#'     Alpha = test_params[i, 1],
#'     Beta = test_params[i, 2],
#'     Gamma = test_params[i, 3],
#'     Delta = test_params[i, 4],
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
#' ## Example 5: Curvature Visualization (Selected pairs)
#'
#' # Create grids around MLE with wider range (±1.5)
#' alpha_grid <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = 25)
#' beta_grid <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = 25)
#' gamma_grid <- seq(mle[3] - 1.5, mle[3] + 1.5, length.out = 25)
#' delta_grid <- seq(mle[4] - 1.5, mle[4] + 1.5, length.out = 25)
#'
#' alpha_grid <- alpha_grid[alpha_grid > 0]
#' beta_grid <- beta_grid[beta_grid > 0]
#' gamma_grid <- gamma_grid[gamma_grid > 0]
#' delta_grid <- delta_grid[delta_grid > 0]
#'
#' # Compute curvature measures for selected pairs
#' determinant_surface_ab <- matrix(NA, nrow = length(alpha_grid), ncol = length(beta_grid))
#' trace_surface_ab <- matrix(NA, nrow = length(alpha_grid), ncol = length(beta_grid))
#'
#' determinant_surface_ag <- matrix(NA, nrow = length(alpha_grid), ncol = length(gamma_grid))
#' trace_surface_ag <- matrix(NA, nrow = length(alpha_grid), ncol = length(gamma_grid))
#'
#' determinant_surface_bd <- matrix(NA, nrow = length(beta_grid), ncol = length(delta_grid))
#' trace_surface_bd <- matrix(NA, nrow = length(beta_grid), ncol = length(delta_grid))
#'
#' # Alpha vs Beta
#' for (i in seq_along(alpha_grid)) {
#'   for (j in seq_along(beta_grid)) {
#'     H <- hsbkw(c(alpha_grid[i], beta_grid[j], mle[3], mle[4]), data)
#'     determinant_surface_ab[i, j] <- det(H)
#'     trace_surface_ab[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Alpha vs Gamma
#' for (i in seq_along(alpha_grid)) {
#'   for (j in seq_along(gamma_grid)) {
#'     H <- hsbkw(c(alpha_grid[i], mle[2], gamma_grid[j], mle[4]), data)
#'     determinant_surface_ag[i, j] <- det(H)
#'     trace_surface_ag[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Beta vs Delta
#' for (i in seq_along(beta_grid)) {
#'   for (j in seq_along(delta_grid)) {
#'     H <- hsbkw(c(mle[1], beta_grid[i], mle[3], delta_grid[j]), data)
#'     determinant_surface_bd[i, j] <- det(H)
#'     trace_surface_bd[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Plot selected curvature surfaces
#'
#' # Determinant plots
#' contour(alpha_grid, beta_grid, determinant_surface_ab,
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Determinant: Alpha vs Beta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(alpha_grid, gamma_grid, determinant_surface_ag,
#'   xlab = expression(alpha), ylab = expression(gamma),
#'   main = "Determinant: Alpha vs Gamma", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(beta_grid, delta_grid, determinant_surface_bd,
#'   xlab = expression(beta), ylab = expression(delta),
#'   main = "Determinant: Beta vs Delta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[2], mle[4], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[4], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Trace plots
#' contour(alpha_grid, beta_grid, trace_surface_ab,
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Trace: Alpha vs Beta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(alpha_grid, gamma_grid, trace_surface_ag,
#'   xlab = expression(alpha), ylab = expression(gamma),
#'   main = "Trace: Alpha vs Gamma", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(beta_grid, delta_grid, trace_surface_bd,
#'   xlab = expression(beta), ylab = expression(delta),
#'   main = "Trace: Beta vs Delta", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[2], mle[4], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[4], pch = 17, col = "#006400", cex = 1.5)
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
#' ## Example 6: Confidence Ellipses (Selected pairs)
#'
#' # Extract selected 2x2 submatrices
#' vcov_ab <- vcov_matrix[1:2, 1:2]
#' vcov_ag <- vcov_matrix[c(1, 3), c(1, 3)]
#' vcov_bd <- vcov_matrix[c(2, 4), c(2, 4)]
#'
#' # Create confidence ellipses
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' # Alpha vs Beta ellipse
#' eig_decomp_ab <- eigen(vcov_ab)
#' ellipse_ab <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_ab[i, ] <- mle[1:2] + sqrt(chi2_val) *
#'     (eig_decomp_ab$vectors %*% diag(sqrt(eig_decomp_ab$values)) %*% v)
#' }
#'
#' # Alpha vs Gamma ellipse
#' eig_decomp_ag <- eigen(vcov_ag)
#' ellipse_ag <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_ag[i, ] <- mle[c(1, 3)] + sqrt(chi2_val) *
#'     (eig_decomp_ag$vectors %*% diag(sqrt(eig_decomp_ag$values)) %*% v)
#' }
#'
#' # Beta vs Delta ellipse
#' eig_decomp_bd <- eigen(vcov_bd)
#' ellipse_bd <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_bd[i, ] <- mle[c(2, 4)] + sqrt(chi2_val) *
#'     (eig_decomp_bd$vectors %*% diag(sqrt(eig_decomp_bd$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_ab <- sqrt(diag(vcov_ab))
#' ci_alpha_ab <- mle[1] + c(-1, 1) * 1.96 * se_ab[1]
#' ci_beta_ab <- mle[2] + c(-1, 1) * 1.96 * se_ab[2]
#'
#' se_ag <- sqrt(diag(vcov_ag))
#' ci_alpha_ag <- mle[1] + c(-1, 1) * 1.96 * se_ag[1]
#' ci_gamma_ag <- mle[3] + c(-1, 1) * 1.96 * se_ag[2]
#'
#' se_bd <- sqrt(diag(vcov_bd))
#' ci_beta_bd <- mle[2] + c(-1, 1) * 1.96 * se_bd[1]
#' ci_delta_bd <- mle[4] + c(-1, 1) * 1.96 * se_bd[2]
#'
#' # Plot selected ellipses side by side
#'
#' # Alpha vs Beta
#' plot(ellipse_ab[, 1], ellipse_ab[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "Alpha vs Beta", las = 1, xlim = range(ellipse_ab[, 1], ci_alpha_ab),
#'   ylim = range(ellipse_ab[, 2], ci_beta_ab)
#' )
#' abline(v = ci_alpha_ab, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_beta_ab, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Alpha vs Gamma
#' plot(ellipse_ag[, 1], ellipse_ag[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = expression(gamma),
#'   main = "Alpha vs Gamma", las = 1, xlim = range(ellipse_ag[, 1], ci_alpha_ag),
#'   ylim = range(ellipse_ag[, 2], ci_gamma_ag)
#' )
#' abline(v = ci_alpha_ag, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_gamma_ag, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[1], mle[3], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[3], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' # Beta vs Delta
#' plot(ellipse_bd[, 1], ellipse_bd[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(beta), ylab = expression(delta),
#'   main = "Beta vs Delta", las = 1, xlim = range(ellipse_bd[, 1], ci_beta_bd),
#'   ylim = range(ellipse_bd[, 2], ci_delta_bd)
#' )
#' abline(v = ci_beta_bd, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_delta_bd, col = "#808080", lty = 3, lwd = 1.5)
#' points(mle[2], mle[4], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[2], true_params[4], pch = 17, col = "#006400", cex = 1.5)
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
hsbkw <- function(par, data) {
  # Input validation
  if (!is.numeric(par) || length(par) != 4) {
    stop("'par' must be a numeric vector of length 4")
  }
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (length(data) < 1) {
    stop("'data' must have at least one observation")
  }

  # Call C++ implementation
  .Call("_gkwdist_hsbkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}
