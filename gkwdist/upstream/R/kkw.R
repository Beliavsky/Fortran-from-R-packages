# ============================================================================#
# KUMARASWAMY-KUMARASWAMY (KKw) DISTRIBUTION
# ============================================================================#
#
# Wrapper functions for the Kumaraswamy-Kumaraswamy distribution.
# C++ implementations are in src/kkw.cpp
#
# Functions:
#   - dkkw: Probability density function (PDF)
#   - pkkw: Cumulative distribution function (CDF)
#   - qkkw: Quantile function (inverse CDF)
#   - rkkw: Random number generation
#   - llkkw: Negative log-likelihood
#   - grkkw: Gradient of negative log-likelihood
#   - hskkw: Hessian of negative log-likelihood
# ============================================================================#


# ----------------------------------------------------------------------------#
# 1. DENSITY FUNCTION (dkkw)
# ----------------------------------------------------------------------------#

#' @title Density of the Kumaraswamy-Kumaraswamy (kkw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution density
#'
#' @description
#' Computes the probability density function (PDF) for the Kumaraswamy-Kumaraswamy
#' (kkw) distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}).
#' This distribution is defined on the interval (0, 1).
#'
#' @param x Vector of quantiles (values between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
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
#'   rule applied to the arguments (\code{x}, \code{alpha}, \code{beta},
#'   \code{delta}, \code{lambda}). Returns \code{0} (or \code{-Inf} if
#'   \code{log = TRUE}) for \code{x} outside the interval (0, 1), or
#'   \code{NaN} if parameters are invalid (e.g., \code{alpha <= 0}, \code{beta <= 0},
#'   \code{delta < 0}, \code{lambda <= 0}).
#'
#' @details
#' The Kumaraswamy-Kumaraswamy (kkw) distribution is a special case of the
#' five-parameter Generalized Kumaraswamy distribution (\code{\link{dgkw}})
#' obtained by setting the parameter \eqn{\gamma = 1}.
#'
#' The probability density function is given by:
#' \deqn{
#' f(x; \alpha, \beta, \delta, \lambda) = (\delta + 1) \lambda \alpha \beta x^{\alpha - 1} (1 - x^\alpha)^{\beta - 1} \bigl[1 - (1 - x^\alpha)^\beta\bigr]^{\lambda - 1} \bigl\{1 - \bigl[1 - (1 - x^\alpha)^\beta\bigr]^\lambda\bigr\}^{\delta}
#' }
#' for \eqn{0 < x < 1}. Note that \eqn{1/(\delta+1)} corresponds to the Beta function
#' term \eqn{B(1, \delta+1)} when \eqn{\gamma=1}.
#'
#' Numerical evaluation follows similar stability considerations as \code{\link{dgkw}}.
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
#' \code{\link{pkkw}}, \code{\link{qkkw}}, \code{\link{rkkw}} (if they exist),
#' \code{\link[stats]{dbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' x_vals <- c(0.2, 0.5, 0.8)
#' alpha_par <- 2.0
#' beta_par <- 3.0
#' delta_par <- 0.5
#' lambda_par <- 1.5
#'
#' # Calculate density
#' densities <- dkkw(x_vals, alpha_par, beta_par, delta_par, lambda_par)
#' print(densities)
#'
#' # Calculate log-density
#' log_densities <- dkkw(x_vals, alpha_par, beta_par, delta_par, lambda_par,
#'   log = TRUE
#' )
#' print(log_densities)
#' # Check: should match log(densities)
#' print(log(densities))
#'
#' # Compare with dgkw setting gamma = 1
#' densities_gkw <- dgkw(x_vals, alpha_par, beta_par,
#'   gamma = 1.0,
#'   delta_par, lambda_par
#' )
#' print(paste("Max difference:", max(abs(densities - densities_gkw)))) # Should be near zero
#'
#' # Plot the density
#' curve_x <- seq(0.01, 0.99, length.out = 200)
#' curve_y <- dkkw(curve_x, alpha_par, beta_par, delta_par, lambda_par)
#' plot(curve_x, curve_y,
#'   type = "l", main = "kkw Density Example",
#'   xlab = "x", ylab = "f(x)", col = "blue"
#' )
#' }
#'
#' @export
dkkw <- function(x, alpha = 1, beta = 1, delta = 0, lambda = 1, log = FALSE) {
  # Input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
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
  .Call("_gkwdist_dkkw",
    as.numeric(x),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(delta),
    as.numeric(lambda),
    as.logical(log),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 2. DISTRIBUTION FUNCTION (pkkw)
# ----------------------------------------------------------------------------#

#' @title Cumulative Distribution Function (CDF) of the kkw Distribution
#' @author Lopes, J. E.
#' @keywords distribution cumulative
#'
#' @description
#' Computes the cumulative distribution function (CDF), \eqn{P(X \le q)}, for the
#' Kumaraswamy-Kumaraswamy (kkw) distribution with parameters \code{alpha}
#' (\eqn{\alpha}), \code{beta} (\eqn{\beta}), \code{delta} (\eqn{\delta}),
#' and \code{lambda} (\eqn{\lambda}). This distribution is defined on the
#' interval (0, 1).
#'
#' @param q Vector of quantiles (values generally between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
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
#'   \code{alpha}, \code{beta}, \code{delta}, \code{lambda}). Returns \code{0}
#'   (or \code{-Inf} if \code{log.p = TRUE}) for \code{q <= 0} and \code{1}
#'   (or \code{0} if \code{log.p = TRUE}) for \code{q >= 1}. Returns \code{NaN}
#'   for invalid parameters.
#'
#' @details
#' The Kumaraswamy-Kumaraswamy (kkw) distribution is a special case of the
#' five-parameter Generalized Kumaraswamy distribution (\code{\link{pgkw}})
#' obtained by setting the shape parameter \eqn{\gamma = 1}.
#'
#' The CDF of the GKw distribution is \eqn{F_{GKw}(q) = I_{y(q)}(\gamma, \delta+1)},
#' where \eqn{y(q) = [1-(1-q^{\alpha})^{\beta}]^{\lambda}} and \eqn{I_x(a,b)}
#' is the regularized incomplete beta function (\code{\link[stats]{pbeta}}).
#' Setting \eqn{\gamma=1} utilizes the property \eqn{I_x(1, b) = 1 - (1-x)^b},
#' yielding the kkw CDF:
#' \deqn{
#' F(q; \alpha, \beta, \delta, \lambda) = 1 - \bigl\{1 - \bigl[1 - (1 - q^\alpha)^\beta\bigr]^\lambda\bigr\}^{\delta + 1}
#' }
#' for \eqn{0 < q < 1}.
#'
#' The implementation uses this closed-form expression for efficiency and handles
#' \code{lower.tail} and \code{log.p} arguments appropriately.
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
#' \code{\link{dkkw}}, \code{\link{qkkw}}, \code{\link{rkkw}},
#' \code{\link[stats]{pbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' q_vals <- c(0.2, 0.5, 0.8)
#' alpha_par <- 2.0
#' beta_par <- 3.0
#' delta_par <- 0.5
#' lambda_par <- 1.5
#'
#' # Calculate CDF P(X <= q)
#' probs <- pkkw(q_vals, alpha_par, beta_par, delta_par, lambda_par)
#' print(probs)
#'
#' # Calculate upper tail P(X > q)
#' probs_upper <- pkkw(q_vals, alpha_par, beta_par, delta_par, lambda_par,
#'   lower.tail = FALSE
#' )
#' print(probs_upper)
#' # Check: probs + probs_upper should be 1
#' print(probs + probs_upper)
#'
#' # Calculate log CDF
#' logs <- pkkw(q_vals, alpha_par, beta_par, delta_par, lambda_par,
#'   log.p = TRUE
#' )
#' print(logs)
#' # Check: should match log(probs)
#' print(log(probs))
#'
#' # Compare with pgkw setting gamma = 1
#' probs_gkw <- pgkw(q_vals, alpha_par, beta_par,
#'   gamma = 1.0,
#'   delta_par, lambda_par
#' )
#' print(paste("Max difference:", max(abs(probs - probs_gkw)))) # Should be near zero
#'
#' # Plot the CDF
#' curve_q <- seq(0.01, 0.99, length.out = 200)
#' curve_p <- pkkw(curve_q, alpha_par, beta_par, delta_par, lambda_par)
#' plot(curve_q, curve_p,
#'   type = "l", main = "kkw CDF Example",
#'   xlab = "q", ylab = "F(q)", col = "blue", ylim = c(0, 1)
#' )
#' }
#'
#' @export
pkkw <- function(q, alpha = 1, beta = 1, delta = 0, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(q)) stop("'q' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
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
  .Call("_gkwdist_pkkw",
    as.numeric(q),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(delta),
    as.numeric(lambda),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 3. QUANTILE FUNCTION (qkkw)
# ----------------------------------------------------------------------------#

#' @title Quantile Function of the Kumaraswamy-Kumaraswamy (kkw) Distribution
#' @author Lopes, J. E.
#' @keywords distribution quantile
#'
#' @description
#' Computes the quantile function (inverse CDF) for the Kumaraswamy-Kumaraswamy
#' (kkw) distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}).
#' It finds the value \code{q} such that \eqn{P(X \le q) = p}. This distribution
#' is a special case of the Generalized Kumaraswamy (GKw) distribution where
#' the parameter \eqn{\gamma = 1}.
#'
#' @param p Vector of probabilities (values between 0 and 1).
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
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
#'   the arguments (\code{p}, \code{alpha}, \code{beta}, \code{delta},
#'   \code{lambda}). Returns:
#'   \itemize{
#'     \item \code{0} for \code{p = 0} (or \code{p = -Inf} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{1} for \code{p = 1} (or \code{p = 0} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{NaN} for \code{p < 0} or \code{p > 1} (or corresponding log scale).
#'     \item \code{NaN} for invalid parameters (e.g., \code{alpha <= 0},
#'           \code{beta <= 0}, \code{delta < 0}, \code{lambda <= 0}).
#'   }
#'   Boundary return values are adjusted accordingly for \code{lower.tail = FALSE}.
#'
#' @details
#' The quantile function \eqn{Q(p)} is the inverse of the CDF \eqn{F(q)}. The CDF
#' for the kkw (\eqn{\gamma=1}) distribution is (see \code{\link{pkkw}}):
#' \deqn{
#' F(q) = 1 - \bigl\{1 - \bigl[1 - (1 - q^\alpha)^\beta\bigr]^\lambda\bigr\}^{\delta + 1}
#' }
#' Inverting this equation for \eqn{q} yields the quantile function:
#' \deqn{
#' Q(p) = \left[ 1 - \left\{ 1 - \left[ 1 - (1 - p)^{1/(\delta+1)} \right]^{1/\lambda} \right\}^{1/\beta} \right]^{1/\alpha}
#' }
#' The function uses this closed-form expression and correctly handles the
#' \code{lower.tail} and \code{log.p} arguments by transforming \code{p}
#' appropriately before applying the formula.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' @seealso
#' \code{\link{qgkw}} (parent distribution quantile function),
#' \code{\link{dkkw}}, \code{\link{pkkw}}, \code{\link{rkkw}},
#' \code{\link[stats]{qbeta}}
#'
#' @examples
#' \donttest{
#' # Example values
#' p_vals <- c(0.1, 0.5, 0.9)
#' alpha_par <- 2.0
#' beta_par <- 3.0
#' delta_par <- 0.5
#' lambda_par <- 1.5
#'
#' # Calculate quantiles
#' quantiles <- qkkw(p_vals, alpha_par, beta_par, delta_par, lambda_par)
#' print(quantiles)
#'
#' # Calculate quantiles for upper tail probabilities P(X > q) = p
#' # e.g., for p=0.1, find q such that P(X > q) = 0.1 (90th percentile)
#' quantiles_upper <- qkkw(p_vals, alpha_par, beta_par, delta_par, lambda_par,
#'   lower.tail = FALSE
#' )
#' print(quantiles_upper)
#' # Check: qkkw(p, ..., lt=F) == qkkw(1-p, ..., lt=T)
#' print(qkkw(1 - p_vals, alpha_par, beta_par, delta_par, lambda_par))
#'
#' # Calculate quantiles from log probabilities
#' log.p_vals <- log(p_vals)
#' quantiles_logp <- qkkw(log.p_vals, alpha_par, beta_par, delta_par, lambda_par,
#'   log.p = TRUE
#' )
#' print(quantiles_logp)
#' # Check: should match original quantiles
#' print(quantiles)
#'
#' # Compare with qgkw setting gamma = 1
#' quantiles_gkw <- qgkw(p_vals, alpha_par, beta_par,
#'   gamma = 1.0,
#'   delta_par, lambda_par
#' )
#' print(paste("Max difference:", max(abs(quantiles - quantiles_gkw)))) # Should be near zero
#'
#' # Verify inverse relationship with pkkw
#' p_check <- 0.75
#' q_calc <- qkkw(p_check, alpha_par, beta_par, delta_par, lambda_par)
#' p_recalc <- pkkw(q_calc, alpha_par, beta_par, delta_par, lambda_par)
#' print(paste("Original p:", p_check, " Recalculated p:", p_recalc))
#' # abs(p_check - p_recalc) < 1e-9 # Should be TRUE
#'
#' # Boundary conditions
#' print(qkkw(c(0, 1), alpha_par, beta_par, delta_par, lambda_par)) # Should be 0, 1
#' print(qkkw(c(-Inf, 0), alpha_par, beta_par, delta_par, lambda_par, log.p = TRUE)) # Should be 0, 1
#' }
#'
#' @export
qkkw <- function(p, alpha = 1, beta = 1, delta = 0, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(p)) stop("'p' must be numeric")
  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive")
  }
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("'beta' must be positive")
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
  .Call("_gkwdist_qkkw",
    as.numeric(p),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(delta),
    as.numeric(lambda),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 4. RANDOM GENERATION (rkkw)
# ----------------------------------------------------------------------------#

#' @title Random Number Generation for the kkw Distribution
#' @author Lopes, J. E.
#' @keywords distribution random
#'
#' @description
#' Generates random deviates from the Kumaraswamy-Kumaraswamy (kkw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}).
#' This distribution is a special case of the Generalized Kumaraswamy (GKw)
#' distribution where the parameter \eqn{\gamma = 1}.
#'
#' @param n Number of observations. If \code{length(n) > 1}, the length is
#'   taken to be the number required. Must be a non-negative integer.
#' @param alpha Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param beta Shape parameter \code{beta} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#' @param delta Shape parameter \code{delta} >= 0. Can be a scalar or a vector.
#'   Default: 0.0.
#' @param lambda Shape parameter \code{lambda} > 0. Can be a scalar or a vector.
#'   Default: 1.0.
#'
#' @return A vector of length \code{n} containing random deviates from the kkw
#'   distribution. The length of the result is determined by \code{n} and the
#'   recycling rule applied to the parameters (\code{alpha}, \code{beta},
#'   \code{delta}, \code{lambda}). Returns \code{NaN} if parameters
#'   are invalid (e.g., \code{alpha <= 0}, \code{beta <= 0}, \code{delta < 0},
#'   \code{lambda <= 0}).
#'
#' @details
#' The generation method uses the inverse transform method based on the quantile
#' function (\code{\link{qkkw}}). The kkw quantile function is:
#' \deqn{
#' Q(p) = \left[ 1 - \left\{ 1 - \left[ 1 - (1 - p)^{1/(\delta+1)} \right]^{1/\lambda} \right\}^{1/\beta} \right]^{1/\alpha}
#' }
#' Random deviates are generated by evaluating \eqn{Q(p)} where \eqn{p} is a
#' random variable following the standard Uniform distribution on (0, 1)
#' (\code{\link[stats]{runif}}).
#'
#' This is equivalent to the general method for the GKw distribution
#' (\code{\link{rgkw}}) specialized for \eqn{\gamma=1}. The GKw method generates
#' \eqn{W \sim \mathrm{Beta}(\gamma, \delta+1)} and then applies transformations.
#' When \eqn{\gamma=1}, \eqn{W \sim \mathrm{Beta}(1, \delta+1)}, which can be
#' generated via \eqn{W = 1 - V^{1/(\delta+1)}} where \eqn{V \sim \mathrm{Unif}(0,1)}.
#' Substituting this \eqn{W} into the GKw transformation yields the same result
#' as evaluating \eqn{Q(1-V)} above (noting \eqn{p = 1-V} is also Uniform).
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer-Verlag.
#' (General methods for random variate generation).
#'
#' @seealso
#' \code{\link{rgkw}} (parent distribution random generation),
#' \code{\link{dkkw}}, \code{\link{pkkw}}, \code{\link{qkkw}},
#' \code{\link[stats]{runif}}, \code{\link[stats]{rbeta}}
#'
#' @examples
#' \donttest{
#' set.seed(2025) # for reproducibility
#'
#' # Generate 1000 random values from a specific kkw distribution
#' alpha_par <- 2.0
#' beta_par <- 3.0
#' delta_par <- 0.5
#' lambda_par <- 1.5
#'
#' x_sample_kkw <- rkkw(1000,
#'   alpha = alpha_par, beta = beta_par,
#'   delta = delta_par, lambda = lambda_par
#' )
#' summary(x_sample_kkw)
#'
#' # Histogram of generated values compared to theoretical density
#' hist(x_sample_kkw,
#'   breaks = 30, freq = FALSE, # freq=FALSE for density
#'   main = "Histogram of kkw Sample", xlab = "x", ylim = c(0, 3.5)
#' )
#' curve(
#'   dkkw(x,
#'     alpha = alpha_par, beta = beta_par, delta = delta_par,
#'     lambda = lambda_par
#'   ),
#'   add = TRUE, col = "red", lwd = 2, n = 201
#' )
#' legend("topright", legend = "Theoretical PDF", col = "red", lwd = 2, bty = "n")
#'
#' # Comparing empirical and theoretical quantiles (Q-Q plot)
#' prob_points <- seq(0.01, 0.99, by = 0.01)
#' theo_quantiles <- qkkw(prob_points,
#'   alpha = alpha_par, beta = beta_par,
#'   delta = delta_par, lambda = lambda_par
#' )
#' emp_quantiles <- quantile(x_sample_kkw, prob_points, type = 7) # type 7 is default
#'
#' plot(theo_quantiles, emp_quantiles,
#'   pch = 16, cex = 0.8,
#'   main = "Q-Q Plot for kkw Distribution",
#'   xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles (n=1000)"
#' )
#' abline(a = 0, b = 1, col = "blue", lty = 2)
#'
#' # Compare summary stats with rgkw(..., gamma=1, ...)
#' # Note: individual values will differ due to randomness
#' x_sample_gkw <- rgkw(1000,
#'   alpha = alpha_par, beta = beta_par, gamma = 1.0,
#'   delta = delta_par, lambda = lambda_par
#' )
#' print("Summary stats for rkkw sample:")
#' print(summary(x_sample_kkw))
#' print("Summary stats for rgkw(gamma=1) sample:")
#' print(summary(x_sample_gkw)) # Should be similar
#' }
#'
#' @export
rkkw <- function(n, alpha = 1, beta = 1, delta = 0, lambda = 1) {
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
  if (!is.numeric(delta) || any(delta < 0)) {
    stop("'delta' must be non-negative")
  }
  if (!is.numeric(lambda) || any(lambda <= 0)) {
    stop("'lambda' must be positive")
  }

  # Call C++ implementation
  .Call("_gkwdist_rkkw",
    as.integer(n),
    as.numeric(alpha),
    as.numeric(beta),
    as.numeric(delta),
    as.numeric(lambda),
    PACKAGE = "gkwdist"
  )
}


# ============================================================================#
# MAXIMUM LIKELIHOOD ESTIMATION FUNCTIONS
# ============================================================================#

# ----------------------------------------------------------------------------#
# 5. NEGATIVE LOG-LIKELIHOOD (llkkw)
# ----------------------------------------------------------------------------#

#' @title Negative Log-Likelihood for the kkw Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize
#'
#' @description
#' Computes the negative log-likelihood function for the Kumaraswamy-Kumaraswamy
#' (kkw) distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}),
#' given a vector of observations. This distribution is a special case of the
#' Generalized Kumaraswamy (GKw) distribution where \eqn{\gamma = 1}.
#'
#' @param par A numeric vector of length 4 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
#'   \code{delta} (\eqn{\delta \ge 0}), \code{lambda} (\eqn{\lambda > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a single \code{double} value representing the negative
#'   log-likelihood (\eqn{-\ell(\theta|\mathbf{x})}). Returns \code{Inf}
#'   if any parameter values in \code{par} are invalid according to their
#'   constraints, or if any value in \code{data} is not in the interval (0, 1).
#'
#' @details
#' The kkw distribution is the GKw distribution (\code{\link{dgkw}}) with \eqn{\gamma=1}.
#' Its probability density function (PDF) is:
#' \deqn{
#' f(x | \theta) = (\delta + 1) \lambda \alpha \beta x^{\alpha - 1} (1 - x^\alpha)^{\beta - 1} \bigl[1 - (1 - x^\alpha)^\beta\bigr]^{\lambda - 1} \bigl\{1 - \bigl[1 - (1 - x^\alpha)^\beta\bigr]^\lambda\bigr\}^{\delta}
#' }
#' for \eqn{0 < x < 1} and \eqn{\theta = (\alpha, \beta, \delta, \lambda)}.
#' The log-likelihood function \eqn{\ell(\theta | \mathbf{x})} for a sample
#' \eqn{\mathbf{x} = (x_1, \dots, x_n)} is \eqn{\sum_{i=1}^n \ln f(x_i | \theta)}:
#' \deqn{
#' \ell(\theta | \mathbf{x}) = n[\ln(\delta+1) + \ln(\lambda) + \ln(\alpha) + \ln(\beta)]
#' + \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta-1)\ln(v_i) + (\lambda-1)\ln(w_i) + \delta\ln(z_i)]
#' }
#' where:
#' \itemize{
#'   \item \eqn{v_i = 1 - x_i^{\alpha}}
#'   \item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
#'   \item \eqn{z_i = 1 - w_i^{\lambda} = 1 - [1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
#' }
#' This function computes and returns the *negative* log-likelihood, \eqn{-\ell(\theta|\mathbf{x})},
#' suitable for minimization using optimization routines like \code{\link[stats]{optim}}.
#' Numerical stability is maintained similarly to \code{\link{llgkw}}.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' @seealso
#' \code{\link{llgkw}} (parent distribution negative log-likelihood),
#' \code{\link{dkkw}}, \code{\link{pkkw}}, \code{\link{qkkw}}, \code{\link{rkkw}},
#' \code{\link{grkkw}} (gradient, if available),
#' \code{\link{hskkw}} (Hessian, if available),
#' \code{\link[stats]{optim}}
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Log-Likelihood Evaluation
#'
#' # Generate sample data
#' set.seed(123)
#' n <- 1000
#' true_params <- c(alpha = 2.0, beta = 3.0, delta = 1.5, lambda = 2.0)
#' data <- rkkw(n,
#'   alpha = true_params[1], beta = true_params[2],
#'   delta = true_params[3], lambda = true_params[4]
#' )
#'
#' # Evaluate negative log-likelihood at true parameters
#' nll_true <- llkkw(par = true_params, data = data)
#' cat("Negative log-likelihood at true parameters:", nll_true, "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 2.5, 1.0, 1.5),
#'   c(2.0, 3.0, 1.5, 2.0),
#'   c(2.5, 3.5, 2.0, 2.5)
#' )
#'
#' nll_values <- apply(test_params, 1, function(p) llkkw(p, data))
#' results <- data.frame(
#'   Alpha = test_params[, 1],
#'   Beta = test_params[, 2],
#'   Delta = test_params[, 3],
#'   Lambda = test_params[, 4],
#'   NegLogLik = nll_values
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Maximum Likelihood Estimation
#'
#' # Optimization using BFGS with analytical gradient
#' fit <- optim(
#'   par = c(1.5, 2.5, 1.0, 1.5),
#'   fn = llkkw,
#'   gr = grkkw,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("alpha", "beta", "delta", "lambda")
#' se <- sqrt(diag(solve(fit$hessian)))
#'
#' results <- data.frame(
#'   Parameter = c("alpha", "beta", "delta", "lambda"),
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
#' start_params <- c(1.5, 2.5, 1.0, 1.5)
#'
#' comparison <- data.frame(
#'   Method = character(),
#'   Alpha = numeric(),
#'   Beta = numeric(),
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
#'       fn = llkkw,
#'       gr = grkkw,
#'       data = data,
#'       method = method
#'     )
#'   } else if (method == "L-BFGS-B") {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llkkw,
#'       gr = grkkw,
#'       data = data,
#'       method = method,
#'       lower = c(0.01, 0.01, 0.01, 0.01),
#'       upper = c(100, 100, 100, 100)
#'     )
#'   } else {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llkkw,
#'       data = data,
#'       method = method
#'     )
#'   }
#'
#'   comparison <- rbind(comparison, data.frame(
#'     Method = method,
#'     Alpha = fit_temp$par[1],
#'     Beta = fit_temp$par[2],
#'     Delta = fit_temp$par[3],
#'     Lambda = fit_temp$par[4],
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
#' # Test H0: delta = 1.5 vs H1: delta free
#' loglik_full <- -fit$value
#'
#' restricted_ll <- function(params_restricted, data, delta_fixed) {
#'   llkkw(par = c(
#'     params_restricted[1], params_restricted[2],
#'     delta_fixed, params_restricted[3]
#'   ), data = data)
#' }
#'
#' fit_restricted <- optim(
#'   par = c(mle[1], mle[2], mle[4]),
#'   fn = restricted_ll,
#'   data = data,
#'   delta_fixed = 1.5,
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
#' # Profile for alpha
#' alpha_grid <- seq(mle[1] - 1, mle[1] + 1, length.out = 40)
#' alpha_grid <- alpha_grid[alpha_grid > 0]
#' profile_ll_alpha <- numeric(length(alpha_grid))
#'
#' for (i in seq_along(alpha_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-1],
#'     fn = function(p) llkkw(c(alpha_grid[i], p), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_alpha[i] <- -profile_fit$value
#' }
#'
#' # Profile for beta
#' beta_grid <- seq(mle[2] - 1, mle[2] + 1, length.out = 40)
#' beta_grid <- beta_grid[beta_grid > 0]
#' profile_ll_beta <- numeric(length(beta_grid))
#'
#' for (i in seq_along(beta_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-2],
#'     fn = function(p) llkkw(c(p[1], beta_grid[i], p[2], p[3]), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_beta[i] <- -profile_fit$value
#' }
#'
#' # Profile for delta
#' delta_grid <- seq(mle[3] - 0.8, mle[3] + 0.8, length.out = 40)
#' delta_grid <- delta_grid[delta_grid > 0]
#' profile_ll_delta <- numeric(length(delta_grid))
#'
#' for (i in seq_along(delta_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-3],
#'     fn = function(p) llkkw(c(p[1], p[2], delta_grid[i], p[3]), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_delta[i] <- -profile_fit$value
#' }
#'
#' # Profile for lambda
#' lambda_grid <- seq(mle[4] - 1, mle[4] + 1, length.out = 40)
#' lambda_grid <- lambda_grid[lambda_grid > 0]
#' profile_ll_lambda <- numeric(length(lambda_grid))
#'
#' for (i in seq_along(lambda_grid)) {
#'   profile_fit <- optim(
#'     par = mle[-4],
#'     fn = function(p) llkkw(c(p[1], p[2], p[3], lambda_grid[i]), data),
#'     method = "Nelder-Mead"
#'   )
#'   profile_ll_lambda[i] <- -profile_fit$value
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
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.7
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
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.7
#' )
#' grid(col = "gray90")
#'
#' plot(delta_grid, profile_ll_delta,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(delta), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", delta)), las = 1
#' )
#' abline(v = mle[3], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[3], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3, lwd = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#808080"),
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.7
#' )
#' grid(col = "gray90")
#'
#' plot(lambda_grid, profile_ll_lambda,
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(lambda), ylab = "Profile Log-Likelihood",
#'   main = expression(paste("Profile: ", lambda)), las = 1
#' )
#' abline(v = mle[4], col = "#8B0000", lty = 2, lwd = 2)
#' abline(v = true_params[4], col = "#006400", lty = 2, lwd = 2)
#' abline(h = threshold, col = "#808080", lty = 3, lwd = 1.5)
#' legend("topright",
#'   legend = c("MLE", "True", "95% CI"),
#'   col = c("#8B0000", "#006400", "#808080"),
#'   lty = c(2, 2, 3), lwd = 2, bty = "n", cex = 0.7
#' )
#' grid(col = "gray90")
#'
#' ## Example 6: 2D Log-Likelihood Surface (Alpha vs Beta)
#'
#' # Create 2D grid
#' alpha_2d <- seq(mle[1] - 0.8, mle[1] + 0.8, length.out = round(n / 25))
#' beta_2d <- seq(mle[2] - 0.8, mle[2] + 0.8, length.out = round(n / 25))
#' alpha_2d <- alpha_2d[alpha_2d > 0]
#' beta_2d <- beta_2d[beta_2d > 0]
#'
#' # Compute log-likelihood surface
#' ll_surface <- matrix(NA, nrow = length(alpha_2d), ncol = length(beta_2d))
#'
#' for (i in seq_along(alpha_2d)) {
#'   for (j in seq_along(beta_2d)) {
#'     ll_surface[i, j] <- -llkkw(c(alpha_2d[i], beta_2d[j], mle[3], mle[4]), data)
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
#'   main = "2D Log-Likelihood: Alpha vs Beta",
#'   levels = seq(min(ll_surface, na.rm = TRUE), max_ll, length.out = 20),
#'   col = "#2E4057", las = 1, lwd = 1
#' )
#'
#' contour(alpha_2d, beta_2d, ll_surface,
#'   levels = c(levels_90, levels_95, levels_99),
#'   col = c("#FFA07A", "#FF6347", "#8B0000"),
#'   lwd = c(2, 2.5, 3), lty = c(3, 2, 1),
#'   add = TRUE, labcex = 0.8
#' )
#'
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
#'
#' ## Example 7: 2D Log-Likelihood Surface (Delta vs Lambda)
#'
#' # Create 2D grid
#' delta_2d <- seq(mle[3] - 0.6, mle[3] + 0.6, length.out = round(n / 25))
#' lambda_2d <- seq(mle[4] - 0.8, mle[4] + 0.8, length.out = round(n / 25))
#' delta_2d <- delta_2d[delta_2d > 0]
#' lambda_2d <- lambda_2d[lambda_2d > 0]
#'
#' # Compute log-likelihood surface
#' ll_surface2 <- matrix(NA, nrow = length(delta_2d), ncol = length(lambda_2d))
#'
#' for (i in seq_along(delta_2d)) {
#'   for (j in seq_along(lambda_2d)) {
#'     ll_surface2[i, j] <- -llkkw(c(mle[1], mle[2], delta_2d[i], lambda_2d[j]), data)
#'   }
#' }
#'
#' # Confidence region levels
#' max_ll2 <- max(ll_surface2, na.rm = TRUE)
#' levels2_90 <- max_ll2 - qchisq(0.90, df = 2) / 2
#' levels2_95 <- max_ll2 - qchisq(0.95, df = 2) / 2
#' levels2_99 <- max_ll2 - qchisq(0.99, df = 2) / 2
#'
#' # Plot contour
#' contour(delta_2d, lambda_2d, ll_surface2,
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "2D Log-Likelihood: Delta vs Lambda",
#'   levels = seq(min(ll_surface2, na.rm = TRUE), max_ll2, length.out = 20),
#'   col = "#2E4057", las = 1, lwd = 1
#' )
#'
#' contour(delta_2d, lambda_2d, ll_surface2,
#'   levels = c(levels2_90, levels2_95, levels2_99),
#'   col = c("#FFA07A", "#FF6347", "#8B0000"),
#'   lwd = c(2, 2.5, 3), lty = c(3, 2, 1),
#'   add = TRUE, labcex = 0.8
#' )
#'
#' points(mle[3], mle[4], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[3], true_params[4], pch = 17, col = "#006400", cex = 1.5)
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
#' }
#'
#' @export
llkkw <- function(par, data) {
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
  .Call("_gkwdist_llkkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 6. GRADIENT (grkkw)
# ----------------------------------------------------------------------------#

#' @title Gradient of the Negative Log-Likelihood for the kkw Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize gradient
#'
#' @description
#' Computes the gradient vector (vector of first partial derivatives) of the
#' negative log-likelihood function for the Kumaraswamy-Kumaraswamy (kkw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}).
#' This distribution is the special case of the Generalized Kumaraswamy (GKw)
#' distribution where \eqn{\gamma = 1}. The gradient is typically used in
#' optimization algorithms for maximum likelihood estimation.
#'
#' @param par A numeric vector of length 4 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
#'   \code{delta} (\eqn{\delta \ge 0}), \code{lambda} (\eqn{\lambda > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a numeric vector of length 4 containing the partial derivatives
#'   of the negative log-likelihood function \eqn{-\ell(\theta | \mathbf{x})} with
#'   respect to each parameter:
#'   \eqn{(-\partial \ell/\partial \alpha, -\partial \ell/\partial \beta, -\partial \ell/\partial \delta, -\partial \ell/\partial \lambda)}.
#'   Returns a vector of \code{NaN} if any parameter values are invalid according
#'   to their constraints, or if any value in \code{data} is not in the
#'   interval (0, 1).
#'
#' @details
#' The components of the gradient vector of the negative log-likelihood
#' (\eqn{-\nabla \ell(\theta | \mathbf{x})}) for the kkw (\eqn{\gamma=1}) model are:
#'
#' \deqn{
#' -\frac{\partial \ell}{\partial \alpha} = -\frac{n}{\alpha} - \sum_{i=1}^{n}\ln(x_i)
#' + (\beta-1)\sum_{i=1}^{n}\frac{x_i^{\alpha}\ln(x_i)}{v_i}
#' - (\lambda-1)\sum_{i=1}^{n}\frac{\beta v_i^{\beta-1} x_i^{\alpha}\ln(x_i)}{w_i}
#' + \delta\sum_{i=1}^{n}\frac{\lambda w_i^{\lambda-1} \beta v_i^{\beta-1} x_i^{\alpha}\ln(x_i)}{z_i}
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \beta} = -\frac{n}{\beta} - \sum_{i=1}^{n}\ln(v_i)
#' + (\lambda-1)\sum_{i=1}^{n}\frac{v_i^{\beta}\ln(v_i)}{w_i}
#' - \delta\sum_{i=1}^{n}\frac{\lambda w_i^{\lambda-1} v_i^{\beta}\ln(v_i)}{z_i}
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \delta} = -\frac{n}{\delta+1} - \sum_{i=1}^{n}\ln(z_i)
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \lambda} = -\frac{n}{\lambda} - \sum_{i=1}^{n}\ln(w_i)
#' + \delta\sum_{i=1}^{n}\frac{w_i^{\lambda}\ln(w_i)}{z_i}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{v_i = 1 - x_i^{\alpha}}
#'   \item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
#'   \item \eqn{z_i = 1 - w_i^{\lambda} = 1 - [1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
#' }
#' These formulas represent the derivatives of \eqn{-\ell(\theta)}, consistent with
#' minimizing the negative log-likelihood. They correspond to the general GKw
#' gradient (\code{\link{grgkw}}) components for \eqn{\alpha, \beta, \delta, \lambda}
#' evaluated at \eqn{\gamma=1}. Note that the component for \eqn{\gamma} is omitted.
#' Numerical stability is maintained through careful implementation.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#' @seealso
#' \code{\link{grgkw}} (parent distribution gradient),
#' \code{\link{llkkw}} (negative log-likelihood for kkw),
#' \code{\link{hskkw}} (Hessian for kkw),
#' \code{\link{dkkw}} (density for kkw),
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
#' true_params <- c(alpha = 2.0, beta = 3.0, delta = 1.5, lambda = 2.0)
#' data <- rkkw(n,
#'   alpha = true_params[1], beta = true_params[2],
#'   delta = true_params[3], lambda = true_params[4]
#' )
#'
#' # Evaluate gradient at true parameters
#' grad_true <- grkkw(par = true_params, data = data)
#' cat("Gradient at true parameters:\n")
#' print(grad_true)
#' cat("Norm:", sqrt(sum(grad_true^2)), "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 2.5, 1.0, 1.5),
#'   c(2.0, 3.0, 1.5, 2.0),
#'   c(2.5, 3.5, 2.0, 2.5)
#' )
#'
#' grad_norms <- apply(test_params, 1, function(p) {
#'   g <- grkkw(p, data)
#'   sqrt(sum(g^2))
#' })
#'
#' results <- data.frame(
#'   Alpha = test_params[, 1],
#'   Beta = test_params[, 2],
#'   Delta = test_params[, 3],
#'   Lambda = test_params[, 4],
#'   Grad_Norm = grad_norms
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Gradient in Optimization
#'
#' # Optimization with analytical gradient
#' fit_with_grad <- optim(
#'   par = c(1.5, 2.5, 1.0, 1.5),
#'   fn = llkkw,
#'   gr = grkkw,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE,
#'   control = list(trace = 0)
#' )
#'
#' # Optimization without gradient
#' fit_no_grad <- optim(
#'   par = c(1.5, 2.5, 1.0, 1.5),
#'   fn = llkkw,
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
#'   Delta = c(fit_with_grad$par[3], fit_no_grad$par[3]),
#'   Lambda = c(fit_with_grad$par[4], fit_no_grad$par[4]),
#'   NegLogLik = c(fit_with_grad$value, fit_no_grad$value),
#'   Iterations = c(fit_with_grad$counts[1], fit_no_grad$counts[1])
#' )
#' print(comparison, digits = 4, row.names = FALSE)
#'
#'
#' ## Example 3: Verifying Gradient at MLE
#'
#' mle <- fit_with_grad$par
#' names(mle) <- c("alpha", "beta", "delta", "lambda")
#'
#' # At MLE, gradient should be approximately zero
#' gradient_at_mle <- grkkw(par = mle, data = data)
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
#' grad_analytical <- grkkw(par = mle, data = data)
#' grad_numerical <- numerical_gradient(llkkw, mle, data)
#'
#' comparison_grad <- data.frame(
#'   Parameter = c("alpha", "beta", "delta", "lambda"),
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
#' theta0 <- c(1.8, 2.8, 1.3, 1.8)
#' score_theta0 <- -grkkw(par = theta0, data = data)
#'
#' # Fisher information at theta0
#' fisher_info <- hskkw(par = theta0, data = data)
#'
#' # Score test statistic
#' score_stat <- t(score_theta0) %*% solve(fisher_info) %*% score_theta0
#' p_value <- pchisq(score_stat, df = 4, lower.tail = FALSE)
#'
#' cat("\nScore Test:\n")
#' cat("H0: alpha=1.8, beta=2.8, delta=1.3, lambda=1.8\n")
#' cat("Test statistic:", score_stat, "\n")
#' cat("P-value:", format.pval(p_value, digits = 4), "\n")
#'
#'
#' ## Example 6: Confidence Ellipse with Gradient Information
#'
#' # For visualization, use first two parameters (alpha, beta)
#' # Observed information
#' obs_info <- hskkw(par = mle, data = data)
#' vcov_full <- solve(obs_info)
#' vcov_2d <- vcov_full[1:2, 1:2]
#'
#' # Create confidence ellipse
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' eig_decomp <- eigen(vcov_2d)
#' ellipse <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse[i, ] <- mle[1:2] + sqrt(chi2_val) *
#'     (eig_decomp$vectors %*% diag(sqrt(eig_decomp$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_2d <- sqrt(diag(vcov_2d))
#' ci_alpha <- mle[1] + c(-1, 1) * 1.96 * se_2d[1]
#' ci_beta <- mle[2] + c(-1, 1) * 1.96 * se_2d[2]
#'
#' # Plot
#' plot(ellipse[, 1], ellipse[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "95% Confidence Region (Alpha vs Beta)", las = 1
#' )
#'
#' # Add marginal CIs
#' abline(v = ci_alpha, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_beta, col = "#808080", lty = 3, lwd = 1.5)
#'
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#'
#' legend("topright",
#'   legend = c("MLE", "True", "95% CR", "Marginal 95% CI"),
#'   col = c("#8B0000", "#006400", "#2E4057", "#808080"),
#'   pch = c(19, 17, NA, NA),
#'   lty = c(NA, NA, 1, 3),
#'   lwd = c(NA, NA, 2, 1.5),
#'   bty = "n"
#' )
#' grid(col = "gray90")
#' }
#'
#' @export
grkkw <- function(par, data) {
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
  .Call("_gkwdist_grkkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 7. HESSIAN (hskkw)
# ----------------------------------------------------------------------------#

#' @title Hessian Matrix of the Negative Log-Likelihood for the kkw Distribution
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize hessian
#'
#' @description
#' Computes the analytic 4x4 Hessian matrix (matrix of second partial derivatives)
#' of the negative log-likelihood function for the Kumaraswamy-Kumaraswamy (kkw)
#' distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
#' (\eqn{\beta}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}).
#' This distribution is the special case of the Generalized Kumaraswamy (GKw)
#' distribution where \eqn{\gamma = 1}. The Hessian is useful for estimating
#' standard errors and in optimization algorithms.
#'
#' @param par A numeric vector of length 4 containing the distribution parameters
#'   in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
#'   \code{delta} (\eqn{\delta \ge 0}), \code{lambda} (\eqn{\lambda > 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a 4x4 numeric matrix representing the Hessian matrix of the
#'   negative log-likelihood function, \eqn{-\partial^2 \ell / (\partial \theta_i \partial \theta_j)},
#'   where \eqn{\theta = (\alpha, \beta, \delta, \lambda)}.
#'   Returns a 4x4 matrix populated with \code{NaN} if any parameter values are
#'   invalid according to their constraints, or if any value in \code{data} is
#'   not in the interval (0, 1).
#'
#' @details
#' This function calculates the analytic second partial derivatives of the
#' negative log-likelihood function based on the kkw log-likelihood
#' (\eqn{\gamma=1} case of GKw, see \code{\link{llkkw}}):
#' \deqn{
#' \ell(\theta | \mathbf{x}) = n[\ln(\delta+1) + \ln(\lambda) + \ln(\alpha) + \ln(\beta)]
#' + \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta-1)\ln(v_i) + (\lambda-1)\ln(w_i) + \delta\ln(z_i)]
#' }
#' where \eqn{\theta = (\alpha, \beta, \delta, \lambda)} and intermediate terms are:
#' \itemize{
#'   \item \eqn{v_i = 1 - x_i^{\alpha}}
#'   \item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
#'   \item \eqn{z_i = 1 - w_i^{\lambda} = 1 - [1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
#' }
#' The Hessian matrix returned contains the elements \eqn{- \frac{\partial^2 \ell(\theta | \mathbf{x})}{\partial \theta_i \partial \theta_j}}
#' for \eqn{\theta_i, \theta_j \in \{\alpha, \beta, \delta, \lambda\}}.
#'
#' Key properties of the returned matrix:
#' \itemize{
#'   \item Dimensions: 4x4.
#'   \item Symmetry: The matrix is symmetric.
#'   \item Ordering: Rows and columns correspond to the parameters in the order
#'     \eqn{\alpha, \beta, \delta, \lambda}.
#'   \item Content: Analytic second derivatives of the *negative* log-likelihood.
#' }
#' This corresponds to the relevant submatrix of the 5x5 GKw Hessian (\code{\link{hsgkw}})
#' evaluated at \eqn{\gamma=1}. The exact analytical formulas are implemented directly.
#'
#' @references
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, *46*(1-2), 79-88.
#'
#' @seealso
#' \code{\link{hsgkw}} (parent distribution Hessian),
#' \code{\link{llkkw}} (negative log-likelihood for kkw),
#' \code{\link{grkkw}} (gradient for kkw),
#' \code{\link{dkkw}} (density for kkw),
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
#' true_params <- c(alpha = 2.0, beta = 3.0, delta = 1.5, lambda = 2.0)
#' data <- rkkw(n,
#'   alpha = true_params[1], beta = true_params[2],
#'   delta = true_params[3], lambda = true_params[4]
#' )
#'
#' # Evaluate Hessian at true parameters
#' hess_true <- hskkw(par = true_params, data = data)
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
#'   par = c(1.5, 2.5, 1.0, 1.5),
#'   fn = llkkw,
#'   gr = grkkw,
#'   data = data,
#'   method = "BFGS",
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("alpha", "beta", "delta", "lambda")
#'
#' # Hessian at MLE
#' hessian_at_mle <- hskkw(par = mle, data = data)
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
#' names(se) <- c("alpha", "beta", "delta", "lambda")
#'
#' # Correlation matrix
#' corr_matrix <- cov2cor(vcov_matrix)
#' cat("\nCorrelation Matrix:\n")
#' print(corr_matrix, digits = 4)
#'
#' # Confidence intervals
#' z_crit <- qnorm(0.975)
#' results <- data.frame(
#'   Parameter = c("alpha", "beta", "delta", "lambda"),
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
#'   c(1.5, 2.5, 1.0, 1.5),
#'   c(2.0, 3.0, 1.5, 2.0),
#'   mle,
#'   c(2.5, 3.5, 2.0, 2.5)
#' )
#'
#' hess_properties <- data.frame(
#'   Alpha = numeric(),
#'   Beta = numeric(),
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
#'   H <- hskkw(par = test_params[i, ], data = data)
#'   eigs <- eigen(H, only.values = TRUE)$values
#'
#'   hess_properties <- rbind(hess_properties, data.frame(
#'     Alpha = test_params[i, 1],
#'     Beta = test_params[i, 2],
#'     Delta = test_params[i, 3],
#'     Lambda = test_params[i, 4],
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
#' ## Example 5: Curvature Visualization (Alpha vs Beta)
#'
#' # Create grid around MLE
#' alpha_grid <- seq(mle[1] - 1, mle[1] + 1, length.out = round(n / 4))
#' beta_grid <- seq(mle[2] - 1, mle[2] + 1, length.out = round(n / 4))
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
#'     H <- hskkw(c(alpha_grid[i], beta_grid[j], mle[3], mle[4]), data)
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
#' ## Example 6: Confidence Ellipse (Alpha vs Beta)
#'
#' # Extract 2x2 submatrix for alpha and beta
#' vcov_2d <- vcov_matrix[1:2, 1:2]
#'
#' # Create confidence ellipse
#' theta <- seq(0, 2 * pi, length.out = round(n / 2))
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' eig_decomp <- eigen(vcov_2d)
#' ellipse <- matrix(NA, nrow = round(n / 2), ncol = 2)
#' for (i in 1:round(n / 2)) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse[i, ] <- mle[1:2] + sqrt(chi2_val) *
#'     (eig_decomp$vectors %*% diag(sqrt(eig_decomp$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_2d <- sqrt(diag(vcov_2d))
#' ci_alpha <- mle[1] + c(-1, 1) * 1.96 * se_2d[1]
#' ci_beta <- mle[2] + c(-1, 1) * 1.96 * se_2d[2]
#'
#' # Plot
#' plot(ellipse[, 1], ellipse[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(alpha), ylab = expression(beta),
#'   main = "95% Confidence Ellipse (Alpha vs Beta)", las = 1
#' )
#'
#' # Add marginal CIs
#' abline(v = ci_alpha, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_beta, col = "#808080", lty = 3, lwd = 1.5)
#'
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#'
#' legend("topright",
#'   legend = c("MLE", "True", "95% CR", "Marginal 95% CI"),
#'   col = c("#8B0000", "#006400", "#2E4057", "#808080"),
#'   pch = c(19, 17, NA, NA),
#'   lty = c(NA, NA, 1, 3),
#'   lwd = c(NA, NA, 2, 1.5),
#'   bty = "n"
#' )
#' grid(col = "gray90")
#'
#'
#' ## Example 7: Confidence Ellipse (Delta vs Lambda)
#'
#' # Extract 2x2 submatrix for delta and lambda
#' vcov_2d_dl <- vcov_matrix[3:4, 3:4]
#'
#' # Create confidence ellipse
#' eig_decomp_dl <- eigen(vcov_2d_dl)
#' ellipse_dl <- matrix(NA, nrow = round(n / 2), ncol = 2)
#' for (i in 1:round(n / 2)) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse_dl[i, ] <- mle[3:4] + sqrt(chi2_val) *
#'     (eig_decomp_dl$vectors %*% diag(sqrt(eig_decomp_dl$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_2d_dl <- sqrt(diag(vcov_2d_dl))
#' ci_delta <- mle[3] + c(-1, 1) * 1.96 * se_2d_dl[1]
#' ci_lambda <- mle[4] + c(-1, 1) * 1.96 * se_2d_dl[2]
#'
#' # Plot
#' plot(ellipse_dl[, 1], ellipse_dl[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(delta), ylab = expression(lambda),
#'   main = "95% Confidence Ellipse (Delta vs Lambda)", las = 1
#' )
#'
#' # Add marginal CIs
#' abline(v = ci_delta, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_lambda, col = "#808080", lty = 3, lwd = 1.5)
#'
#' points(mle[3], mle[4], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[3], true_params[4], pch = 17, col = "#006400", cex = 1.5)
#'
#' legend("topright",
#'   legend = c("MLE", "True", "95% CR", "Marginal 95% CI"),
#'   col = c("#8B0000", "#006400", "#2E4057", "#808080"),
#'   pch = c(19, 17, NA, NA),
#'   lty = c(NA, NA, 1, 3),
#'   lwd = c(NA, NA, 2, 1.5),
#'   bty = "n"
#' )
#' grid(col = "gray90")
#' }
#' @export
hskkw <- function(par, data) {
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
  .Call("_gkwdist_hskkw",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}
