# ============================================================================#
# BETA DISTRIBUTION (gamma, delta+1 Parameterization)
# ============================================================================#
#
# Wrapper functions for the Beta distribution using the parameterization
# common in the Generalized Kumaraswamy family.
#
# IMPORTANT PARAMETERIZATION NOTE:
# This implementation uses gamma (γ) and delta (δ) parameters, where:
#   shape1 = gamma
#   shape2 = delta + 1
#
# This differs from stats::dbeta(x, shape1, shape2) where shape2 is used directly.
#
# Relationship to stats::Beta:
#   dbeta_(x, gamma, delta) ≡ dbeta(x, shape1 = gamma, shape2 = delta + 1)
#
# C++ implementations are in src/beta.cpp
#
# Functions:
#   - dbeta_: Probability density function (PDF)
#   - pbeta_: Cumulative distribution function (CDF)
#   - qbeta_: Quantile function (inverse CDF)
#   - rbeta_: Random number generation
#   - llbeta: Negative log-likelihood
#   - grbeta: Gradient of negative log-likelihood
#   - hsbeta: Hessian of negative log-likelihood
# ============================================================================#


# ----------------------------------------------------------------------------#
# 1. DENSITY FUNCTION (dbeta_)
# ----------------------------------------------------------------------------#

#' @title Density of the Beta Distribution (gamma, delta+1 Parameterization)
#' @author Lopes, J. E.
#' @keywords distribution density beta
#'
#' @description
#' Computes the probability density function (PDF) for the standard Beta
#' distribution, using a parameterization common in generalized distribution
#' families. The distribution is parameterized by \code{gamma} (\eqn{\gamma}) and
#' \code{delta} (\eqn{\delta}), corresponding to the standard Beta distribution
#' with shape parameters \code{shape1 = gamma} and \code{shape2 = delta + 1}.
#' The distribution is defined on the interval (0, 1).
#'
#' @param x Vector of quantiles (values between 0 and 1).
#' @param gamma First shape parameter (\code{shape1}), \eqn{\gamma > 0}. Can be a
#'   scalar or a vector. Default: 1.0.
#' @param delta Second shape parameter is \code{delta + 1} (\code{shape2}), requires
#'   \eqn{\delta \ge 0} so that \code{shape2 >= 1}. Can be a scalar or a vector.
#'   Default: 0.0 (leading to \code{shape2 = 1}).
#' @param log Logical; if \code{TRUE}, the logarithm of the density is
#'   returned (\eqn{\log(f(x))}). Default: \code{FALSE}.
#'
#' @return A vector of density values (\eqn{f(x)}) or log-density values
#'   (\eqn{\log(f(x))}). The length of the result is determined by the recycling
#'   rule applied to the arguments (\code{x}, \code{gamma}, \code{delta}).
#'   Returns \code{0} (or \code{-Inf} if \code{log = TRUE}) for \code{x}
#'   outside the interval (0, 1), or \code{NaN} if parameters are invalid
#'   (e.g., \code{gamma <= 0}, \code{delta < 0}).
#'
#' @details
#' The probability density function (PDF) calculated by this function corresponds
#' to a standard Beta distribution \eqn{Beta(\gamma, \delta+1)}:
#' \deqn{
#' f(x; \gamma, \delta) = \frac{x^{\gamma-1} (1-x)^{(\delta+1)-1}}{B(\gamma, \delta+1)} = \frac{x^{\gamma-1} (1-x)^{\delta}}{B(\gamma, \delta+1)}
#' }
#' for \eqn{0 < x < 1}, where \eqn{B(a,b)} is the Beta function
#' (\code{\link[base]{beta}}).
#'
#' This specific parameterization arises as a special case of the five-parameter
#' Generalized Kumaraswamy (GKw) distribution (\code{\link{dgkw}}) obtained
#' by setting the parameters \eqn{\alpha = 1}, \eqn{\beta = 1}, and \eqn{\lambda = 1}.
#' It is therefore equivalent to the McDonald (Mc)/Beta Power distribution
#' (\code{\link{dmc}}) with \eqn{\lambda = 1}.
#'
#' Note the difference in the second parameter compared to \code{\link[stats]{dbeta}},
#' where \code{dbeta(x, shape1, shape2)} uses \code{shape2} directly. Here,
#' \code{shape1 = gamma} and \code{shape2 = delta + 1}.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous Univariate
#' Distributions, Volume 2* (2nd ed.). Wiley.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' @seealso
#' \code{\link[stats]{dbeta}} (standard R implementation),
#' \code{\link{dgkw}} (parent distribution density),
#' \code{\link{dmc}} (McDonald/Beta Power density),
#' \code{pbeta_}, \code{qbeta_}, \code{rbeta_} (other functions for this parameterization, if they exist).
#'
#' @examples
#' \donttest{
#' # Example values
#' x_vals <- c(0.2, 0.5, 0.8)
#' gamma_par <- 2.0 # Corresponds to shape1
#' delta_par <- 3.0 # Corresponds to shape2 - 1
#' shape1 <- gamma_par
#' shape2 <- delta_par + 1
#'
#' # Calculate density using dbeta_
#' densities <- dbeta_(x_vals, gamma_par, delta_par)
#' print(densities)
#'
#' # Compare with stats::dbeta
#' densities_stats <- stats::dbeta(x_vals, shape1 = shape1, shape2 = shape2)
#' print(paste("Max difference vs stats::dbeta:", max(abs(densities - densities_stats))))
#'
#' # Compare with dgkw setting alpha=1, beta=1, lambda=1
#' densities_gkw <- dgkw(x_vals,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print(paste("Max difference vs dgkw:", max(abs(densities - densities_gkw))))
#'
#' # Compare with dmc setting lambda=1
#' densities_mc <- dmc(x_vals, gamma = gamma_par, delta = delta_par, lambda = 1.0)
#' print(paste("Max difference vs dmc:", max(abs(densities - densities_mc))))
#'
#' # Calculate log-density
#' log_densities <- dbeta_(x_vals, gamma_par, delta_par, log = TRUE)
#' print(log_densities)
#' print(stats::dbeta(x_vals, shape1 = shape1, shape2 = shape2, log = TRUE))
#'
#' # Plot the density
#' curve_x <- seq(0.001, 0.999, length.out = 200)
#' curve_y <- dbeta_(curve_x, gamma = 2, delta = 3) # Beta(2, 4)
#' plot(curve_x, curve_y,
#'   type = "l", main = "Beta(2, 4) Density via dbeta_",
#'   xlab = "x", ylab = "f(x)", col = "blue"
#' )
#' curve(stats::dbeta(x, 2, 4), add = TRUE, col = "red", lty = 2)
#' legend("topright",
#'   legend = c("dbeta_(gamma=2, delta=3)", "stats::dbeta(shape1=2, shape2=4)"),
#'   col = c("blue", "red"), lty = c(1, 2), bty = "n"
#' )
#' }
#'
#' @export
dbeta_ <- function(x, gamma = 1, delta = 0, log = FALSE) {
  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }
  if (!is.numeric(gamma) || any(gamma <= 0, na.rm = TRUE)) {
    stop("'gamma' must be positive (gamma > 0)")
  }
  if (!is.numeric(delta) || any(delta < 0, na.rm = TRUE)) {
    stop("'delta' must be non-negative (delta >= 0)")
  }
  if (!is.logical(log) || length(log) != 1) {
    stop("'log' must be a single logical value")
  }

  # Informative warning about parameterization
  if (getOption("gkwdist.warn.beta.param", default = FALSE)) {
    message(
      "Note: dbeta_(x, gamma, delta) uses shape1=gamma, shape2=delta+1\n",
      "Equivalent to: dbeta(x, shape1=", gamma[1], ", shape2=", delta[1] + 1, ")\n",
      "Set options(gkwdist.warn.beta.param=FALSE) to suppress this message."
    )
  }

  # Call C++ implementation
  .Call("_gkwdist_dbeta_",
    as.numeric(x),
    as.numeric(gamma),
    as.numeric(delta),
    as.logical(log),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 2. DISTRIBUTION FUNCTION (pbeta_)
# ----------------------------------------------------------------------------#

#' @title CDF of the Beta Distribution (gamma, delta+1 Parameterization)
#' @author Lopes, J. E.
#' @keywords distribution cumulative beta
#'
#' @description
#' Computes the cumulative distribution function (CDF), \eqn{F(q) = P(X \le q)},
#' for the standard Beta distribution, using a parameterization common in
#' generalized distribution families. The distribution is parameterized by
#' \code{gamma} (\eqn{\gamma}) and \code{delta} (\eqn{\delta}), corresponding to
#' the standard Beta distribution with shape parameters \code{shape1 = gamma}
#' and \code{shape2 = delta + 1}.
#'
#' @param q Vector of quantiles (values generally between 0 and 1).
#' @param gamma First shape parameter (\code{shape1}), \eqn{\gamma > 0}. Can be a
#'   scalar or a vector. Default: 1.0.
#' @param delta Second shape parameter is \code{delta + 1} (\code{shape2}), requires
#'   \eqn{\delta \ge 0} so that \code{shape2 >= 1}. Can be a scalar or a vector.
#'   Default: 0.0 (leading to \code{shape2 = 1}).
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}, otherwise, \eqn{P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of probabilities, \eqn{F(q)}, or their logarithms/complements
#'   depending on \code{lower.tail} and \code{log.p}. The length of the result
#'   is determined by the recycling rule applied to the arguments (\code{q},
#'   \code{gamma}, \code{delta}). Returns \code{0} (or \code{-Inf} if
#'   \code{log.p = TRUE}) for \code{q <= 0} and \code{1} (or \code{0} if
#'   \code{log.p = TRUE}) for \code{q >= 1}. Returns \code{NaN} for invalid
#'   parameters.
#'
#' @details
#' This function computes the CDF of a Beta distribution with parameters
#' \code{shape1 = gamma} and \code{shape2 = delta + 1}. It is equivalent to
#' calling \code{stats::pbeta(q, shape1 = gamma, shape2 = delta + 1,
#' lower.tail = lower.tail, log.p = log.p)}.
#'
#' This distribution arises as a special case of the five-parameter
#' Generalized Kumaraswamy (GKw) distribution (\code{\link{pgkw}}) obtained
#' by setting \eqn{\alpha = 1}, \eqn{\beta = 1}, and \eqn{\lambda = 1}.
#' It is therefore also equivalent to the McDonald (Mc)/Beta Power distribution
#' (\code{\link{pmc}}) with \eqn{\lambda = 1}.
#'
#' The function likely calls R's underlying \code{pbeta} function but ensures
#' consistent parameter recycling and handling within the C++ environment,
#' matching the style of other functions in the related families.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous Univariate
#' Distributions, Volume 2* (2nd ed.). Wiley.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' @seealso
#' \code{\link[stats]{pbeta}} (standard R implementation),
#' \code{\link{pgkw}} (parent distribution CDF),
#' \code{\link{pmc}} (McDonald/Beta Power CDF),
#' \code{dbeta_}, \code{qbeta_}, \code{rbeta_} (other functions for this parameterization, if they exist).
#'
#' @examples
#' \donttest{
#' # Example values
#' q_vals <- c(0.2, 0.5, 0.8)
#' gamma_par <- 2.0 # Corresponds to shape1
#' delta_par <- 3.0 # Corresponds to shape2 - 1
#' shape1 <- gamma_par
#' shape2 <- delta_par + 1
#'
#' # Calculate CDF using pbeta_
#' probs <- pbeta_(q_vals, gamma_par, delta_par)
#' print(probs)
#'
#' # Compare with stats::pbeta
#' probs_stats <- stats::pbeta(q_vals, shape1 = shape1, shape2 = shape2)
#' print(paste("Max difference vs stats::pbeta:", max(abs(probs - probs_stats))))
#'
#' # Compare with pgkw setting alpha=1, beta=1, lambda=1
#' probs_gkw <- pgkw(q_vals,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print(paste("Max difference vs pgkw:", max(abs(probs - probs_gkw))))
#'
#' # Compare with pmc setting lambda=1
#' probs_mc <- pmc(q_vals, gamma = gamma_par, delta = delta_par, lambda = 1.0)
#' print(paste("Max difference vs pmc:", max(abs(probs - probs_mc))))
#'
#' # Calculate upper tail P(X > q)
#' probs_upper <- pbeta_(q_vals, gamma_par, delta_par, lower.tail = FALSE)
#' print(probs_upper)
#' print(stats::pbeta(q_vals, shape1, shape2, lower.tail = FALSE))
#'
#' # Calculate log CDF
#' log.probs <- pbeta_(q_vals, gamma_par, delta_par, log.p = TRUE)
#' print(log.probs)
#' print(stats::pbeta(q_vals, shape1, shape2, log.p = TRUE))
#'
#' # Plot the CDF
#' curve_q <- seq(0.001, 0.999, length.out = 200)
#' curve_p <- pbeta_(curve_q, gamma = 2, delta = 3) # Beta(2, 4)
#' plot(curve_q, curve_p,
#'   type = "l", main = "Beta(2, 4) CDF via pbeta_",
#'   xlab = "q", ylab = "F(q)", col = "blue"
#' )
#' curve(stats::pbeta(x, 2, 4), add = TRUE, col = "red", lty = 2)
#' legend("bottomright",
#'   legend = c("pbeta_(gamma=2, delta=3)", "stats::pbeta(shape1=2, shape2=4)"),
#'   col = c("blue", "red"), lty = c(1, 2), bty = "n"
#' )
#' }
#'
#' @export
pbeta_ <- function(q, gamma = 1, delta = 0, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(q)) {
    stop("'q' must be numeric")
  }
  if (!is.numeric(gamma) || any(gamma <= 0, na.rm = TRUE)) {
    stop("'gamma' must be positive (gamma > 0)")
  }
  if (!is.numeric(delta) || any(delta < 0, na.rm = TRUE)) {
    stop("'delta' must be non-negative (delta >= 0)")
  }
  if (!is.logical(lower.tail) || length(lower.tail) != 1) {
    stop("'lower.tail' must be a single logical value")
  }
  if (!is.logical(log.p) || length(log.p) != 1) {
    stop("'log.p' must be a single logical value")
  }

  # Call C++ implementation
  .Call("_gkwdist_pbeta_",
    as.numeric(q),
    as.numeric(gamma),
    as.numeric(delta),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 3. QUANTILE FUNCTION (qbeta_)
# ----------------------------------------------------------------------------#

#' @title Quantile Function of the Beta Distribution (gamma, delta+1 Parameterization)
#' @author Lopes, J. E.
#' @keywords distribution quantile beta
#'
#' @description
#' Computes the quantile function (inverse CDF) for the standard Beta
#' distribution, using a parameterization common in generalized distribution
#' families. It finds the value \code{q} such that \eqn{P(X \le q) = p}. The
#' distribution is parameterized by \code{gamma} (\eqn{\gamma}) and \code{delta}
#' (\eqn{\delta}), corresponding to the standard Beta distribution with shape
#' parameters \code{shape1 = gamma} and \code{shape2 = delta + 1}.
#'
#' @param p Vector of probabilities (values between 0 and 1).
#' @param gamma First shape parameter (\code{shape1}), \eqn{\gamma > 0}. Can be a
#'   scalar or a vector. Default: 1.0.
#' @param delta Second shape parameter is \code{delta + 1} (\code{shape2}), requires
#'   \eqn{\delta \ge 0} so that \code{shape2 >= 1}. Can be a scalar or a vector.
#'   Default: 0.0 (leading to \code{shape2 = 1}).
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{p = P(X \le q)},
#'   otherwise, probabilities are \eqn{p = P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as
#'   \eqn{\log(p)}. Default: \code{FALSE}.
#'
#' @return A vector of quantiles corresponding to the given probabilities \code{p}.
#'   The length of the result is determined by the recycling rule applied to
#'   the arguments (\code{p}, \code{gamma}, \code{delta}).
#'   Returns:
#'   \itemize{
#'     \item \code{0} for \code{p = 0} (or \code{p = -Inf} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{1} for \code{p = 1} (or \code{p = 0} if \code{log.p = TRUE},
#'           when \code{lower.tail = TRUE}).
#'     \item \code{NaN} for \code{p < 0} or \code{p > 1} (or corresponding log scale).
#'     \item \code{NaN} for invalid parameters (e.g., \code{gamma <= 0},
#'           \code{delta < 0}).
#'   }
#'   Boundary return values are adjusted accordingly for \code{lower.tail = FALSE}.
#'
#' @details
#' This function computes the quantiles of a Beta distribution with parameters
#' \code{shape1 = gamma} and \code{shape2 = delta + 1}. It is equivalent to
#' calling \code{stats::qbeta(p, shape1 = gamma, shape2 = delta + 1,
#' lower.tail = lower.tail, log.p = log.p)}.
#'
#' This distribution arises as a special case of the five-parameter
#' Generalized Kumaraswamy (GKw) distribution (\code{\link{qgkw}}) obtained
#' by setting \eqn{\alpha = 1}, \eqn{\beta = 1}, and \eqn{\lambda = 1}.
#' It is therefore also equivalent to the McDonald (Mc)/Beta Power distribution
#' (\code{\link{qmc}}) with \eqn{\lambda = 1}.
#'
#' The function likely calls R's underlying \code{qbeta} function but ensures
#' consistent parameter recycling and handling within the C++ environment,
#' matching the style of other functions in the related families. Boundary
#' conditions (p=0, p=1) are handled explicitly.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous Univariate
#' Distributions, Volume 2* (2nd ed.). Wiley.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' @seealso
#' \code{\link[stats]{qbeta}} (standard R implementation),
#' \code{\link{qgkw}} (parent distribution quantile function),
#' \code{\link{qmc}} (McDonald/Beta Power quantile function),
#' \code{dbeta_}, \code{pbeta_}, \code{rbeta_} (other functions for this parameterization, if they exist).
#'
#' @examples
#' \donttest{
#' # Example values
#' p_vals <- c(0.1, 0.5, 0.9)
#' gamma_par <- 2.0 # Corresponds to shape1
#' delta_par <- 3.0 # Corresponds to shape2 - 1
#' shape1 <- gamma_par
#' shape2 <- delta_par + 1
#'
#' # Calculate quantiles using qbeta_
#' quantiles <- qbeta_(p_vals, gamma_par, delta_par)
#' print(quantiles)
#'
#' # Compare with stats::qbeta
#' quantiles_stats <- stats::qbeta(p_vals, shape1 = shape1, shape2 = shape2)
#' print(paste("Max difference vs stats::qbeta:", max(abs(quantiles - quantiles_stats))))
#'
#' # Compare with qgkw setting alpha=1, beta=1, lambda=1
#' quantiles_gkw <- qgkw(p_vals,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print(paste("Max difference vs qgkw:", max(abs(quantiles - quantiles_gkw))))
#'
#' # Compare with qmc setting lambda=1
#' quantiles_mc <- qmc(p_vals, gamma = gamma_par, delta = delta_par, lambda = 1.0)
#' print(paste("Max difference vs qmc:", max(abs(quantiles - quantiles_mc))))
#'
#' # Calculate quantiles for upper tail
#' quantiles_upper <- qbeta_(p_vals, gamma_par, delta_par, lower.tail = FALSE)
#' print(quantiles_upper)
#' print(stats::qbeta(p_vals, shape1, shape2, lower.tail = FALSE))
#'
#' # Calculate quantiles from log probabilities
#' log.p_vals <- log(p_vals)
#' quantiles_logp <- qbeta_(log.p_vals, gamma_par, delta_par, log.p = TRUE)
#' print(quantiles_logp)
#' print(stats::qbeta(log.p_vals, shape1, shape2, log.p = TRUE))
#'
#' # Verify inverse relationship with pbeta_
#' p_check <- 0.75
#' q_calc <- qbeta_(p_check, gamma_par, delta_par)
#' p_recalc <- pbeta_(q_calc, gamma_par, delta_par)
#' print(paste("Original p:", p_check, " Recalculated p:", p_recalc))
#' # abs(p_check - p_recalc) < 1e-9 # Should be TRUE
#'
#' # Boundary conditions
#' print(qbeta_(c(0, 1), gamma_par, delta_par)) # Should be 0, 1
#' print(qbeta_(c(-Inf, 0), gamma_par, delta_par, log.p = TRUE)) # Should be 0, 1
#' }
#'
#' @export
qbeta_ <- function(p, gamma = 1, delta = 0, lower.tail = TRUE, log.p = FALSE) {
  # Input validation
  if (!is.numeric(p)) {
    stop("'p' must be numeric")
  }
  if (!is.numeric(gamma) || any(gamma <= 0, na.rm = TRUE)) {
    stop("'gamma' must be positive (gamma > 0)")
  }
  if (!is.numeric(delta) || any(delta < 0, na.rm = TRUE)) {
    stop("'delta' must be non-negative (delta >= 0)")
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
  .Call("_gkwdist_qbeta_",
    as.numeric(p),
    as.numeric(gamma),
    as.numeric(delta),
    as.logical(lower.tail),
    as.logical(log.p),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 4. RANDOM GENERATION (rbeta_)
# ----------------------------------------------------------------------------#

#' @title Random Generation for the Beta Distribution (gamma, delta+1 Parameterization)
#' @author Lopes, J. E.
#' @keywords distribution random beta
#'
#' @description
#' Generates random deviates from the standard Beta distribution, using a
#' parameterization common in generalized distribution families. The distribution
#' is parameterized by \code{gamma} (\eqn{\gamma}) and \code{delta} (\eqn{\delta}),
#' corresponding to the standard Beta distribution with shape parameters
#' \code{shape1 = gamma} and \code{shape2 = delta + 1}. This is a special case
#' of the Generalized Kumaraswamy (GKw) distribution where \eqn{\alpha = 1},
#' \eqn{\beta = 1}, and \eqn{\lambda = 1}.
#'
#' @param n Number of observations. If \code{length(n) > 1}, the length is
#'   taken to be the number required. Must be a non-negative integer.
#' @param gamma First shape parameter (\code{shape1}), \eqn{\gamma > 0}. Can be a
#'   scalar or a vector. Default: 1.0.
#' @param delta Second shape parameter is \code{delta + 1} (\code{shape2}), requires
#'   \eqn{\delta \ge 0} so that \code{shape2 >= 1}. Can be a scalar or a vector.
#'   Default: 0.0 (leading to \code{shape2 = 1}, i.e., Uniform).
#'
#' @return A numeric vector of length \code{n} containing random deviates from the
#'   Beta(\eqn{\gamma, \delta+1}) distribution, with values in (0, 1). The length
#'   of the result is determined by \code{n} and the recycling rule applied to
#'   the parameters (\code{gamma}, \code{delta}). Returns \code{NaN} if parameters
#'   are invalid (e.g., \code{gamma <= 0}, \code{delta < 0}).
#'
#' @details
#' This function generates samples from a Beta distribution with parameters
#' \code{shape1 = gamma} and \code{shape2 = delta + 1}. It is equivalent to
#' calling \code{stats::rbeta(n, shape1 = gamma, shape2 = delta + 1)}.
#'
#' This distribution arises as a special case of the five-parameter
#' Generalized Kumaraswamy (GKw) distribution (\code{\link{rgkw}}) obtained
#' by setting \eqn{\alpha = 1}, \eqn{\beta = 1}, and \eqn{\lambda = 1}.
#' It is therefore also equivalent to the McDonald (Mc)/Beta Power distribution
#' (\code{\link{rmc}}) with \eqn{\lambda = 1}.
#'
#' The function likely calls R's underlying \code{rbeta} function but ensures
#' consistent parameter recycling and handling within the C++ environment,
#' matching the style of other functions in the related families.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous Univariate
#' Distributions, Volume 2* (2nd ed.). Wiley.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer-Verlag.
#'
#' @seealso
#' \code{\link[stats]{rbeta}} (standard R implementation),
#' \code{\link{rgkw}} (parent distribution random generation),
#' \code{\link{rmc}} (McDonald/Beta Power random generation),
#' \code{dbeta_}, \code{pbeta_}, \code{qbeta_} (other functions for this parameterization, if they exist).
#'
#' @examples
#' \donttest{
#' set.seed(2030) # for reproducibility
#'
#' # Generate 1000 samples using rbeta_
#' gamma_par <- 2.0 # Corresponds to shape1
#' delta_par <- 3.0 # Corresponds to shape2 - 1
#' shape1 <- gamma_par
#' shape2 <- delta_par + 1
#'
#' x_sample <- rbeta_(1000, gamma = gamma_par, delta = delta_par)
#' summary(x_sample)
#'
#' # Compare with stats::rbeta
#' x_sample_stats <- stats::rbeta(1000, shape1 = shape1, shape2 = shape2)
#' # Visually compare histograms or QQ-plots
#' hist(x_sample, main = "rbeta_ Sample", freq = FALSE, breaks = 30)
#' curve(dbeta_(x, gamma_par, delta_par), add = TRUE, col = "red", lwd = 2)
#' hist(x_sample_stats, main = "stats::rbeta Sample", freq = FALSE, breaks = 30)
#' curve(stats::dbeta(x, shape1, shape2), add = TRUE, col = "blue", lwd = 2)
#' # Compare summary stats (should be similar due to randomness)
#' print(summary(x_sample))
#' print(summary(x_sample_stats))
#'
#' # Compare summary stats with rgkw(alpha=1, beta=1, lambda=1)
#' x_sample_gkw <- rgkw(1000,
#'   alpha = 1.0, beta = 1.0, gamma = gamma_par,
#'   delta = delta_par, lambda = 1.0
#' )
#' print("Summary stats for rgkw(a=1,b=1,l=1) sample:")
#' print(summary(x_sample_gkw))
#'
#' # Compare summary stats with rmc(lambda=1)
#' x_sample_mc <- rmc(1000, gamma = gamma_par, delta = delta_par, lambda = 1.0)
#' print("Summary stats for rmc(l=1) sample:")
#' print(summary(x_sample_mc))
#' }
#'
#' @export
rbeta_ <- function(n, gamma = 1, delta = 0) {
  # Input validation
  if (length(n) > 1) n <- length(n)
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("'n' must be a positive integer")
  }
  n <- as.integer(n)

  if (!is.numeric(gamma) || any(gamma <= 0, na.rm = TRUE)) {
    stop("'gamma' must be positive (gamma > 0)")
  }
  if (!is.numeric(delta) || any(delta < 0, na.rm = TRUE)) {
    stop("'delta' must be non-negative (delta >= 0)")
  }

  # Call C++ implementation
  .Call("_gkwdist_rbeta_",
    as.integer(n),
    as.numeric(gamma),
    as.numeric(delta),
    PACKAGE = "gkwdist"
  )
}


# ============================================================================#
# MAXIMUM LIKELIHOOD ESTIMATION FUNCTIONS
# ============================================================================#

# ----------------------------------------------------------------------------#
# 5. NEGATIVE LOG-LIKELIHOOD (llbeta)
# ----------------------------------------------------------------------------#

#' @title Negative Log-Likelihood for the Beta Distribution (gamma, delta+1 Parameterization)
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize beta
#'
#' @description
#' Computes the negative log-likelihood function for the standard Beta
#' distribution, using a parameterization common in generalized distribution
#' families. The distribution is parameterized by \code{gamma} (\eqn{\gamma}) and
#' \code{delta} (\eqn{\delta}), corresponding to the standard Beta distribution
#' with shape parameters \code{shape1 = gamma} and \code{shape2 = delta + 1}.
#' This function is suitable for maximum likelihood estimation.
#'
#' @param par A numeric vector of length 2 containing the distribution parameters
#'   in the order: \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a single \code{double} value representing the negative
#'   log-likelihood (\eqn{-\ell(\theta|\mathbf{x})}). Returns \code{Inf}
#'   if any parameter values in \code{par} are invalid according to their
#'   constraints, or if any value in \code{data} is not in the interval (0, 1).
#'
#' @details
#' This function calculates the negative log-likelihood for a Beta distribution
#' with parameters \code{shape1 = gamma} (\eqn{\gamma}) and \code{shape2 = delta + 1} (\eqn{\delta+1}).
#' The probability density function (PDF) is:
#' \deqn{
#' f(x | \gamma, \delta) = \frac{x^{\gamma-1} (1-x)^{\delta}}{B(\gamma, \delta+1)}
#' }
#' for \eqn{0 < x < 1}, where \eqn{B(a,b)} is the Beta function (\code{\link[base]{beta}}).
#' The log-likelihood function \eqn{\ell(\theta | \mathbf{x})} for a sample
#' \eqn{\mathbf{x} = (x_1, \dots, x_n)} is \eqn{\sum_{i=1}^n \ln f(x_i | \theta)}:
#' \deqn{
#' \ell(\theta | \mathbf{x}) = \sum_{i=1}^{n} [(\gamma-1)\ln(x_i) + \delta\ln(1-x_i)] - n \ln B(\gamma, \delta+1)
#' }
#' where \eqn{\theta = (\gamma, \delta)}.
#' This function computes and returns the *negative* log-likelihood, \eqn{-\ell(\theta|\mathbf{x})},
#' suitable for minimization using optimization routines like \code{\link[stats]{optim}}.
#' It is equivalent to the negative log-likelihood of the GKw distribution
#' (\code{\link{llgkw}}) evaluated at \eqn{\alpha=1, \beta=1, \lambda=1}, and also
#' to the negative log-likelihood of the McDonald distribution (\code{\link{llmc}})
#' evaluated at \eqn{\lambda=1}. The term \eqn{\ln B(\gamma, \delta+1)} is typically
#' computed using log-gamma functions (\code{\link[base]{lgamma}}) for numerical stability.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous Univariate
#' Distributions, Volume 2* (2nd ed.). Wiley.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#'
#' @seealso
#' \code{\link{llgkw}}, \code{\link{llmc}} (related negative log-likelihoods),
#' \code{dbeta_}, \code{pbeta_}, \code{qbeta_}, \code{rbeta_},
#' \code{grbeta} (gradient, if available),
#' \code{hsbeta} (Hessian, if available),
#' \code{\link[stats]{optim}}, \code{\link[base]{lbeta}}.
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Log-Likelihood Evaluation
#'
#' # Generate sample data
#' set.seed(123)
#' n <- 1000
#' true_params <- c(gamma = 2.0, delta = 3.0)
#' data <- rbeta_(n, gamma = true_params[1], delta = true_params[2])
#'
#' # Evaluate negative log-likelihood at true parameters
#' nll_true <- llbeta(par = true_params, data = data)
#' cat("Negative log-likelihood at true parameters:", nll_true, "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 2.5),
#'   c(2.0, 3.0),
#'   c(2.5, 3.5)
#' )
#'
#' nll_values <- apply(test_params, 1, function(p) llbeta(p, data))
#' results <- data.frame(
#'   Gamma = test_params[, 1],
#'   Delta = test_params[, 2],
#'   NegLogLik = nll_values
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Maximum Likelihood Estimation
#'
#' # Optimization using L-BFGS-B with bounds
#' fit <- optim(
#'   par = c(1.5, 2.5),
#'   fn = llbeta,
#'   gr = grbeta,
#'   data = data,
#'   method = "L-BFGS-B",
#'   lower = c(0.01, 0.01),
#'   upper = c(100, 100),
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("gamma", "delta")
#' se <- sqrt(diag(solve(fit$hessian)))
#'
#' results <- data.frame(
#'   Parameter = c("gamma", "delta"),
#'   True = true_params,
#'   MLE = mle,
#'   SE = se,
#'   CI_Lower = mle - 1.96 * se,
#'   CI_Upper = mle + 1.96 * se
#' )
#' print(results, digits = 4)
#'
#' cat(sprintf(
#'   "\nMLE corresponds approx to Beta(%.2f, %.2f)\n",
#'   mle[1], mle[2] + 1
#' ))
#' cat(
#'   "True corresponds to Beta(%.2f, %.2f)\n",
#'   true_params[1], true_params[2] + 1
#' )
#'
#' cat("\nNegative log-likelihood at MLE:", fit$value, "\n")
#' cat("AIC:", 2 * fit$value + 2 * length(mle), "\n")
#' cat("BIC:", 2 * fit$value + length(mle) * log(n), "\n")
#'
#'
#' ## Example 3: Comparing Optimization Methods
#'
#' methods <- c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG")
#' start_params <- c(1.5, 2.5)
#'
#' comparison <- data.frame(
#'   Method = character(),
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
#'       fn = llbeta,
#'       gr = grbeta,
#'       data = data,
#'       method = method
#'     )
#'   } else if (method == "L-BFGS-B") {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llbeta,
#'       gr = grbeta,
#'       data = data,
#'       method = method,
#'       lower = c(0.01, 0.01),
#'       upper = c(100, 100)
#'     )
#'   } else {
#'     fit_temp <- optim(
#'       par = start_params,
#'       fn = llbeta,
#'       data = data,
#'       method = method
#'     )
#'   }
#'
#'   comparison <- rbind(comparison, data.frame(
#'     Method = method,
#'     Gamma = fit_temp$par[1],
#'     Delta = fit_temp$par[2],
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
#' # Test H0: delta = 3 vs H1: delta free
#' loglik_full <- -fit$value
#'
#' restricted_ll <- function(params_restricted, data, delta_fixed) {
#'   llbeta(par = c(params_restricted[1], delta_fixed), data = data)
#' }
#'
#' fit_restricted <- optim(
#'   par = mle[1],
#'   fn = restricted_ll,
#'   data = data,
#'   delta_fixed = 3,
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
#'     par = mle[2],
#'     fn = function(p) llbeta(c(gamma_grid[i], p), data),
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
#'     par = mle[1],
#'     fn = function(p) llbeta(c(p, delta_grid[i]), data),
#'     method = "BFGS"
#'   )
#'   profile_ll_delta[i] <- -profile_fit$value
#' }
#'
#' # 95% confidence threshold
#' chi_crit <- qchisq(0.95, df = 1)
#' threshold <- max(profile_ll_gamma) - chi_crit / 2
#'
#' # Plot
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
#'
#' ## Example 6: 2D Log-Likelihood Surface (Gamma vs Delta)
#'
#' # Create 2D grid with wider range (±1.5)
#' gamma_2d <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = round(n / 25))
#' delta_2d <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = round(n / 25))
#' gamma_2d <- gamma_2d[gamma_2d > 0]
#' delta_2d <- delta_2d[delta_2d > 0]
#'
#' # Compute log-likelihood surface
#' ll_surface_gd <- matrix(NA, nrow = length(gamma_2d), ncol = length(delta_2d))
#'
#' for (i in seq_along(gamma_2d)) {
#'   for (j in seq_along(delta_2d)) {
#'     ll_surface_gd[i, j] <- -llbeta(c(gamma_2d[i], delta_2d[j]), data)
#'   }
#' }
#'
#' # Confidence region levels
#' max_ll_gd <- max(ll_surface_gd, na.rm = TRUE)
#' levels_90_gd <- max_ll_gd - qchisq(0.90, df = 2) / 2
#' levels_95_gd <- max_ll_gd - qchisq(0.95, df = 2) / 2
#' levels_99_gd <- max_ll_gd - qchisq(0.99, df = 2) / 2
#'
#' # Plot contour
#'
#' contour(gamma_2d, delta_2d, ll_surface_gd,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "2D Log-Likelihood: Gamma vs Delta",
#'   levels = seq(min(ll_surface_gd, na.rm = TRUE), max_ll_gd, length.out = 20),
#'   col = "#2E4057", las = 1, lwd = 1
#' )
#'
#' contour(gamma_2d, delta_2d, ll_surface_gd,
#'   levels = c(levels_90_gd, levels_95_gd, levels_99_gd),
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
#' }
#'
#' @export
llbeta <- function(par, data) {
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
  if (any(data <= 0 | data >= 1, na.rm = TRUE)) {
    warning("'data' contains values outside (0, 1)")
  }

  # Call C++ implementation
  .Call("_gkwdist_llbeta",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 6. GRADIENT (grbeta)
# ----------------------------------------------------------------------------#

#' @title Gradient of the Negative Log-Likelihood for the Beta Distribution (gamma, delta+1 Parameterization)
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize gradient beta
#'
#' @description
#' Computes the gradient vector (vector of first partial derivatives) of the
#' negative log-likelihood function for the standard Beta distribution, using
#' a parameterization common in generalized distribution families. The
#' distribution is parameterized by \code{gamma} (\eqn{\gamma}) and \code{delta}
#' (\eqn{\delta}), corresponding to the standard Beta distribution with shape
#' parameters \code{shape1 = gamma} and \code{shape2 = delta + 1}.
#' The gradient is useful for optimization algorithms.
#'
#' @param par A numeric vector of length 2 containing the distribution parameters
#'   in the order: \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a numeric vector of length 2 containing the partial derivatives
#'   of the negative log-likelihood function \eqn{-\ell(\theta | \mathbf{x})} with
#'   respect to each parameter: \eqn{(-\partial \ell/\partial \gamma, -\partial \ell/\partial \delta)}.
#'   Returns a vector of \code{NaN} if any parameter values are invalid according
#'   to their constraints, or if any value in \code{data} is not in the
#'   interval (0, 1).
#'
#' @details
#' This function calculates the gradient of the negative log-likelihood for a
#' Beta distribution with parameters \code{shape1 = gamma} (\eqn{\gamma}) and
#' \code{shape2 = delta + 1} (\eqn{\delta+1}). The components of the gradient
#' vector (\eqn{-\nabla \ell(\theta | \mathbf{x})}) are:
#'
#' \deqn{
#' -\frac{\partial \ell}{\partial \gamma} = n[\psi(\gamma) - \psi(\gamma+\delta+1)] -
#' \sum_{i=1}^{n}\ln(x_i)
#' }
#' \deqn{
#' -\frac{\partial \ell}{\partial \delta} = n[\psi(\delta+1) - \psi(\gamma+\delta+1)] -
#' \sum_{i=1}^{n}\ln(1-x_i)
#' }
#'
#' where \eqn{\psi(\cdot)} is the digamma function (\code{\link[base]{digamma}}).
#' These formulas represent the derivatives of \eqn{-\ell(\theta)}, consistent with
#' minimizing the negative log-likelihood. They correspond to the relevant components
#' of the general GKw gradient (\code{\link{grgkw}}) evaluated at \eqn{\alpha=1, \beta=1, \lambda=1}.
#' Note the parameterization: the standard Beta shape parameters are \eqn{\gamma} and \eqn{\delta+1}.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous Univariate
#' Distributions, Volume 2* (2nd ed.). Wiley.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' (Note: Specific gradient formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{grgkw}}, \code{\link{grmc}} (related gradients),
#' \code{\link{llbeta}} (negative log-likelihood function),
#' \code{hsbeta} (Hessian, if available),
#' \code{dbeta_}, \code{pbeta_}, \code{qbeta_}, \code{rbeta_},
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{grad}} (for numerical gradient comparison),
#' \code{\link[base]{digamma}}.
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Gradient Evaluation
#'
#' # Generate sample data
#' set.seed(123)
#' n <- 1000
#' true_params <- c(gamma = 2.0, delta = 3.0)
#' data <- rbeta_(n, gamma = true_params[1], delta = true_params[2])
#'
#' # Evaluate gradient at true parameters
#' grad_true <- grbeta(par = true_params, data = data)
#' cat("Gradient at true parameters:\n")
#' print(grad_true)
#' cat("Norm:", sqrt(sum(grad_true^2)), "\n")
#'
#' # Evaluate at different parameter values
#' test_params <- rbind(
#'   c(1.5, 2.5),
#'   c(2.0, 3.0),
#'   c(2.5, 3.5)
#' )
#'
#' grad_norms <- apply(test_params, 1, function(p) {
#'   g <- grbeta(p, data)
#'   sqrt(sum(g^2))
#' })
#'
#' results <- data.frame(
#'   Gamma = test_params[, 1],
#'   Delta = test_params[, 2],
#'   Grad_Norm = grad_norms
#' )
#' print(results, digits = 4)
#'
#'
#' ## Example 2: Gradient in Optimization
#'
#' # Optimization with analytical gradient
#' fit_with_grad <- optim(
#'   par = c(1.5, 2.5),
#'   fn = llbeta,
#'   gr = grbeta,
#'   data = data,
#'   method = "L-BFGS-B",
#'   lower = c(0.01, 0.01),
#'   upper = c(100, 100),
#'   hessian = TRUE,
#'   control = list(trace = 0)
#' )
#'
#' # Optimization without gradient
#' fit_no_grad <- optim(
#'   par = c(1.5, 2.5),
#'   fn = llbeta,
#'   data = data,
#'   method = "L-BFGS-B",
#'   lower = c(0.01, 0.01),
#'   upper = c(100, 100),
#'   hessian = TRUE,
#'   control = list(trace = 0)
#' )
#'
#' comparison <- data.frame(
#'   Method = c("With Gradient", "Without Gradient"),
#'   Gamma = c(fit_with_grad$par[1], fit_no_grad$par[1]),
#'   Delta = c(fit_with_grad$par[2], fit_no_grad$par[2]),
#'   NegLogLik = c(fit_with_grad$value, fit_no_grad$value),
#'   Iterations = c(fit_with_grad$counts[1], fit_no_grad$counts[1])
#' )
#' print(comparison, digits = 4, row.names = FALSE)
#'
#'
#' ## Example 3: Verifying Gradient at MLE
#'
#' mle <- fit_with_grad$par
#' names(mle) <- c("gamma", "delta")
#'
#' # At MLE, gradient should be approximately zero
#' gradient_at_mle <- grbeta(par = mle, data = data)
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
#' grad_analytical <- grbeta(par = mle, data = data)
#' grad_numerical <- numerical_gradient(llbeta, mle, data)
#'
#' comparison_grad <- data.frame(
#'   Parameter = c("gamma", "delta"),
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
#' theta0 <- c(1.8, 2.8)
#' score_theta0 <- -grbeta(par = theta0, data = data)
#'
#' # Fisher information at theta0
#' fisher_info <- hsbeta(par = theta0, data = data)
#'
#' # Score test statistic
#' score_stat <- t(score_theta0) %*% solve(fisher_info) %*% score_theta0
#' p_value <- pchisq(score_stat, df = 2, lower.tail = FALSE)
#'
#' cat("\nScore Test:\n")
#' cat("H0: gamma=1.8, delta=2.8\n")
#' cat("Test statistic:", score_stat, "\n")
#' cat("P-value:", format.pval(p_value, digits = 4), "\n")
#'
#'
#' ## Example 6: Confidence Ellipse (Gamma vs Delta)
#'
#' # Observed information
#' obs_info <- hsbeta(par = mle, data = data)
#' vcov_full <- solve(obs_info)
#'
#' # Create confidence ellipse
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' eig_decomp <- eigen(vcov_full)
#' ellipse <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse[i, ] <- mle + sqrt(chi2_val) *
#'     (eig_decomp$vectors %*% diag(sqrt(eig_decomp$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_2d <- sqrt(diag(vcov_full))
#' ci_gamma <- mle[1] + c(-1, 1) * 1.96 * se_2d[1]
#' ci_delta <- mle[2] + c(-1, 1) * 1.96 * se_2d[2]
#'
#' # Plot
#'
#' plot(ellipse[, 1], ellipse[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "95% Confidence Region (Gamma vs Delta)", las = 1
#' )
#'
#' # Add marginal CIs
#' abline(v = ci_gamma, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_delta, col = "#808080", lty = 3, lwd = 1.5)
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
grbeta <- function(par, data) {
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
  .Call("_gkwdist_grbeta",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}


# ----------------------------------------------------------------------------#
# 7. HESSIAN (hsbeta)
# ----------------------------------------------------------------------------#

#' @title Hessian Matrix of the Negative Log-Likelihood for the Beta Distribution (gamma, delta+1 Parameterization)
#' @author Lopes, J. E.
#' @keywords distribution likelihood optimize hessian beta
#'
#' @description
#' Computes the analytic 2x2 Hessian matrix (matrix of second partial derivatives)
#' of the negative log-likelihood function for the standard Beta distribution,
#' using a parameterization common in generalized distribution families. The
#' distribution is parameterized by \code{gamma} (\eqn{\gamma}) and \code{delta}
#' (\eqn{\delta}), corresponding to the standard Beta distribution with shape
#' parameters \code{shape1 = gamma} and \code{shape2 = delta + 1}. The Hessian
#' is useful for estimating standard errors and in optimization algorithms.
#'
#' @param par A numeric vector of length 2 containing the distribution parameters
#'   in the order: \code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}).
#' @param data A numeric vector of observations. All values must be strictly
#'   between 0 and 1 (exclusive).
#'
#' @return Returns a 2x2 numeric matrix representing the Hessian matrix of the
#'   negative log-likelihood function, \eqn{-\partial^2 \ell / (\partial \theta_i \partial \theta_j)},
#'   where \eqn{\theta = (\gamma, \delta)}.
#'   Returns a 2x2 matrix populated with \code{NaN} if any parameter values are
#'   invalid according to their constraints, or if any value in \code{data} is
#'   not in the interval (0, 1).
#'
#' @details
#' This function calculates the analytic second partial derivatives of the
#' negative log-likelihood function (\eqn{-\ell(\theta|\mathbf{x})}) for a Beta
#' distribution with parameters \code{shape1 = gamma} (\eqn{\gamma}) and
#' \code{shape2 = delta + 1} (\eqn{\delta+1}). The components of the Hessian
#' matrix (\eqn{-\mathbf{H}(\theta)}) are:
#'
#' \deqn{
#' -\frac{\partial^2 \ell}{\partial \gamma^2} = n[\psi'(\gamma) - \psi'(\gamma+\delta+1)]
#' }
#' \deqn{
#' -\frac{\partial^2 \ell}{\partial \gamma \partial \delta} = -n\psi'(\gamma+\delta+1)
#' }
#' \deqn{
#' -\frac{\partial^2 \ell}{\partial \delta^2} = n[\psi'(\delta+1) - \psi'(\gamma+\delta+1)]
#' }
#'
#' where \eqn{\psi'(\cdot)} is the trigamma function (\code{\link[base]{trigamma}}).
#' These formulas represent the second derivatives of \eqn{-\ell(\theta)},
#' consistent with minimizing the negative log-likelihood. They correspond to
#' the relevant 2x2 submatrix of the general GKw Hessian (\code{\link{hsgkw}})
#' evaluated at \eqn{\alpha=1, \beta=1, \lambda=1}. Note the parameterization
#' difference from the standard Beta distribution (\code{shape2 = delta + 1}).
#'
#' The returned matrix is symmetric.
#'
#' @references
#' Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous Univariate
#' Distributions, Volume 2* (2nd ed.). Wiley.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
#' distributions. *Journal of Statistical Computation and Simulation*,
#'
#' (Note: Specific Hessian formulas might be derived or sourced from additional references).
#'
#' @seealso
#' \code{\link{hsgkw}}, \code{\link{hsmc}} (related Hessians),
#' \code{\link{llbeta}} (negative log-likelihood function),
#' \code{grbeta} (gradient, if available),
#' \code{dbeta_}, \code{pbeta_}, \code{qbeta_}, \code{rbeta_},
#' \code{\link[stats]{optim}},
#' \code{\link[numDeriv]{hessian}} (for numerical Hessian comparison),
#' \code{\link[base]{trigamma}}.
#'
#' @examples
#' \donttest{
#' ## Example 1: Basic Hessian Evaluation
#'
#' # Generate sample data
#' set.seed(123)
#' n <- 1000
#' true_params <- c(gamma = 2.0, delta = 3.0)
#' data <- rbeta_(n, gamma = true_params[1], delta = true_params[2])
#'
#' # Evaluate Hessian at true parameters
#' hess_true <- hsbeta(par = true_params, data = data)
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
#'   par = c(1.5, 2.5),
#'   fn = llbeta,
#'   gr = grbeta,
#'   data = data,
#'   method = "L-BFGS-B",
#'   lower = c(0.01, 0.01),
#'   upper = c(100, 100),
#'   hessian = TRUE
#' )
#'
#' mle <- fit$par
#' names(mle) <- c("gamma", "delta")
#'
#' # Hessian at MLE
#' hessian_at_mle <- hsbeta(par = mle, data = data)
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
#' names(se) <- c("gamma", "delta")
#'
#' # Correlation matrix
#' corr_matrix <- cov2cor(vcov_matrix)
#' cat("\nCorrelation Matrix:\n")
#' print(corr_matrix, digits = 4)
#'
#' # Confidence intervals
#' z_crit <- qnorm(0.975)
#' results <- data.frame(
#'   Parameter = c("gamma", "delta"),
#'   True = true_params,
#'   MLE = mle,
#'   SE = se,
#'   CI_Lower = mle - z_crit * se,
#'   CI_Upper = mle + z_crit * se
#' )
#' print(results, digits = 4)
#'
#' cat(sprintf(
#'   "\nMLE corresponds approx to Beta(%.2f, %.2f)\n",
#'   mle[1], mle[2] + 1
#' ))
#' cat(
#'   "True corresponds to Beta(%.2f, %.2f)\n",
#'   true_params[1], true_params[2] + 1
#' )
#'
#'
#' ## Example 4: Determinant and Trace Analysis
#'
#' # Compute at different points
#' test_params <- rbind(
#'   c(1.5, 2.5),
#'   c(2.0, 3.0),
#'   mle,
#'   c(2.5, 3.5)
#' )
#'
#' hess_properties <- data.frame(
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
#'   H <- hsbeta(par = test_params[i, ], data = data)
#'   eigs <- eigen(H, only.values = TRUE)$values
#'
#'   hess_properties <- rbind(hess_properties, data.frame(
#'     Gamma = test_params[i, 1],
#'     Delta = test_params[i, 2],
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
#' ## Example 5: Curvature Visualization (Gamma vs Delta)
#'
#' # Create grid around MLE
#' gamma_grid <- seq(mle[1] - 1.5, mle[1] + 1.5, length.out = 25)
#' delta_grid <- seq(mle[2] - 1.5, mle[2] + 1.5, length.out = 25)
#' gamma_grid <- gamma_grid[gamma_grid > 0]
#' delta_grid <- delta_grid[delta_grid > 0]
#'
#' # Compute curvature measures
#' determinant_surface <- matrix(NA,
#'   nrow = length(gamma_grid),
#'   ncol = length(delta_grid)
#' )
#' trace_surface <- matrix(NA,
#'   nrow = length(gamma_grid),
#'   ncol = length(delta_grid)
#' )
#'
#' for (i in seq_along(gamma_grid)) {
#'   for (j in seq_along(delta_grid)) {
#'     H <- hsbeta(c(gamma_grid[i], delta_grid[j]), data)
#'     determinant_surface[i, j] <- det(H)
#'     trace_surface[i, j] <- sum(diag(H))
#'   }
#' }
#'
#' # Plot
#'
#' contour(gamma_grid, delta_grid, determinant_surface,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Hessian Determinant", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#' contour(gamma_grid, delta_grid, trace_surface,
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "Hessian Trace", las = 1,
#'   col = "#2E4057", lwd = 1.5, nlevels = 15
#' )
#' points(mle[1], mle[2], pch = 19, col = "#8B0000", cex = 1.5)
#' points(true_params[1], true_params[2], pch = 17, col = "#006400", cex = 1.5)
#' grid(col = "gray90")
#'
#'
#' ## Example 6: Confidence Ellipse (Gamma vs Delta)
#'
#' # Extract 2x2 submatrix (full matrix in this case)
#' vcov_2d <- vcov_matrix
#'
#' # Create confidence ellipse
#' theta <- seq(0, 2 * pi, length.out = 100)
#' chi2_val <- qchisq(0.95, df = 2)
#'
#' eig_decomp <- eigen(vcov_2d)
#' ellipse <- matrix(NA, nrow = 100, ncol = 2)
#' for (i in 1:100) {
#'   v <- c(cos(theta[i]), sin(theta[i]))
#'   ellipse[i, ] <- mle + sqrt(chi2_val) *
#'     (eig_decomp$vectors %*% diag(sqrt(eig_decomp$values)) %*% v)
#' }
#'
#' # Marginal confidence intervals
#' se_2d <- sqrt(diag(vcov_2d))
#' ci_gamma <- mle[1] + c(-1, 1) * 1.96 * se_2d[1]
#' ci_delta <- mle[2] + c(-1, 1) * 1.96 * se_2d[2]
#'
#' # Plot
#'
#' plot(ellipse[, 1], ellipse[, 2],
#'   type = "l", lwd = 2, col = "#2E4057",
#'   xlab = expression(gamma), ylab = expression(delta),
#'   main = "95% Confidence Ellipse (Gamma vs Delta)", las = 1
#' )
#'
#' # Add marginal CIs
#' abline(v = ci_gamma, col = "#808080", lty = 3, lwd = 1.5)
#' abline(h = ci_delta, col = "#808080", lty = 3, lwd = 1.5)
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
hsbeta <- function(par, data) {
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
  .Call("_gkwdist_hsbeta",
    as.numeric(par),
    as.numeric(data),
    PACKAGE = "gkwdist"
  )
}
