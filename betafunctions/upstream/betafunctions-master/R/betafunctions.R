#' Compute Moments of Two-to-Four Parameter Beta Probability Density Distributions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Beta probability density distributions.
#' @param alpha The alpha shape parameter.
#' @param beta The beta shape parameter.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the first four raw, central, and standardized
#' # moments of this distribution using betamoments():
#' betamoments(alpha = 5, beta = 3, l = 0.25, u = 0.75,
#' types = c("raw", "central", "standardized"), orders = 4)
#' @references Hanson, B. A (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series.
#' @return A list of moment types, each a list of moment orders.
#' @export
betamoments <- function(alpha, beta, l = 0, u = 1, types = c("raw", "central", "standardized"), orders = 4) {
  moments <- base::list()
  if (any(types == "raw")) {
    for (i in 1:orders) {
      moments[["raw"]][[i]] <- stats::integrate(function(x) { dBeta.4P(x, l, u, alpha, beta) * x^i },
                                            lower = l, upper = u)$value
    }
  }
  if (any(types == "central")) {
    Mu <- stats::integrate(function(x) { dBeta.4P(x, l, u, alpha, beta) * x^1 }, lower = l, upper = u)$value
    for (i in 1:orders) {
      moments[["central"]][[i]] <- stats::integrate(function(x) { dBeta.4P(x, l, u, alpha, beta) * (x - Mu)^i },
                                            lower = l, upper = u)$value
    }
  }
  if (base::any(types == "standardized")) {
    Mu <- stats::integrate(function(x) { dBeta.4P(x, l, u, alpha, beta) * x^1 }, lower = l, upper = u)$value
    SigmaSquared <- stats::integrate(function(x) { dBeta.4P(x, l, u, alpha, beta) * (x - Mu)^2 },
                       lower = l, upper = u)$value
    for (i in 1:orders) {
      moments[["standardized"]][[i]] <- stats::integrate(function(x) { dBeta.4P(x, l, u, alpha, beta) * ((x - Mu)^i / sqrt(SigmaSquared)^i) },
                                            lower = l, upper = u)$value
    }
  }
  return(moments)
}

#' Compute Mode of Two- and Four-Parameter Beta Probability Density distribution.
#'
#' @description Computes the mode of a Beta distribution with specified shape- and location parameters.
#' @param alpha The alpha shape parameter of the Probability Density Distribution.
#' @param beta The beta shape parameter of the Probability Density Distribution.
#' @param l The first (lower) location parameter of a four-parameter distribution. Default set to \code{0}.
#' @param u The second (upper) location parameter of a four-parameter distribution. Default set to \code{1}.
#' @examples
#' # To calculate the mode of a two-parameter (standard) Beta distribution with
#' # shape parameters alpha = 5 and beta = 3:
#' betamode(alpha = 5, beta = 3)
#'
#' # To calculate the mode of a four-parameter Beta distribution with shape
#' # parameters alpha = 5 and beta = 3, and location parameters l = 25 and
#' # u = 150:
#' betamode(alpha = 5, beta = 3, l = 25, u = 150)
#' @export
betamode <- function(alpha, beta, l = 0, u = 1) {
  ((alpha - 1) / (alpha + beta - 2)) * (u - l) + l
}

#' Compute Median of Two- and Four-Parameter Beta Probability Density distribution.
#'
#' @description Computes the median of a Beta distribution with specified shape- and location parameters.
#' @param alpha The alpha shape parameter.
#' @param beta The beta shape parameter.
#' @param l The first (lower) location parameter. Default set to \code{0}.
#' @param u The second (upper) location parameter. Default set to \code{1}.
#' @examples
#' # To calculate the median of a two-parameter (standard) Beta distribution with
#' # shape parameters alpha = 5 and beta = 3:
#' betamedian(alpha = 5, beta = 3)
#'
#' # To calculate the median of a four-parameter Beta distribution with shape
#' # parameters alpha = 5 and beta = 3, and location parameters l = 25 and
#' # u = 150:
#' betamedian(alpha = 5, beta = 3, l = 25, u = 150)
#' @export
betamedian <- function(alpha, beta, l = 0, u = 1) {
  (alpha - (1/3)) / (alpha + beta - (2/3)) * (u - l) + l
}

#' Compute Moments of Binomial Probability Mass Functions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Binomial probability mass functions.
#' @param n Number of Binomial trials
#' @param p Probability of success per trial.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @examples
#' # Assume some variable follows a Binomial distribution with number of trials
#' # equal to 100 and a probability of success on each trial of 0.75. To compute
#' # the first four raw, central, and standardized moments of this distribution
#' # using binomialmoments():
#' binomialmoments(n = 100, p = 0.75, types = c("raw", "central",
#' "standardized"), orders = 4)
#'
#' # To only compute the (e.g.) standardized moments:
#' binomialmoments(n = 100, p = 0.75, types = "standardized")
#'
#' # To compute moments beyond the fourth order (e.g., the sixth):
#' binomialmoments(n = 100, p = 0.75, orders = 6)
#' @return A list of moment types, each a list of moment orders.
#' @export
binomialmoments <- function(n, p, types = c("raw", "central", "standardized"), orders = 4) {
  moments <- base::list()
  if (any(types == "raw")) {
    for (i in 1:orders) {
      moments[["raw"]][[i]] <- base::sum(stats::dbinom(0:n, n, p) * (0:n)^i)
    }
  }
  if (any(types == "central")) {
    Mu <- base::sum(stats::dbinom(0:n, n, p) * (0:n))
    for (i in 1:orders) {
      moments[["central"]][[i]] <- base::sum(stats::dbinom(0:n, n, p) * ((0:n - Mu)^i))
    }
  }
  if (base::any(types == "standardized")) {
    Mu <- base::sum(stats::dbinom(0:n, n, p) * (0:n))
    SigmaSquared <- base::sum(stats::dbinom(0:n, n, p) * ((0:n - Mu)^2))
    for (i in 1:orders) {
      moments[["standardized"]][[i]] <- base::sum(stats::dbinom(0:n, n, p) * ((0:n - Mu)^i / sqrt(SigmaSquared)^i))
    }
  }
  return(moments)
}

#' Compute Moments of Beta-Binomial Probability Mass Functions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Beta-Binomial probability mass functions.
#' @param N Number of trials.
#' @param l The first (lower) location-parameter of the Beta distribution.
#' @param u The second (upper) location-parameter of the Beta distribution.
#' @param alpha The alpha (first) shape-parameter of the Beta distribution.
#' @param beta The beta (second) shape-parameter of the Beta-distribution.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @examples
#' # Assume 100 observations of a discrete variable with probabilities of
#' # positive outcomes adhering to a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = .95, and shape parameters a = 5 and
#' # b = 3. To compute the first four raw, central, and standardized moments of
#' # this distrubution using betabinomialmoments():
#' betabinomialmoments(N = 100, l = .25, u = .95, alpha = 5, beta = 3,
#' types = c("raw", "central", "standardized"), orders = 4)
#' @references Hanson, B. A (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series.
#' @return A list of moment types, each a list of moment orders.
#' @export
betabinomialmoments <- function(N, l, u, alpha, beta, types = c("raw", "central", "standardized"), orders = 4) {
  if (l < 0 | u > 1) {
    base::warning("Beta-distribution location parameter(s) out of bounds (l < 0 or u > 1).")
  }
  weights <- NULL
  for (i in 0:N) {
    weights[i + 1] <- stats::integrate(function(x) { stats::dbinom(i, N, x) * dBeta.4P(x, l, u, alpha, beta) }, lower = 0, upper = 1)$value
  }
  moments <- base::list()
  if (base::any(types == "raw")) {
    for (i in 1:orders) {
      moments[["raw"]][[i]] <- sum((0:N)^i * weights)
    }
  }
  if (base::any(types == "central")) {
    mu1 <- sum((0:N)^1 * weights)
    for (i in 1:orders) {
      moments[["central"]][[i]] <- sum((0:N - mu1)^i * weights)
    }
  }
  if (base::any(types == "standardized")) {
    mu1 <- base::sum((0:N)^1 * weights)
    sigma2 <- base::sum((0:N - mu1)^2 * weights)
    for (i in 1:orders) {
      moments[["standardized"]][[i]] <- base::sum((0:N - mu1)^i * weights) / base::sqrt(sigma2)^i
    }
  }
  return(moments)
}

#' Compute Moments of Observed Value Distribution.
#'
#' @description Computes Raw, Central, or Standardized moment properties of a vector of observed scores.
#' @param x A vector of values, the distribution of which moments are to be calculated.
#' @param type A character vector determining which moment-types are to be calculated. Permissible values are \code{"raw"}, \code{"central"}, and \code{"standardized"}.
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @param correct Logical. Whether to include bias correction in estimation of orders. Default is \code{TRUE}.
#' @return A list of moment types, each a list of moment orders.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # To compute the first four raw, central, and standardized moments for this
#' # distribution of observed scores using observedmoments():
#' observedmoments(x = testdata, type = c("raw", "central", "standardized"),
#' orders = 4, correct = TRUE)
#' @export
observedmoments <- function(x, type = c("raw", "central", "standardized"),  orders = 4, correct = TRUE) {
  x <- stats::na.omit(x)
  moments <- base::list()
  if (base::any(type == "raw")) {
    for (i in 1:orders) {
      moments[["raw"]][[i]] <- base::sum(x^i)/base::length(x)
    }
  }
  if (base::any(type == "central")) {
    for (i in 1:orders) {
      if (correct) {
        moments[["central"]][[i]] <- base::sum((x - base::mean(x))^i)/(base::length(x) - 1)
      }
      else {
        moments[["central"]][[i]] <- base::sum((x - base::mean(x))^i)/(base::length(x))
      }
    }
  }
  if (base::any(type == "standardized")) {
    for (i in 1:orders) {
      if (correct) {
        moments[["standardized"]][[i]] <- 1 / base::length(x) * base::sum(((x - base::mean(x))^i)/sqrt(stats::var(x))^i)
      }
      else {
        moments[["standardized"]][[i]] <- 1 / base::length(x) * base::sum(((x - base::mean(x))^i)/sqrt(stats::var(x))^i)
      }
    }
  }
  return(moments)
}


#' Alpha Shape-Parameter Given Location-Parameters, Mean, and Variance a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the Beta value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the lower-location parameter will take on which ever value it must, and cannot be specified).
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param l The lower-bound location parameter of the Beta distribution. Default is 0 (as it is for the Standard Beta distribution).
#' @param u The upper-bound location parameter of the Beta distribution. Default is 1 (as it is for the Standard Beta distribution).
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Alpha shape-parameter in order to produce a  Beta probability density distribution with the target mean and variance, given specified lower- and upper bounds of the Beta distribution.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0, rescaled to proportion
#' # of maximum.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 1))
#'
#' # To find the alpha shape-parameter of a Standard (two-parameter) Beta
#' # distribution with the same mean and variance as the observed-score
#' # distribution using AMS():
#' AMS(mean(testdata), var(testdata))
#' @export
AMS <- function(mean, variance, l = 0, u = 1, sd = NULL) {
  if (!is.null(sd)) {
    variance <- sd^2
  }
  return(((l - mean) * (l * (mean - u) - mean^2 + mean * u - variance)) / (variance * (l - u)))
}

#' Beta Shape-Parameter Given Location-Parameters, Mean, and Variance of a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the Beta value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the lower-location parameter will take on which ever value it must, and cannot be specified).
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param l The lower-bound location parameter of the Beta distribution. Default is 0 (as it is for the Standard Beta distribution).
#' @param u The upper-bound location parameter of the Beta distribution. Default is 1 (as it is for the Standard Beta distribution).
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance, given specified lower- and upper bounds of the Beta distribution.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0, rescaled to proportion
#' # of maximum.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 1))
#'
#' # To find the beta shape-parameter of a Standard (two-parameter) Beta
#' # distribution with the same mean and variance as the observed-score
#' # distribution using BMS():
#' BMS(mean(testdata), var(testdata))
#'
#' # To find the beta shape-parameter of a four-parameter Beta
#' # distribution with specified lower- and upper-bounds of l = 0.25 and
#' # u = 0.75 using BMS:
#' BMS(mean(testdata), var(testdata), 0.25, 0.75)
#' @export
BMS <- function(mean, variance, l = 0, u = 1, sd = NULL) {
  if (!is.null(sd)) {
    variance <- sd^2
  }
  return(((mean - u) * (l * (u - mean) + mean^2 - mean * u + variance)) / (variance * (u - l)))
}

#' Lower Location Parameter Given Shape Parameters, Mean, Variance, and Upper Location Parameter of a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the lower-bound value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the lower-location parameter will take on which ever value it must, and cannot be specified).
#' @param alpha The alpha (first) shape-parameter of the target Beta probability density distribution.
#' @param beta The beta (second) shape-parameter of the target Beta probability density distribution.
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param u The upper-bound of the Beta distribution. Default is NULL (i.e., does not take a specified u-parameter into account).
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta lower location-parameter (\code{l}) in order to produce a Beta probability density distribution with the target moments and parameters.
#' @examples
#' # Generate some fictional data.
#' set.seed(1234)
#' testdata <- rBeta.4P(100000, 0.25, 0.75, 5, 3)
#' hist(testdata, xlim = c(0, 1), freq = FALSE)
#'
#' # Suppose you know three of the four necessary parameters to fit a four-
#' # parameter Beta distribution (i. e., u = 0.75, alpha = 5, beta = 3) to this
#' # data. To find the value for the necessary l parameter, estimate the mean
#' # and variance of the distribution:
#' M <- mean(testdata)
#' S2 <- var(testdata)
#'
#' # To find the l parameter necessary to produce a four-parameter Beta
#' # distribution with the target mean, variance, and u, alpha, and beta
#' # parameters using the LMSBAU() function:
#' (l <- LABMSU(alpha = 5, beta = 3, mean = M, variance = S2, u = 0.75))
#' curve(dBeta.4P(x, l, .75, 5, 3), add = TRUE, lwd = 2)
#' @export
LABMSU <- function(alpha = NULL, beta = NULL, u = NULL, mean = NULL, variance = NULL,  sd = NULL) {
  l <- NULL
  if (!is.null(sd)) {
    variance <- sd^2
  }
  if (is.null(l) & !is.null(alpha) & !is.null(beta) & !is.null(mean) & !is.null(u)) {
    l <- (alpha * mean - alpha * u + beta * mean) / beta
  }
  if (is.null(l) & is.null(u) & !is.null(alpha) & !is.null(beta) & !is.null(mean) & !is.null(variance)) {
    l <- mean - ((alpha * sqrt(variance * (alpha + beta + 1))) / sqrt(alpha * beta))
  }
  if (is.null(l) & !is.null(u) & !is.null(alpha) & !is.null(beta) & !is.null(mean) & !is.null(variance)) {
    l <- (alpha * (mean - u)^3 + beta * variance * (beta * u - mean + u)) / (beta^2 * variance)
  }
  return(l)
}

#' Upper Location Parameter Given Shape Parameters, Mean, Variance, and Lower Location Parameter of a Four-Parameter Beta Probability Density Distribution.
#'
#' @description Calculates the upper-bound value required to produce a Beta probability density distribution with defined moments and parameters. Be advised that not all combinations of moments and parameters can be satisfied (e.g., specifying mean, variance, skewness and kurtosis uniquely determines both location-parameters, meaning that the value of the upper-location parameter will take on which ever value it must, and cannot be specified).
#' @param alpha The alpha shape-parameter of the target Beta probability density distribution.
#' @param beta The beta shape-parameter of the target Beta probability density distribution.
#' @param mean The mean (first raw moment) of the target Standard Beta probability density distribution.
#' @param variance The variance (second central moment) of the target Standard Beta probability density distribution.
#' @param l The lower-bound of the Beta distribution. Default is NULL (i.e., does not take a specified l-parameter into account).
#' @param sd Optional alternative to specifying \code{var}. The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta upper location-parameter (\code{u}) in order to produce a Beta probability density distribution with the target moments and parameters.
#' @examples
#' # Generate some fictional data.
#' set.seed(1234)
#' testdata <- rBeta.4P(100000, 0.25, 0.75, 5, 3)
#' hist(testdata, xlim = c(0, 1), freq = FALSE)
#'
#' # Suppose you know three of the four necessary parameters to fit a four-
#' # parameter Beta distribution (i. e., l = 0.25, alpha = 5, beta = 3) to this
#' # data. To find the value for the necessary u parameter, estimate the mean
#' # and variance of the distribution:
#' M <- mean(testdata)
#' S2 <- var(testdata)
#'
#' # To find the l parameter necessary to produce a four-parameter Beta
#' # distribution with the target mean, variance, and u, alpha, and beta
#' # parameters using the LMSBAU() function:
#' (u <- UABMSL(alpha = 5, beta = 3, mean = M, variance = S2, l = 0.25))
#' curve(dBeta.4P(x, 0.25, u, 5, 3), add = TRUE, lwd = 2)
#' @export
UABMSL <- function(alpha = NULL, beta = NULL, mean = NULL, variance = NULL, l = NULL, sd = NULL) {
  u <- NULL
  if (!is.null(sd)) {
    variance <- sd^2
  }
  if (is.null(u) & !is.null(alpha) & !is.null(beta) & !is.null(mean) & !is.null(l)) {
    u <- (beta * (mean - l) / alpha) + mean
  }
  if (is.null(u) & is.null(l) & !is.null(alpha) & !is.null(beta) & !is.null(mean) & !is.null(variance)) {
    u <- mean + ((beta * sqrt(variance * (alpha + beta + 1))) / sqrt(alpha * beta))
  }
  if (is.null(u) & !is.null(l) & !is.null(alpha) & !is.null(beta) & !is.null(mean) & !is.null(variance)) {
    u <- (alpha * variance * (alpha * l + l - mean) - beta * (l - mean)^3) / (alpha^2 * variance)
  }
  if (is.null(u)) {
    warning("Not enough information.")
  }
  return(u)
}

#' Probability of Some Specific Observation under the Beta Probability Density Distribution with Specific Location Parameters, Mean, and Variance.
#'
#' @description Calculates the probability of some specific observation falling under a specified interval  ([0, x] or [x, 1]) under the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param q A specific point on the x-axis of the Standard Beta probability density distribution with a defined mean and variance.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param lower.tail Whether the density that should be considered is between the lower-end (i.e., [0 -> x]) or the higher-end of the distribution (i.e., [x -> 1]).
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return A value representing the probability of a random draw from the Standard Beta probability density distribution with a defined mean and variance being from one of two defined intervals (i.e., [0 -> x] or [x -> 1]).
#' @examples
#' # To compute the proportion of the density under the lower-end tail of a
#' # point along the Standard (two-parameter) Probability Density Distribution
#' # (e.g., 0.5) with mean of 0.6 and variance of 0.04:
#' pBetaMS(q = 0.5, mean = 0.6, variance = 0.04)
#'
#' # To compute the proportion of the density under the lower-end tail of a
#' # point along the Four-Parameter Beta Probability Density Distribution
#' # (e.g., 50) with mean of 60 and variance of 400, and lower-bound of 0 and
#' # upper-bound of 100:
#' pBetaMS(q = 50, mean = 60, variance = 400, l = 0, u = 100)
#' @export
pBetaMS <- function(q, mean, variance = NULL, sd = NULL, lower.tail = TRUE, l = 0, u = 1) {
  if ((!is.null(variance) & !is.null(sd))) {
    if (variance != sd^2) {
      warning("Nonequivalent values of variance and sd specified. Using variance.")
    }
  }
  if (base::is.null(variance) & !base::is.null(sd)) variance <- sd^2
  pBeta.4P(q, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u), lower.tail = lower.tail)
}

#' Density Under a Specific Point of the Beta Probability Density Distribution with Specific Location Parameters, Mean, and Variance.
#'
#' @description Calculates the density under specific points of the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param x A specific point on the x-axis of the Standard Beta Probability Density Distribution.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return A numeric value representing the required value for the beta Shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @examples
#' # To compute the density at a specific point (e.g., 0.5) along the Standard
#' # (two-parameter) Probability Density Distribution with mean of 0.6 and variance of 0.04:
#' dBetaMS(x = 0.5, mean = 0.6, variance = 0.04)
#'
#' # To compute the density at a specific point (e.g., 50) along the four-
#' # parameter Beta distribution with a mean of 60, variance of 400, and lower-
#' # bound of 0 and upper-bound of 100:
#' dBetaMS(x = 50, mean = 60, variance = 400, l = 0, u = 100)
#' @export
dBetaMS <- function(x, mean, variance = NULL, sd = NULL, l = 0, u = 1) {
  if ((!base::is.null(variance) & !base::is.null(sd))) {
    if (variance != sd^2) {
      warning("Nonequivalent values of variance and sd specified. Using variance.")
    }
  }
  if (base::is.null(variance) & !base::is.null(sd)) variance <- sd^2
  dBeta.4P(x, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u))
}

#' Quantile Containing Specific Proportion of the Distribution, Given a Specific Probability of the Beta Probability Density Distribution with Specific Mean and Variance.
#'
#' @description Calculates the quantile corresponding to a specific probability of some observation falling within the [0, x] (\code{lt = TRUE}) or [x, 1] (\code{lt = FALSE}) interval under the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param p A value of probability marking the point of the Y-axis to correspond to the X-axis.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param lower.tail Logical. Specifies which end of the tail for which to calculate quantile. Default is \code{TRUE} (meaning, find q for lower tail.)
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return A numeric value representing the quantile for which the specified proportion of observations fall within.
#' @examples
#' # To compute the quantile at a specific point (e.g., 0.5) along the Standard
#' # (two-parameter) Probability Density Distribution with mean of 0.6 and variance of 0.04:
#' qBetaMS(p = 0.5, mean = 0.6, variance = 0.04)
#'
#' # To compute the quantile at a specific points(e.g., 0.5) along the four-
#' # parameter Beta distribution with a mean of 60, variance of 400, and lower-
#' # bound of 0 and upper-bound of 100:
#' qBetaMS(p = 0.5, mean = 60, variance = 400, l = 0, u = 100)
#' @export
qBetaMS <- function(p, mean, variance = NULL, sd = NULL, lower.tail = TRUE, l = 0, u = 1) {
  if ((!base::is.null(variance) & !base::is.null(sd))) {
    if (variance != sd^2) {
      warning("Nonequivalent values of variance and sd specified. Using var.")
    }
  }
  if (base::is.null(variance) & !base::is.null(sd)) variance <- sd^2
  qBeta.4P(p, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u), lower.tail = lower.tail)
}

#' Random Draw from the Beta Probability Density Distribution With Specific Mean and Variance.
#'
#' @description Draws random samples of observations from the Standard Beta probability density distribution with defined mean and variance.
#' @param n Number of observations to be drawn from under the Standard Beta Probability Density Distribution.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param variance The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard probability density distribution.
#' @param l The lower-bound location parameter. Default set to 0 (the standard Beta distribution).
#' @param u The upper-bound location parameter. Default set to 1 (the standard Beta distribution).
#' @return A vector of length \code{n}, each value representing a random draw from the Standard Beta probability density distribution with defined mean and variance.
#' @examples
#' # To draw a random sample of 100 values from a Standard Beta distribution
#' # with a mean of 0.6 and variance = 0.04:
#' rBetaMS(n = 100, mean = 0.6, variance = 0.04)
#'
#' # To draw a random sample of 100 values from a Four-Parameter Beta
#* # distribution with a mean of 60 and variance = 400, and a lower-bound of
#' # 0 and an upper-bound of 100:
#' rBetaMS(n = 100, mean = 60, variance = 400, l = 0, u = 100)
#' @export
rBetaMS <- function(n, mean, variance = NULL, sd = NULL, l = 0, u = 1) {
  if ((!base::is.null(variance) & !base::is.null(sd))) {
    if (variance != sd^2) {
      warning("Nonequivalent values of var and sd specified. Using var.")
    }
  }
  if (is.null(variance) & !is.null(sd)) variance <- sd^2
  rBeta.4P(n, l, u, AMS(mean, variance, l, u), BMS(mean, variance, l, u))
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Probability Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Beta probability density functions.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-density values.
#' @param to The point of the x-axis to where y-density values are to be produced.
#' @param by The resolution (or spacing) at which to produce y-density values.
#' @param alpha The alpha (first) shape-parameter value for the Standard Beta probability density distribution.
#' @param beta The beta (second) shape-parameter for the Standard Beta probability density distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with density-values of y to plot against corresponding location values of x.
#' @examples
#' # To box in an area under a four-parameter Beta distribution with location
#' # parameters l = .25 and u = .75, and shape parameters alpha = 5 and
#' # rbeta = 3, from 0.4 to 0.6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 7))
#' coords <- Beta.gfx.poly.pdf(from = 0.4, to = 0.6, by = 0.001, alpha = 5,
#' beta = 3, l = 0.25, u = 0.75)
#' polygon(coords)
#' @export
Beta.gfx.poly.pdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  x <- base::c(from, base::seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- 0
    } else {
      y[i] <- dBeta.4P(x[i], l, u, alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Quantile Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Beta probability quantile functions.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-quantile values.
#' @param to The point of the \code{x}-axis to where \code{y}-quantile values are to be produced.
#' @param by The resolution (or spacing) at which to produce \code{y}-density values.
#' @param alpha The alpha shape-parameter value for the Standard Beta probability distribution.
#' @param beta The beta shape-parameter for the Standard Beta probability distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with quantile-values of \code{y} to plot against corresponding location values of \code{x}.
#' @examples
#' # To box in an area under a four-parameter Beta quantile distribution with
#' # location parameters l = .25 and u = 75, and shape parameters alpha = 5 and
#' # beta = 3, from .4 to .6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1))
#' coords <- Beta.gfx.poly.qdf(from = 0.4, to = 0.6, by = 0.001, alpha = 5,
#' beta = 3, l = 0.25, u = 0.75)
#' polygon(coords)
#' @export
Beta.gfx.poly.qdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  x <- base::c(from, base::seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- l
    } else {
      y[i] <- qBeta.4P(x[i], l, u, alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Cumulative Probability Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Beta cumulative probability functions.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-density values.
#' @param to The point of the \code{x}-axis to where \code{y}-density values are to be produced.
#' @param by The resolution (or spacing) at which to produce \code{y}-density values.
#' @param alpha The alpha shape-parameter value for the Standard Beta cumulative probability distribution.
#' @param beta The beta shape-parameter for the Standard Beta cumulative probability distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with cumulative probability-values of y to plot against corresponding location values of \code{x}.
#' @examples
#' # To box in an area under a four-parameter Beta cumulative distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters
#' # alpha = 5 and beta = 3, from 0.4 to 0.6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1))
#' coords <- Beta.gfx.poly.cdf(from = 0.4, to = 0.6, by = 0.001, alpha = 5,
#' beta = 3, l = 0.25, u = 0.75)
#' polygon(coords)
#' @export
Beta.gfx.poly.cdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  x <- base::c(from, seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- 0
    } else {
      y[i] <- pBeta.4P(x[i], l, u, alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Most Likely True Alpha Value Given Observed Outcome.
#'
#' @description Given a fitted Standard (two-parameter) Beta Distribution, return the alpha shape-parameter value where the observed mean becomes the mode.
#' @param alpha Observed alpha-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param beta Observed beta-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @return The Alpha shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-alpha parameter most likely to have produced the observations:
#' MLA(a = 10, b = 8)
#' @export
MLA <- function(alpha, beta, x = NULL, n = NULL) {
  if (base::is.null(x) | base::is.null(n)) {
    n <- alpha + beta
    A <- (alpha*(n - 2) + n) / n
    A
  } else {
    x*(n - 2) + 1
  }
}

#' Most Likely True Beta Value Given Observed Outcome.
#'
#' @description Assuming a prior standard (two-parameter) Beta Distribution, return the beta shape-parameter value where the observed mean becomes the mode.
#' @param alpha Observed alpha-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param beta Observed beta-parameter value for fitted Standard Beta Probability Density Distribution.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-beta parameter most likely to have produced the observations:
#' MLB(a = 10, b = 8)
#' @return The Beta shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @export
MLB <- function(alpha, beta, x = NULL, n = NULL) {
  if (base::is.null(x) | base::is.null(n)) {
    n <- alpha + beta
    A <- (alpha*(n - 2) + n) / n
    B <- n - A
    B
  } else {
    n - (x*(n - 2) + 1)
  }
}

#' Most Likely Mean of the Standard Beta Probability Density Distribution, Given that the Observation is Considered the Most Likely Observation of the Standard Beta Probability Density Distribution (i.e., the mode).
#'
#' @description Assuming a prior Standard (two-parameter) Beta Distribution, returns the expected mean of the distribution under the assumption that the observed value is the most likely value of the distribution.
#' @param alpha Observed alpha value for fitted Standard Beta Probability Density Distribution.
#' @param beta Observed beta value for fitted Standard Beta Probability Density Distribution.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-mean most likely to have produced the observations:
#' MLM(a = 10, b = 8)
#' @return The expected mean of the Standard Beta probability density distribution, for which the observed mean is the most likely value.
#' @export
MLM <- function(alpha, beta, x = NULL, n = NULL) {
  if (base::is.null(x) | base::is.null(n)) {
    n <- alpha + beta
    A <- (alpha*(n - 2) + n) / n
    B <- (n - A)
    A / (A + B)
  } else {
    (x*(n - 2) + 1) / n
  }
}

#' Probability Density under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Gives the density at desired values of \code{x} under the Four-Parameter Beta Probability Density Distribution.
#' @param x Value of \code{x}.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return The value for the probability density at specified values of \code{x}.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the probability density at a specific point of
#' # the distribution (e.g., 0.5) using dBeta.4P():
#' dBeta.4P(x = 0.5, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
dBeta.4P <- function(x, l, u, alpha, beta) {
 bfunc <- base::beta(alpha, beta)
  base::sapply(x, function(x) {
    if (x < l | x > u) {
      0
    } else {
      (1 / bfunc) * ((x - l)^(alpha - 1) * (u - x)^(beta - 1)) / (u - l)^(alpha + beta -1)
    }
    }
  )
}

#' Random Number Generation under the Beta-Binomial Probability Mass Distribution.
#'
#' @param n Number of draws.
#' @param N Number of trials.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The alpha (first) shape parameter.
#' @param beta The beta (second) shape parameter.
#' @return A vector with length \code{n} of random values drawn from the Beta-Binomial Distribution.
#' @examples
#' # To draw a sample of 50 values from a Beta-Binomial distribution with
#' # number of trials = 100, and with success-probabilities drawn from a
#' # Four-Parameter Beta distribution with location parameters l = 0.25 and
#' # u = 0.95, and shape-parameters alpha = 5 and beta = 3:
#' rBetaBinom(n = 50, N = 100, l = 0.25, u = 0.95, alpha = 5, beta = 3)
#' @export
rBetaBinom <- function(n, N, l, u, alpha, beta) {
  weights <- sapply(0:N, function(x) { stats::integrate(function(y) { stats::dbinom(x, N, y) * dBeta.4P(y, l, u, alpha, beta) }, lower = 0, upper = 1)$value})
  base::sample(0:N, size = n, prob = weights, replace = TRUE)
}

#' Probability Mass under the Beta-Binomial Probability-Mass Distribution.
#'
#' @description Gives the density at \code{x} under the Beta-Binomial PMF.
#' @param x Value of \code{x} (a specific number of successes).
#' @param N The total number of trials.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return The value for the probability mass at \code{x} given the specified Beta-Binomial distribution.
#' @examples
#' # Assume some variable follows a Beta-Binomial distribution with 100 number
#' # of trials, and with probabilities of successful trials drawn from a four-
#' # parameter Beta distribution with location parameters l = 0.25 and u = 0.75
#' # and shape parameters alpha = 5 and beta = 3. To compute the probability
#' # density at a specific point of the distribution (e.g., 50):
#' dBetaBinom(x = 50, N = 100, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
dBetaBinom <- function(x, N, l, u, alpha, beta) {
  sapply(x, function(x) { stats::integrate(function(y) { stats::dbinom(x, N, y) * dBeta.4P(y, l, u, alpha, beta) }, lower = 0, upper = 1)$value})
}

#' Cumulative Probability Function under the Beta-Binomial Probability Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under the Beta-Binomial Probability Distribution.
#' @param q The quantile or a vector of quantiles for which the proportion is to be calculated.
#' @param N The total number of trials.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lower.tail Whether the proportion to be calculated is to be under the lower or upper tail. Default is \code{TRUE} (lower tail).
#' @return A vector of proportions of observations falling under specified quantiles under the four-parameter Beta distribution.
#' @examples
#' # Assume some variable follows a Beta-Binomial distribution with number of
#' # trials = 50, and probabilities of successful trials are drawn from a four-
#' # parameter Beta distribution with location parameters l = 0.25 and u =
#' # 0.75, and shape parameters alpha = 5 and beta = 3. To compute the
#' # cumulative probability at a specific point of the distribution (e.g., 25):
#' pBetaBinom(q = 25, N = 50, l = .25, u = .75, alpha = 5, beta = 3)
#' @export
pBetaBinom <- function(q, N, l, u, alpha, beta, lower.tail = TRUE) {
  if (lower.tail) {
    1 - sapply(q, function(x) { sapply(x, function(x) { sum(dBetaBinom(x:N, N, l, u, alpha, beta)) })})
  } else {
    sapply(q, function(x) { sapply(x, function(x) { sum(dBetaBinom(x:N, N, l, u, alpha, beta)) })})
  }
}

#' Random Number Generation under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for generating random numbers from a specified Four-Parameter Beta Distribution.
#' @param n Number of draws.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The alpha (first) shape parameter.
#' @param beta The beta (second) shape parameter.
#' @return A vector with length \code{n} of random values drawn from the Four-Parameter Beta Distribution.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To draw a random value from this distribution using
#' # rBeta.4P():
#' rBeta.4P(n = 1, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
rBeta.4P <- function(n, l, u, alpha, beta) {
  stats::rbeta(n, alpha, beta) * (u - l) + l
}

#' Cumulative Probability Function under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under the Four-Parameter Beta Distribution.
#' @param q The quantile or a vector of quantiles for which the proportion is to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lower.tail Whether the proportion to be calculated is to be under the lower or upper tail. Default is \code{TRUE} (lower tail).
#' @return A vector of proportions of observations falling under specified quantiles under the four-parameter Beta distribution.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the cumulative probability at a specific point of
#' # the distribution (e.g., 0.5)
#' # using pBeta.4P():
#' pBeta.4P(q = 0.5, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
pBeta.4P <- function(q, l, u, alpha, beta, lower.tail = TRUE) {
  base::sapply(q, function(x) {
    num <- stats::integrate(function(y) { dBeta.4P(y, l, u, alpha, beta) }, lower = l, upper = x)$value
    den <- stats::integrate(function(y) { dBeta.4P(y, l, u, alpha, beta) }, lower = l, upper = u)$value
    if (lower.tail) {
      num/den
      } else {
        1 - num/den
      }
    }
  )
}

#' Quantile Given Probability Under the Four-Parameter Beta Distribution.
#'
#' @description Function for calculating the quantile (i.e., value of \code{x}) for a given proportion (i.e., the value of \code{y}) under the Four-Parameter Beta Distribution.
#' @param p A vector (or single value) of proportions or probabilities for which the corresponding value of \code{x} (i.e., the quantiles) are to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lower.tail Logical. Whether the quantile(s) to be calculated is to be under the lower or upper tail. Default is \code{TRUE} (lower tail).
#' @return A vector of quantiles for specified probabilities or proportions of observations under the four-parameter Beta distribution.
#' @examples
#' # Assume some variable follows a four-parameter Beta distribution with
#' # location parameters l = 0.25 and u = 0.75, and shape parameters alpha = 5
#' # and beta = 3. To compute the quantile at a specific point of the
#' # distribution (e.g., 0.5) using qBeta.4P():
#' qBeta.4P(p = 0.5, l = 0.25, u = 0.75, alpha = 5, beta = 3)
#' @export
qBeta.4P <- function(p, l, u, alpha, beta, lower.tail = TRUE) {
  if (lower.tail) {
    stats::qbeta(p, alpha, beta) * (u - l) + l
  } else {
    (1 - stats::qbeta(p, alpha, beta)) * (u - l) + l
  }
}

#' Method of Moment Estimates of Shape- and Location Parameters of the Four-Parameter Beta Distribution.
#'
#' @description An implementation of the method of moments estimation of four-parameter Beta distribution parameters presented by Hanson (1991). Given a vector of values, calculates the shape- and location parameters required to produce a four-parameter Beta distribution with the same mean, variance, skewness and kurtosis (i.e., the first four moments) as the observed-score distribution.
#' @param scores A vector of values to which the four-parameter Beta distribution is to be fitted.
#' @param mean If scores are not supplied: specification of the mean for the target four-parameter Beta distribution.
#' @param variance If scores are not supplied: specification of the variance for the target four-parameter Beta distribution.
#' @param skewness If scores are not supplied: specification of the skewness for the target four-parameter Beta distribution.
#' @param kurtosis If scores are not supplied: specification of the kurtosis for the target four-parameter Beta distribution.
#' @return A list of parameter-values required to produce a four-parameter Beta distribution with the same first four moments as the observed distribution.
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100), freq = FALSE)
#'
#' # To fit and retrieve the parameters for a four-parameter Beta distribution
#' # to the observed-score distribution using Beta.4p.fit():
#' (params.4p <- Beta.4p.fit(testdata))
#' curve(dBeta.4P(x, params.4p$l, params.4p$u, params.4p$alpha, params.4p$beta), add = TRUE)
#' @export
Beta.4p.fit <- function(scores, mean = NULL, variance = NULL, skewness = NULL, kurtosis = NULL) {
  if (!base::any(base::is.null(c(mean, variance, skewness, kurtosis)))) {
    m1 <- mean
    s2 <- variance
    g3 <- skewness
    g4 <- kurtosis
  } else {
    m <- observedmoments(scores)
    m1 <- m$raw[[1]]
    s2 <- m$central[[2]]
    g3 <- m$standardized[[3]]
    g4 <- m$standardized[[4]]
  }
  r <- 6 * (g4 - g3^2 - 1) / (6 + 3 * g3^2 - 2 * g4)
    if (g3 < 0) {
      a <- r / 2 * (1 + base::sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
      b <- r / 2 * (1 - base::sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
    } else {
      b <- r / 2 * (1 + base::sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
      a <- r / 2 * (1 - base::sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
    }
    l <- m1 - ((a * base::sqrt(s2 * (a + b + 1))) / base::sqrt(a * b))
    u <- m1 + ((b * base::sqrt(s2 * (a + b + 1))) / base::sqrt(a * b))
  return(base::list("alpha" = a, "beta" = b, "l" = l, "u" = u))
}

#' Method of Moment Estimates of Shape-Parameters of the Two-Parameter (Standard) Beta Distribution.
#'
#' @description An implementation of the method of moments estimation of two-parameter Beta distribution parameters. Given a vector of values, calculates the shape parameters required to produce a two-parameter Beta distribution with the same mean and variance (i.e., the first two moments) as the observed-score distribution.
#' @param scores A vector of values to which the two-parameter Beta distribution is to be fitted. The values ought to fall within the [0, 1] interval.
#' @param mean The mean of the target Beta distribution. Alternative to feeding the function raw scores.
#' @param variance The variance of the target Beta distribution. Alternative to feeding the function raw scores.
#' @param l Optional specification of a lower-bound parameter of the Beta distribution. Default is 0 (i.e., the lower-bound of the Standard two-parameter Beta distribution).
#' @param u Optional specification of an upper-bound parameter of the Beta distribution. Default is 1 (i.e., the lower-bound of the Standard two-parameter Beta distribution).
#' @return A list of parameter-values required to produce a Standard two-parameter Beta distribution with the same first two moments as the observed distribution.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 1), freq = FALSE)
#'
#' # To fit and retrieve the parameters for a two-parameter Beta distribution
#' # to the observed-score distribution using Beta.2p.fit():
#' (params.2p <- Beta.2p.fit(testdata))
#' curve(dbeta(x, params.2p$alpha, params.2p$beta), add = TRUE)
#' @export
Beta.2p.fit <- function(scores = NULL, mean = NULL, variance = NULL, l = 0, u = 1) {
  if (!is.null(scores)) {
    if (base::max(scores) > u | base::min(scores) < l) {
      warning(paste("Input values outside the range of the specified location parameters of the Beta distribution (i.e., there are values falling outside the [", l, ", ", u, "] interval).", sep = ""))
    }
    m1 <- base::mean(scores)
    s2 <- stats::var(scores)
  } else {
    m1 <- mean
    s2 <- variance
  }
  a <- AMS(m1, s2, l, u)
  b <- BMS(m1, s2, l, u)
  return(base::list("alpha" = a, "beta" = b, "l" = l, "u" = u))
}

#' An implementation of the Beta-density Compound Cumulative Binomial Distribution.
#'
#' @description The Beta Compound Binomial distribution: The product of the four-parameter Beta probability density function and the binomial cumulative probability mass function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param alpha The alpha shape-parameter of the Beta distribution.
#' @param beta The beta shape-parameter of the Beta distribution.
#' @param n The number of trials for the Binomial distribution.
#' @param c The "true-cut" (proportion) of the Binomial distribution.
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Binomial distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @note The Binomial distribution cut-point is up-to but not including, unlike the standard behaviour of base-R pbinom() function.
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Binomial error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0 can be calculated as:
#' dBeta.pBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5,
#' lower.tail = TRUE)
#'
#' #By integration, the population proportion of (e.g.) passing scores in some
#' #region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' #calculated as:
#' integrate(function(x) { dBeta.pBinom(x, 0.25, .75, 5, 3, 10, 0.5) },
#' lower = 0.25, upper = 0.5)
#' @export
dBeta.pBinom <- function(x, l, u, alpha, beta, n, c, lower.tail = FALSE) {
  if (!lower.tail) {
    dBeta.4P(x, l, u, alpha, beta) * (1 - stats::pbinom(round(n * c) - 1, round(n), x, lower.tail = TRUE))
  } else {
    dBeta.4P(x, l, u, alpha, beta) * stats::pbinom(round(n * c) - 1, round(n), x, lower.tail = TRUE)
  }
}

#' An implementation of a Beta-density Compound Cumulative Gamma-Binomial Distribution.
#'
#' @description The Beta Compound Binomial distribution: The product of the four-parameter Beta probability density function and the binomial cumulative probability mass function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param alpha The alpha shape-parameter of the four-parameter Beta distribution.
#' @param beta The beta shape-parameter of the four-parameter Beta distribution.
#' @param n The number of "trials" for the Gamma-Binomial distribution.
#' @param c The "true-cut" (proportion) on the Gamma-Binomial distribution. Need not be an integer (unlike Binomial distribution).
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Binomial distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @export
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Binomial error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0 can be calculated as:
#' dBeta.pGammaBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pGammaBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10.1, c = 0.5,
#' lower.tail = TRUE)
#'
#' #By integration, the population proportion of (e.g.) passing scores in some
#' #region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' #calculated as:
#' integrate(function(x) { dBeta.pGammaBinom(x, 0.25, 0.75, 5, 3, 10, 0.5) },
#' lower = 0.25, upper = 0.5)
dBeta.pGammaBinom <- function(x, l, u, alpha, beta, n, c, lower.tail = FALSE) {
  if (!lower.tail) {
    dBeta.4P(x, l, u, alpha, beta) * pGammaBinom(n * c, n, x, FALSE)
  } else {
    dBeta.4P(x, l, u, alpha, beta) * pGammaBinom(n * c, n, x, lower.tail = TRUE)
  }
}

#' An implementation of a Beta-density Compound Cumulative-Beta Distribution.
#'
#' @description The Beta Compound Beta distribution: The product of the four-parameter Beta probability density function and the Beta cumulative probability function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param alpha The alpha shape-parameter of the Beta density distribution.
#' @param beta The beta shape-parameter of the Beta density distribution.
#' @param n The number of trials for the Beta cumulative probability distribution.
#' @param c The "true-cut" (proportion) of on the Beta cumulative probability distribution.
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Beta cumulative probability distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Beta error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0.5 can be calculated as:
#' dBeta.pBeta(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pBeta(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5,
#' lower.tail = TRUE)
#'
#' # By integration, the population proportion of (e.g.) passing scores in some
#' # region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' # calculated as:
#' integrate(function(x) { dBeta.pBeta(x, 0.25, 0.75, 5, 3, 10, 0.5) },
#' lower = 0.25, upper = 0.5)
#' @export
dBeta.pBeta <- function(x, l, u, alpha, beta, n, c, lower.tail = FALSE) {
  if(!lower.tail) {
    dBeta.4P(x, l, u, alpha, beta) * stats::pbeta(c, x * n, (1 - x) * n, lower.tail = FALSE)
  } else {
    dBeta.4P(x, l, u, alpha, beta) * (1 - stats::pbeta(c, x * n, (1 - x) * n, lower.tail = FALSE))
  }
}

#' Gamma-extended Binomial coefficient (choose function).
#'
#' @description Extends the Binomial coefficient for positive non-integers (including 0) by employing the Gamma rather than the factorial function.
#' @param n In Binomial terms, the number of Binomial "trials". Need not be an integer.
#' @param k In Binomial terms, the number of successful "trials". Need not be an integer.
#' @note Not defined for negative integers.
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' # Compare choose function with gchoose function for integers:
#' gchoose(c(8, 9, 10), c(3, 4, 5)) == choose(c(8, 9, 10), c(3, 4, 5))
#'
#' # The gchoose function also works for non-integers:
#' gchoose(10.5, 7.5)
#' @export
gchoose <- function(n, k) {
  gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))
}

#' Cumulative probability density function under the Gamma-extended Binomial distribution.
#'
#' @description Extends the cumulative Binomial probability mass function to positive non-integers, effectively turning the mass-function into a density-function.
#' @param q Vector of quantiles.
#' @param size Number of "trials" (zero or more). Need not be integer.
#' @param prob Probability of "success" on each "trial". Need not be integer.
#' @param lower.tail Logical. If TRUE (default), probabilities are P[X<x], otherwise, P[X >= x]. Note that this differs from base-R \code{binom()} functions.
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' # Assume some variable follows a Gamma-Binomial  distribution with
#' # "number of trials" = 10.5 and probability of "success" for each "trial"
#' # = 0.75, to compute the cumulative probability to attain a "number of
#' # success" below a specific point (e.g., less than 7.5 "successes":
#' pGammaBinom(q = 7.5, size = 10.5, prob = 0.75)
#'
#' # Conversely, to attain a value at or above 7.5:
#' pGammaBinom(q = 7.5, size = 10.5, prob = 0.75, lower.tail = FALSE)
#' @export
pGammaBinom <- function(q, size, prob, lower.tail = TRUE) {
  base::sapply(prob, function(x) {
    num <- stats::integrate(function(y) { dGammaBinom(y, size, x) }, lower = 0, upper = q)$value
    den <- stats::integrate(function(y) { dGammaBinom(y, size, x) }, lower = 0, upper = size)$value
    if (lower.tail) {
      num/den
    } else {
      1 - num/den
    }
  })
}

#' Probability density function under the Gamma-extended Binomial distribution.
#'
#' @param x Vector of quantiles.
#' @param size Number of "trials" (zero or more). Need not be integer.
#' @param prob Probability of "success" on each "trial". Need not be integer.
#' @param nc Whether to include a normalizing constant making sure that the sum of the distribution's density is 1.
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' #' # Assume some variable follows a Gamma-Binomial distribution with
#' # "number of trials" = 10.5 and probability of "success" for each "trial"
#' # = 0.75, to compute the probability density to attain a "number of success"
#' # at a specific point (e.g., 7.5 "successes"):
#' dGammaBinom(x = 7.5, size = 10.5, prob = 0.75)
#'
#' # Including a normalizing constant (then diverges from binomial dist.):
#' dGammaBinom(x = 7.5, size = 10.5, prob = 0.75, nc = TRUE)
#' dGammaBinom(x = 7, size = 10, prob = 0.75) == dbinom(7, 10, 0.75)
#' dGammaBinom(x = 7, size = 10, prob = 0.75, nc = TRUE) == dbinom(7, 10, 0.75)
#' @export
dGammaBinom <- function(x, size, prob, nc = FALSE) {
  if (nc) {
    den <- stats::integrate(function(z) { (gchoose(size, z) * prob^z * (1 - prob)^(size - z)) }, lower = 0, upper = size)$value
    base::sapply(x, function(y) {
      if (y < 0 | y > size) {
        0
      } else {
        (gchoose(size, y) * prob^y * (1 - prob)^(size - y)) / den
      }
    })
  } else {
    base::sapply(x, function(y) {
      if (y < 0 | y > size) {
        0
      } else {
        (gchoose(size, y) * prob^y * (1 - prob)^(size - y))
      }
    })
  }
}

#' Random number generation under the Gamma-extended Binomial distribution.
#'
#' @param n Number of observations.
#' @param size Number of "trials" (zero or more). Need not be integer.
#' @param prob Probability of "success" on each "trial". Need not be integer.
#' @param precision The precision with which the quantile is to be calculated. Default is 1e-4 (i.e., search terminates when there is no registered change in estimate at the fourth decimal). Tuning this value will impact the time it takes for the search algorithm to arrive at an estimate.
#' @note Calls \code{qGammaBinom()}, which makes the random draw slower than what one might be used to (since \code{qGammaBinom()} calls \code{pGammaBinom()} and employs a search-algorithm to find the appropriate value down to a specifiable level of precision).
#' @examples
#' # Assume some variable follows a Gamma-Binomial distribution with
#' # "number of trials" = 10.5 and probability of "success" for each "trial"
#' # = 0.75 To draw a random value from this distribution:
#' rGammaBinom(n = 1, size = 10, prob = 0.75)
#' @export
rGammaBinom <- function(n, size, prob, precision = 1e-4) {
  qGammaBinom(stats::runif(n, 0, 1), size = size, prob = prob, precision = precision)
}

#' Quantile function for the Gamma-extended Binomial distribution.
#'
#' @param p Vector of probabilities.
#' @param size Number of "trials" (zero or more, including positive non-integers).
#' @param prob Probability of success on each "trial".
#' @param lower.tail Logical. If TRUE (default), probabilities are P[X < x], otherwise P[X > x].
#' @param precision The precision with which the quantile is to be calculated. Default is 1e-7 (i.e., search terminates when there is no registered change in estimate at the seventh decimal). Tuning this value will impact the time it takes for the search algorithm to arrive at an estimate.
#' @note This function uses a bisection search-algorithm to find the number of successes corresponding to the specified quantile(s). This algorithm is inefficient with respect to the number of iterations required to converge on the solution. More efficient algorithms might be added in later versions.
#' @references Loeb, D. E. (1992). A generalization of the binomial coefficients. Discrete Mathematics, 105(1-3).
#' @examples
#' # For a Gamma-extended Binomial distribution with number of trials = 10 and
#' # probability of success per trial of 0.75, calculate the number of success-
#' # ful trials at or below the 25% quantile:
#' qGammaBinom(p = 0.25, size = 10, prob = 0.75)
#'
#' # Conversely, for a Gamma-extended Binomial distribution with number of
#' # trials = 10 and probability of success per trial of 0.75, calculate the
#' # number of successful trials at or above the 25% quantile:
#' qGammaBinom(p = 0.25, size = 10, prob = 0.75, lower.tail = FALSE)
#' @export
qGammaBinom <- function(p, size, prob, lower.tail = TRUE, precision = 1e-7) {
  base::sapply(p, function(c) {
    if (c == 0 | c == 1) {
      if (c == 0) {
        x <- 0
      } else {
        x <- size
      }
      x
    } else {
      x <- size / 2
      y <- pGammaBinom(x, size, prob, lower.tail)
      a <- size
      b <- x
      while(base::abs(y - c) > precision) {
        if (y < c) {
          if (lower.tail) {
            x <- x + abs((a - b)) / 2
          } else {
            x <- x - abs((a - b)) / 2
          }
        } else {
          if (lower.tail) {
            x <- x - base::abs((a - b)) / 2
          } else {
            x <- x + base::abs((a - b)) / 2
          }
        }
        a <- b
        b <- x
        y <- pGammaBinom(x, size, prob, lower.tail)
      }
    }
    x
    })
}


#' Probability Mass function for Lord's Two-Term Approximation to the Compound Binomial Distribution.
#'
#' @description Gives the density at \code{x} under Lord's two-term approximation to the compound Binomial PMF.
#' @param x Value of \code{x} (a specific number of successes).
#' @param N The total number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param p Probability of success for each trial.
#' @export
#' @examples
#' # Assume some variable follows a compound Binomial distribution with 100
#' # trials, a 50% probability of success on each trial, and Lord's k = 1. To
#' # compute the probability density at a specific point of the distribution
#' # (e.g., 50):
#' dcBinom(x = 50, N = 100, k = 1, p = .5)
dcBinom <- function(x, N, k, p) {
  sapply(x, function(x) {
    stats::dbinom(x, N, p) - k*p*(1 - p)*(stats::dbinom(x, N - 2, p) - 2*stats::dbinom(x - 1, N - 2, p) + stats::dbinom(x - 2, N - 2, p))
  })
}

#' Cumulative Probability Mass function for Lord's Two-Term Approximation to the Compound Binomial Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under Lord's two-term approximation to the compound Binomial distribution.
#' @param q The quantile or vector of quantiles for which the proportion is to be calculated.
#' @param N Total number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param p Probability of success for each trial.
#' @param lower.tail Logical. If TRUE (default), probabilities are P[X<x], otherwise, P[X >= x]. Note that this differs from base-R \code{binom()} functions.
#' @export
#' @examples
#' # Assume some variable follows a compound Binomial distribution with 100
#' # trials, a 50% probability of success on each trial, and Lord's k = 1. To
#' # compute the cumulative probability at a specific point of the distribution
#' # (e.g., 50):
#' pcBinom(q = 50, N = 100, k = 1, p = .5)
pcBinom <- function(q, N, k, p, lower.tail = TRUE) {
  if (lower.tail) {
    1 - sapply(q, function(x) {
      sapply(x, function(x) {
        sum(dcBinom(x:N, N, k, p))
      })
    })
  } else {
    sapply(q, function(x) {
      sapply(x, function(x) {
        sum(dcBinom(x:N, N, k, p))
      })
    })
  }
}

#' Random Number Generation under Lord's Two-Term Approximation to the Compound Binomial Distribution.
#'
#' @description Random Number Generation under Lord's Two-Term Approximation to the Compound Binomial Distribution.
#' @param n Number of draws.
#' @param N Number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param p Probability of success for each trial.
#' @note For larger values of \code{k}, the distribution can yield negative probabilities. This function handles such occurrences by adding the absolute value of the minimum probability to all observations if there are any negative probabilities and then normalize the distribution so that the total density is equal to 1.
#' @export
#' @examples
#' # To draw a sample of 50 values from a Compound-Binomial distribution with
#' # number of trials = 100, a 50% probability of success for each trial, and
#' # Lord's k = 1:
#' set.seed(1234)
#' rcBinom(n = 50, N = 100, k = 1, p = .5)
#'
#' # To draw values where the probabilities vary for each draw:
#' rcBinom(n = 50, N = 100, k = 1, p = runif(50))
rcBinom <- function(n, N, k, p) {
  if (length(p) == 1) {
    weights <- dcBinom(0:N, N, k, p)
    if (any(weights < 0)) {
      warning("Density function returned negative probabilities.")
      weights <- weights + abs(min(weights))
    }
    if (sum(weights) != 1) {
      warning("Sum of probabilities did not equal 1.")
      weights <- weights/sum(weights)
    }
    sample(0:N, size = n, prob = weights, replace = TRUE)
  } else {
    if (length(p) < n) {
      p <- as.vector(matrix(p, nrow = n, ncol = 1))
    }
    sapply(p, function(x) {
      weights <- dcBinom(0:N, N, k, x)
      if (any(weights < 0)) {
        warning("Density function returned negative probabilities.")
        weights <- weights + abs(min(weights))
      }
      if (sum(weights != 1)) {
        warning("Sum of probabilities did not equal 1.")
        weights <- weights/sum(weights)
      }
      sample(0:N, size = 1, prob = weights, replace = TRUE)
    })
  }
}

#' Probability Mass function for Lord's Beta Compound Binomial Distribution.
#'
#' @description Gives the density at \code{x} under the Beta Compound-Binomial distribution, where the Compound-Binomial distribution is Lord's two-term approximation.
#' @param x Value of \code{x} (a specific number of successes).
#' @param N Number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param l The lower-bound location parameter of the four-parameter Beta distribution.
#' @param u The upper-bound location parameter of the four-parameter Beta distribution.
#' @param alpha The first shape-parameter of the four-parameter Beta distribution.
#' @param beta The second shape-parameter of the four-parameter Beta distribution.
#' # Assume some variable follows a Beta compound Binomial distribution with 100
#' # trials, Lord's k = 1, and probabilities of successful trials drawn from a
#' # four-parameter Beta distribution with location-parameters l = .15 and u =
#' # .85, and shape parameters alpha = 6 and beta = 4. To compute the
#' # probability density at a specific point of the distribution (e.g., 50):
#' dBetacBinom(x = 50, N = 100, k = 1, l = .15, u = .85, alpha = 6, beta = 4)
dBetacBinom <- function(x, N, k, l, u, alpha, beta) {
  sapply(x, function(x) {
    stats::integrate(function(y) {
      dcBinom(x, N, k, y) * dBeta.4P(y, l, u, alpha, beta)
    }, lower = 0, upper = 1)$value
  })
}

#' Random Number Generation under Lord's Beta Compound-Binomial Distribution.
#'
#' @description Random number generation under Lord's Beta Compound-Binomial distribution, where the Compound-Binomial distribution is Lord's two-term approximation.
#' @param x Number of draws.
#' @param N Number of trials.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param l The lower-bound location parameter of the four-parameter Beta distribution.
#' @param u The upper-bound location parameter of the four-parameter Beta distribution.
#' @param alpha The first shape-parameter of the four-parameter Beta distribution.
#' @param beta The second shape-parameter of the four-parameter Beta distribution.
#' @note For larger values of \code{k}, the distribution can yield negative probabilities which returns an error.
#' @export
#' @examples
#' # To draw a sample of 50 values from a Beta Compound-Binomial distribution
#' # with number of trials = 100, Lord's k = 1, and probabilities of successful
#' # trials drawn from a four-parameter Beta distribution with location-
#' # parameters l = .15 and u = .85, and shape parameters alpha = 6 and
#' # beta = 4:
#' rBetacBinom(x = 50, N = 100, k = 1, l = .15, u = .85, alpha = 6, beta = 4)
rBetacBinom <- function(x, N, k, l, u, alpha, beta) {
  weights <- dBetacBinom(0:N, N, k, l, u, alpha, beta)
  if (any(weights < 0)) {
    warning("Density function returned negative probabilities.")
    weights <- weights + abs(min(weights))
  }
  if (sum(weights != 1)) {
    warning("Sum of probabilities did not equal 1.")
    weights <- weights/sum(weights)
  }
  sample(0:N, size = x, prob = weights, replace = TRUE)
}
