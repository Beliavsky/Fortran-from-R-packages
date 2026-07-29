#' Data-driven Local Polynomial for the Trend's Derivatives in Equidistant Time
#' Series
#'
#'  This function runs through an iterative process in order to find the
#'  optimal bandwidth for the nonparametric estimation of the first or second
#'  derivative of the trend in an equidistant time series (with short-memory
#'  errors) and subsequently employs the obtained bandwidth via local
#'  polynomial regression.
#'
#'@param y a numeric vector that contains the time series ordered from past to
#'present.
#'@param d an integer \code{1} or \code{2} that defines the order of
#'derivative; the default is \code{d = 1}.
#'@param mu an integer \code{0}, ..., \code{3} that represents the smoothness
#'parameter of the kernel weighting function and thus defines the kernel
#'function that will be used within the local polynomial regression; is set to
#'\code{1} by default.
#'
#'\tabular{cl}{
#'\strong{Number} \tab \strong{Kernel}\cr
#'\code{0} \tab Uniform Kernel\cr
#'\code{1} \tab Epanechnikov Kernel\cr
#'\code{2} \tab Bisquare Kernel\cr
#'\code{3} \tab Triweight Kernel
#'}
#'@param pp an integer \code{1} (local linear regression) or \code{3} (local
#'cubic regression) that indicates the order of polynomial upon which
#'\eqn{c_f}, i.e. the variance factor, will be calculated by
#'\code{\link{msmooth}}; the default is \code{pp = 1}.
#'@param bStart.p a numeric object that indicates the starting value of the
#'bandwidth for the iterative process for the calculation of \eqn{c_f}; should
#'be \eqn{> 0}; is set to \code{0.15} by default.
#'@param bStart a numeric object that indicates the starting value of the
#'bandwidth for the iterative process; should be \eqn{> 0}; is set to
#'\code{0.15} by default.
#'
#'@details
#'The trend's derivative is estimated based on the additive
#'nonparametric regression model for an equidistant time series
#'\deqn{y_t = m(x_t) + \epsilon_t,}
#'where \eqn{y_t} is the observed time series, \eqn{x_t} is the rescaled time
#'on the interval \eqn{[0, 1]}, \eqn{m(x_t)} is a smooth and deterministic
#'trend function and \eqn{\epsilon_t} are stationary errors with
#'\eqn{E(\epsilon_t) = 0} and short-range dependence (see also Beran and Feng,
#'2002). With this function, the first or second derivative of \eqn{m(x_t)}
#'can be estimated without a parametric model assumption for the error series.
#'
#'The iterative-plug-in (IPI) algorithm, which numerically minimizes the
#'Asymptotic Mean Squared Error (AMISE), was proposed by Feng, Gries and
#'Fritz (2020).
#'
#'Define \eqn{I[m^{(k)}] = \int_{c_b}^{d_b} [m^{(k)}(x)]^2 dx}{I[m^(k)] =
#'int_[c_b]^[d_b] \{m^(k)\}^2 dx}, \eqn{\beta_{(\nu, k)} = \int_{-1}^{1} u^k
#'K_{(\nu, k)}(u) du}{\beta_(\nu, k) = int_[-1]^[1] u^k K_(\nu, k)(u) du}
#'and \eqn{R(K) = \int_{-1}^{1} K_{(\nu, k)}^{2}(u) du}{R(K) = int_[-1]^[1]
#'\{K_(\nu, k)\}^2 du}, where \eqn{p} is the order of the polynomial,
#'\eqn{k = p + 1} is the order of the asymptotically equivalent kernel,
#'\eqn{\nu} is the order of the trend function's derivative, \eqn{0 \leq c_{b}
#'< d_{b} \leq 1}{0 \le c_b < d_b \le 1}, \eqn{c_f} is the variance factor and
#'\eqn{K_{(\nu, k)}(u)}{K_(\nu, k)(u)} the \eqn{k}-th order equivalent kernel
#'obtained for the estimation of \eqn{m^{(\nu)}}{m^(\nu)} in the interior.
#'\eqn{m^{(\nu)}}{m^(\nu)} is the \eqn{\nu}-th order derivative (\eqn{\nu = 0,
#'1, 2, ...}) of the nonparametric trend.
#'
#'Furthermore, we define
#'\deqn{C_{1} = \frac{I[m^{(k)}] \beta_{(\nu, k)}^2}{(k!)^2}}{C_1 =
#'[I[m^(k)]\{\beta_(\nu, k)\}^2] / (k!)^2}
#'and
#'\deqn{C_{2} = \frac{2 \pi c_{f} (d_b - c_b) R(K)}{nh^{2 \nu + 1}}}{C_2 =
#'2\pi(d_b - c_b)R(K)c_f / (nh^[2\nu + 1])}
#'with \eqn{h} being the bandwidth and \eqn{n} being the number of
#'observations. The AMISE is then
#'\deqn{AMISE(h) = h^{2(k-\nu)}C_{1} + C_{2}.}{AMISE(h) = h^[2(k - \nu)] C_1 +
#'C_2.}
#'
#'The variance factor \eqn{c_f} is first obtained from a pilot-estimation of
#'the time series' nonparametric trend (\eqn{\nu = 0}) with polynomial order
#'\eqn{p_p}. The estimate is then plugged into the iterative procedure for
#'estimating the first or second derivative (\eqn{\nu = 1} or \eqn{\nu = 2}).
#'For further details on the asymptotic theory or the algorithm, we refer the
#'user to Feng, Fritz and Gries (2020) and Feng et al. (2019).
#'
#'The function itself is applicable in the following way: Based on a data input
#'\code{y}, an order of polynomial \code{pp} for the variance factor estimation
#'procedure, a starting value for the relative bandwidth \code{bStart.p} in the
#'variance factor estimation procedure, a kernel function defined by the
#'smoothness parameter \code{mu} and a starting value for the relative
#'bandwidth \code{bStart} in the bandwidth estimation procedure, an optimal
#'bandwidth is numerically calculated for the trend's derivative of order
#'\code{d}. In fact, aside from the input vector \code{y}, every argument has a
#'default setting that can be adjusted for the individual case. However, it is
#'recommended to initially use the default values for the estimation of the
#'first derivative and adjust the argument \code{d} to \code{d = 2} for the
#'estimation of the second derivative. Following Feng, Gries and Fritz (2020),
#'the initial bandwidth does not affect the resulting optimal bandwidth in
#'theory. However in practice, local minima of the AMISE can influence the
#'results. Therefore, the default starting bandwidth is set to \code{0.15}, the
#'suggested starting bandwidth by Feng, Gries and Fritz (2020) for the
#'data-driven  estimation of the first derivative. The recommended initial
#'bandwidth for the second derivative, however, is \code{0.2} and not
#'\code{0.15}. Thus, if the algorithm does not give suitable results
#'(especially for \code{d = 2}), the adjustment of the initial bandwidth might
#'be a good starting point. Analogously, the default starting bandwidth for the
#'trend estimation for the variance factor is \code{bStart.p = 0.15}, although
#'according to Feng, Gries and Fritz (2020), \code{bStart.p = 0.1} is suggested
#'for \code{pp = 1} and \code{bStart.p = 0.2} for \code{pp = 3}. The default is
#'therefore a compromise between the two suggested values. For more specific
#'information on the input arguments consult the section \emph{Arguments}.
#'
#'After the bandwidth estimation, the nonparametric derivative of the series
#'is calculated with respect to the obtained optimal bandwidth by means of a
#'local polynomial regression. The output object is then a list that contains,
#'among other components, the original time series, the estimates of the
#'derivative and the estimated optimal bandwidth.
#'
#'The default print method for this function delivers key numbers such as
#'the iteration steps and the generated optimal bandwidth rounded to the fourth
#'decimal. The exact numbers and results such as the estimated nonparametric
#'trend series are saved within the output object and can be addressed via the
#'\code{$} sign.
#'
#'NOTE:
#'
#'The estimates are obtained for the rescaled time points on the interval
#'\eqn{[0, 1]}. Therefore, the estimated derivatives might not reflect the
#'derivatives for the actual time points. To rescale them, we refer the
#'user to the \code{\link{rescale}} function of the \code{smoots} package.
#'
#'With package version 1.1.0, this function implements C++ code by means
#'of the \code{\link[Rcpp:Rcpp-package]{Rcpp}} and
#'\code{\link[RcppArmadillo:RcppArmadillo-package]{RcppArmadillo}} packages for
#'better performance.
#'
#'@return The function returns a list with different components:
#'
#'\describe{
#'\item{b0}{the optimal bandwidth chosen by the IPI-algorithm.}
#'\item{bStart}{the starting bandwidth for the local polynomial
#'regression based derivative estimation procedure; input argument.}
#'\item{bStart.p}{the starting bandwidth for the nonparametric trend estimation
#'that leads to the variance factor estimate; input argument.}
#'\item{bvc}{indicates whether an enlarged bandwidth was used for the variance
#'factor estimation or not; it is always set to \code{"Y"} (yes) for this
#'function.}
#'\item{cf0}{the estimated variance factor; in contrast to the definitions
#'given in the \emph{Details} section, this object actually contains an
#'estimated value of \eqn{2\pi c_f}, i.e. it corresponds to the estimated sum
#'of autocovariances.}
#'\item{InfR}{the inflation rate setting.}
#'\item{iterations}{the bandwidths of the single iterations steps}
#'\item{Mcf}{the estimation method for the variance factor estimation; it is
#'always estimated nonparametrically (\code{"NP"}) within this function.}
#'\item{mu}{the smoothness parameter of the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{niterations}{the total number of iterations until convergence.}
#'\item{orig}{the original input series; input argument.}
#'\item{p}{the order of polynomial for the local polynomial
#'regression used within derivative estimation procedure.}
#'\item{pp}{the order of polynomial for the local polynomial
#'regression used in the variance factor estimation; input argument.}
#'\item{v}{the considered order of the trend's derivative; input argument
#'\code{d}.}
#'\item{ws}{the weighting system matrix used within the local polynomial
#'regression; this matrix is a condensed version of a complete weighting system
#'matrix; in each row of \code{ws}, the weights for conducting the smoothing
#'procedure at a specific observation time point can be found; the first
#'\eqn{[nb + 0.5]} rows, where \eqn{n} corresponds to the number of
#'observations, \eqn{b} is the bandwidth considered for smoothing and
#'\eqn{[.]} denotes the integer part, contain the weights at the
#'\eqn{[nb + 0.5]} left-hand boundary points; the weights in row
#'\eqn{[nb + 0.5] + 1} are representative for the estimation at all
#'interior points and the remaining rows contain the weights for the right-hand
#'boundary points; each row has exactly \eqn{2[nb + 0.5] + 1} elements,
#'more specifically the weights for observations of the nearest
#'\eqn{2[nb + 0.5] + 1} time points; moreover, the weights are normalized,
#'i.e. the weights are obtained under consideration of the time points
#'\eqn{x_t = t/n}, where \eqn{t = 1, 2, ..., n}.}
#'\item{ye}{the nonparametric estimates of the derivative for the rescaled
#'time points on the interval \eqn{[0, 1]}.}
#'}
#'
#'@md
#'
#'@export
#'
#'@references
#' Feng, Y., Gries, T. and Fritz, M. (2020). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Journal of Nonparametric Statistics, 32:2, 510-533.
#'
#' Feng, Y., Gries, T., Letmathe, S. and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University. Unpublished.
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#'@examples
#'# Logarithm of test data
#'test_data <- gdpUS
#'y <- log(test_data$GDP)
#'t <- seq(from = 1947, to = 2019.25, by = 0.25)
#'
#'# Applied dsmooth function for the trend's first derivative
#'result_d <- dsmooth(y, d = 1, mu = 1, pp = 1, bStart.p = 0.1, bStart = 0.15)
#'estim <- result_d$ye
#'
#'# Plot of the results
#'plot(t, estim, xlab = "Year", ylab = "First derivative", type = "l",
#'  main = paste0("Estimated first derivative of the trend for log-quarterly ",
#'  "US-GDP, Q1 1947 - Q2 2019"), cex.axis = 0.8, cex.main = 0.8,
#'  cex.lab = 0.8, bty = "n")
#'
#'# Print result
#'result_d


# The main function "dsmooth"------------------------------------------

# d = 1 or 2 indicate iterative plug-in for m' or m'', respectively
dsmooth <- function(y, d = c(1, 2), mu = c(0, 1, 2, 3), pp = c(1, 3),
                    bStart.p = 0.15, bStart = 0.15) {

  if (length(y) <= 1 || !all(!is.na(y)) || !is.numeric(y)) {
    stop("The argument 'y' must be a numeric vector with length > 1 and ",
         "without NAs.")
  }

  if (!(length(d) %in% c(1, 2)) || !all(!is.na(d)) || !is.numeric(d)) {
    stop("The argument 'd' must be a single integer value (1 or 2).")
  }
  d <- floor(d)

  if (!(length(mu) %in% c(1, 4)) || !all(!is.na(mu)) || !is.numeric(mu)) {
    stop("The argument 'mu' must be a single integer value (0, 1, 2 or 3).")
  }
  mu <- floor(mu)

  if (!(length(pp) %in% c(1, 2)) || !all(!is.na(pp)) || !is.numeric(pp)) {
    stop("The argument 'pp' must be a single integer value (either 1 or 3).")
  }
  pp <- floor(pp)

  if (length(bStart.p) != 1 || is.na(bStart.p) || !is.numeric(bStart.p) ||
      (bStart <= 0)) {
    stop("The argument 'bStart.p' must be a single non-NA double value with ",
         "bStart.p > 0.")
  }

  if (length(bStart) != 1 || is.na(bStart) || !is.numeric(bStart) ||
      (bStart <= 0)) {
    stop("The argument 'bStart' must be a single non-NA double value with ",
         "bStart > 0.")
  }

  if (all(d == c(1, 2))) d <- 1
  if (all(mu == c(0, 1, 2, 3))) mu <- 1
  if (all(pp == c(1, 3))) pp <- 1

  if (!(d %in% c(1, 2))) {
    stop("Input of argument 'd' incorrect. It must be set to either 1 or 2.")
  }
  if (!(pp %in% c(1, 3))) {
    stop("Input of argument 'pp' incorrect. It must be set to either 1 or 3.")
  }


  if (pp == 1) {result.p <- msmoothCalc(y, 1, mu, bStart.p, "A", "lpr")}
  if (pp == 3) {result.p <- msmoothCalc(y, 3, mu, bStart.p, "B", "lpr")}
  cf0 <- result.p$cf0

  # Using the knn idea, bb=1, or not
  bb <- 1  # default
  if (d == 1) {InfR <- "Nai"} else {InfR <- "Var"}  # default
  cb <- 0.05
  # In this code ccf and bStart will be provided by the user
  # Input parameters
  n <- length(y)
  p <- d + 1
  k <- p + 1
  pd <- p + 2
  runc <- 1
  n1 <- trunc(n * cb)

  # New method for the kernel constants with p = 1 or 3------------------------

  # Kernel constants-----------------------------------------------------------
  m <- 1000000  # for the numerical integral
  u <- (-m:m) / (m + 0.5)
  # For d = 1, the four 3rd-order kernels for m' (Table 5.7, Mueller, 1988)

  wkp <- lookup$lookup_3kerns[[(mu + 1), d]](u)

  Rp <- sum(wkp ** 2) / m
  mukp <- sum((u ** k) * wkp) / m

  # Two constant in the bandwidth
  c1 <- factorial(k) ** 2 * (2 * d + 1) / (2 * (k - d))  # for d = 1 or 2
  c2 <- (1 - 2 * cb) * Rp / (mukp) ** 2

  steps <- rep(NA, 40)
  bd2_func <- lookup$InfR2_lookup[[d, InfR]]
  expoDer1 <- lookup$exponDer1[[d]]
  expoDer2 <- lookup$exponDer2[[d]]

  # The main iteration---------------------------------------------------------

  noi <- 40  # maximal number of iterations

  for (i in 1:noi) {
    if (runc == 1) {
      if (i > 1) {bold1 <- bold}
      if (i == 1) {bold <- bStart} else {bold <- bopt}

      # Look up EIM inflation rates from the lookup table
      bd <- bd2_func(bold)

      if (bd >= 0.49) {bd <- 0.49}
      yed <- c(gsmoothCalcCpp(y, k, pd, mu, bd, bb))
      I2 <- sum(yed[(n1 + 1):(n - n1)] ** 2) / (n - 2 * n1)
	    c3 <- cf0 / I2

	    bopt <- (c1 * c2 * c3)^expoDer1 * n^(-expoDer1)
	    if (bopt < n^expoDer2) bopt <- n^expoDer2

      if (bopt > 0.49) {bopt <- 0.49}
      steps[i] = bopt

      if (i > 2 && abs(bold - bopt) / bopt < 1 / n) {runc <- 0}
      if (i > 3 && abs(bold1 - bopt) / bopt < 1 / n){
        bopt <- (bold + bopt) / 2
        runc <- 0
      }
    }
  }

  # Smooth with the selected bandwidth-----------------------------------------

  if (bopt < n^expoDer2) bopt <- n^expoDer2
  if (bopt > 0.49) {bopt <- 0.49}
  est.opt <- gsmoothCalc2Cpp(y, d, p, mu, bopt, bb)
  ye <- c(est.opt$ye)
  ws <- est.opt$ws
  # Final results
  results <- list(ye = ye, orig = y, b0 = bopt, ws = ws,
                  cf0 = cf0, v = d, n = length(y),
                  iterations = steps[!is.na(steps)],
                  niterations = length(steps[!is.na(steps)]),
                  pp = pp, bStart.p = bStart.p, bStart = bStart,
                  InfR = InfR, mu = mu, p = d + 1, bvc = "Y", Mcf = "NP")

  class(results) <- "smoots"
  attr(results, "function") <- "dsmooth"

  results
}
# End of the function
