#' Data-driven Nonparametric Regression for the Trend in Equidistant Time
#' Series
#'
#' This function runs an iterative plug-in algorithm to find the optimal
#' bandwidth for the estimation of the nonparametric trend in equidistant
#' time series (with short memory errors) and then employs the resulting
#' bandwidth via either local polynomial or kernel regression.
#'
#'@param y a numeric vector that contains the input time series ordered from
#'past to present.
#'@param p an integer \code{1} (local linear regression) or \code{3} (local
#'cubic regression); represents the order of polynomial within the local
#'polynomial regression (see also the 'Details' section); is set to \code{1}
#'by default; is automatically set to \code{1} if \code{method = "kr"}.
#'@param mu an integer \code{0}, ..., \code{3} that represents the smoothness
#'parameter of the kernel weighting function and thus defines the kernel
#'function that will be used within the local polynomial regression; is set
#'to \code{1} by default.
#'
#'\tabular{cl}{
#'\strong{Number} \tab \strong{Kernel}\cr
#'\code{0} \tab Uniform Kernel\cr
#'\code{1} \tab Epanechnikov Kernel\cr
#'\code{2} \tab Bisquare Kernel\cr
#'\code{3} \tab Triweight Kernel
#'}
#'@param bStart a numeric object that indicates the starting value of the
#'bandwidth for the iterative process; should be \eqn{> 0}; is set to
#'\code{0.15} by default.
#'@param alg a control parameter (as character) that indicates the
#'corresponding algorithm used (set to \code{"A"} by default for \code{p = 1}
#'and to \code{"B"} for \code{p = 3}).
#'
#'\tabular{cl}{
#'\strong{Algorithm} \tab \strong{Description}\cr
#'\code{"A"} \tab Nonparametric estimation of the variance factor with an
#'enlarged bandwidth, optimal inflation rate\cr
#'\code{"B"} \tab Nonparametric estimation of the variance factor with an
#'enlarged bandwidth, naive inflation rate \cr
#'\code{"O"} \tab Nonparametric estimation of the variance factor, optimal
#'inflation rate\cr
#'\code{"N"} \tab Nonparametric estimation of the variance factor, naive
#'inflation rate\cr
#'\code{"OAM"} \tab Estimation of the variance factor with
#'ARMA(\eqn{p,q})-models, optimal inflation rate\cr
#'\code{"NAM"} \tab Estimation of the variance factor with
#'ARMA(\eqn{p,q})-models, naive inflation rate\cr
#'\code{"OA"} \tab Estimation of the variance factor with
#'AR(\eqn{p})-models, optimal inflation rate\cr
#'\code{"NA"} \tab Estimation of the variance factor with
#'AR(\eqn{p})-models, naive inflation rate\cr
#'\code{"OM"} \tab Estimation of the variance factor with
#'MA(\eqn{q})-models, optimal inflation rate\cr
#'\code{"NM"} \tab Estimation of the variance factor with
#'MA(\eqn{q})-models, naive inflation rate
#'}
#'
#'It is proposed to use \code{alg = "A"} in combination with \code{p = 1}.
#'If the user finds that the chosen bandwidth by algorithm \code{"A"} is too
#'small, \code{alg = "B"} with preferably \code{p = 3} is suggested. For more
#'information on the components of the different algorithms, please consult
#'\code{\link{tsmooth}}.
#'@param method the smoothing approach; \code{"lpr"} represents the local
#'polynomial regression, whereas \code{"kr"} implements a kernel regression;
#'is set to \code{"lpr"} by default.
#'
#'@export
#'
#'@details
#'
#'The trend is estimated based on the additive
#'nonparametric regression model for an equidistant time series
#'\deqn{y_t = m(x_t) + \epsilon_t,}
#'where \eqn{y_t} is the observed time series, \eqn{x_t} is the rescaled time
#'on the interval \eqn{[0, 1]}, \eqn{m(x_t)} is a smooth and deterministic
#'trend function and \eqn{\epsilon_t} are stationary errors with
#'\eqn{E(\epsilon_t) = 0} and short-range dependence (see also Beran and Feng,
#'2002). With this function \eqn{m(x_t)} can be estimated without a parametric
#'model assumption for the error series. Thus, after estimating and removing
#'the trend, any suitable parametric model, e.g. an ARMA(\eqn{p,q}) model, can
#'be fitted to the residuals (see \code{\link[stats]{arima}}).
#'
#'The iterative-plug-in (IPI) algorithm, which numerically minimizes the
#'Asymptotic Mean Squared Error (AMISE), was proposed by Feng, Gries
#'and Fritz (2020).
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
#'The function calculates suitable estimates for \eqn{c_f}, the variance
#'factor, and \eqn{I[m^{(k)}]}{I[m^(k)]} over different iterations. In each
#'iteration, a bandwidth is obtained in accordance with the AMISE that once
#'more serves as an input for the following iteration. The process repeats
#'until either convergence or the 40th iteration is reached. For further
#'details on the asymptotic theory or the algorithm, please consult Feng,
#'Gries and Fritz (2020) or Feng et al. (2019).
#'
#'To apply the function, only few arguments are needed: a data input \code{y},
#'an order of polynomial \code{p}, a kernel function defined by the smoothness
#'parameter \code{mu}, a starting value for the relative bandwidth
#'\code{bStart} and a final smoothing method \code{method}.
#'In fact, aside from the input vector \code{y}, every argument has a default
#'setting that can be adjusted for the individual case. It is recommended to
#'initially use the default values for \code{p}, \code{alg} and
#'\code{bStart} and adjust them in the rare case of the resulting optimal
#'bandwidth being either too small or too large. Theoretically, the
#'initial bandwidth does not affect the selected optimal bandwidth. However, in
#'practice local minima of the AMISE might exist and influence the selected
#'bandwidth. Therefore, the default setting is \code{bStart = 0.15}, which is a
#'compromise between the starting values \code{bStart = 0.1} for \code{p = 1}
#'and \code{bStart = 0.2} for \code{p = 3} that were proposed by Feng, Gries
#'and Fritz (2020). In the rare case of a clearly unsuitable optimal bandwidth,
#'a starting bandwidth that differs from the default value is a first
#'possible approach to obtain a better result. Other argument adjustments can
#'be tried as well. For more specific information on the input arguments
#'consult the section \emph{Arguments}.
#'
#'When applying the function, an optimal bandwidth is obtained based on the
#'IPI algorithm proposed by Feng, Gries and Fritz (2020). In a second step,
#'the nonparametric trend of the series is calculated with respect
#'to the chosen bandwidth and the selected regression method (\code{lpf} or
#'\code{kr}). It is notable that \code{p} is automatically set to 1 for
#'\code{method = "kr"}. The output object is then a list that contains, among
#'other components, the original time series, the estimated trend values and
#'the series without the trend.
#'
#'The default print method for this function delivers key numbers such as
#'the iteration steps and the generated optimal bandwidth rounded to the fourth
#'decimal. The exact numbers and results such as the estimated nonparametric
#'trend series are saved within the output object and can be addressed via the
#'\code{$} sign.
#'
#'NOTE:
#'
#'With package version 1.1.0, this function implements C++ code by means
#'of the \code{\link[Rcpp:Rcpp-package]{Rcpp}} and
#'\code{\link[RcppArmadillo:RcppArmadillo-package]{RcppArmadillo}} packages for
#'better performance.
#'
#'@md
#'
#'@return The function returns a list with different components:
#'
#'\describe{
#'\item{AR.BIC}{the Bayesian Information Criterion of the optimal AR(\eqn{p})
#'model when estimating the variance factor via autoregressive models
#'(if calculated; calculated for \code{alg = "OA"} and \code{alg = "NA"}).}
#'\item{ARMA.BIC}{the Bayesian Information Criterion of the optimal
#'ARMA(\eqn{p,q}) model when estimating the variance factor via
#'autoregressive-moving-average models (if calculated; calculated for
#'\code{alg = "OAM"} and \code{alg = "NAM"}).}
#'\item{cb}{the percentage of omitted observations on each side of the
#'observation period; always equal to 0.05.}
#'\item{b0}{the optimal bandwidth chosen by the IPI-algorithm.}
#'\item{bb}{the boundary bandwidth method used within the IPI; always equal to
#'1.}
#'\item{bStart}{the starting value of the (relative) bandwidth; input
#'argument.}
#'\item{bvc}{indicates whether an enlarged bandwidth was used for the variance
#'factor estimation or not; depends on the chosen algorithm.}
#'\item{cf0}{the estimated variance factor; in contrast to the definitions
#'given in the \emph{Details} section, this object actually contains an
#'estimated value of \eqn{2\pi c_f}, i.e. it corresponds to the estimated sum
#'of autocovariances.}
#'\item{cf0.AR}{the estimated variance factor obtained by estimation of
#'autoregressive models (if calculated; \code{alg = "OA"} or \code{"NA"}).}
#'\item{cf0.ARMA}{the estimated variance factor obtained by estimation of
#'autoregressive-moving-average models (if calculated; calculated for
#'\code{alg = "OAM"} and \code{alg = "NAM"}).}
#'\item{cf0.LW}{the estimated variance factor obtained by Lag-Window Spectral
#'Density Estimation following Bühlmann (1996) (if calculated; calculated for
#'algorithms \code{"A"}, \code{"B"}, \code{"O"} and \code{"N"}).}
#'\item{cf0.MA}{the estimated variance factor obtained by estimation of
#'moving-average models (if calculated; calculated for \code{alg = "OM"} and
#'\code{alg = "NM"}).}
#'\item{I2}{the estimated value of \eqn{I[m^{(k)}]}{I[m^(k)]}.}
#'\item{InfR}{the setting for the inflation rate according to the chosen
#'algorithm.}
#'\item{iterations}{the bandwidths of the single iterations steps}
#'\item{L0.opt}{the optimal bandwidth for the lag-window spectral density
#'estimation (if calculated; calculated for algorithms \code{"A"}, \code{"B"},
#'\code{"O"} and \code{"N"}).}
#'\item{MA.BIC}{the Bayesian Information Criterion of the optimal MA(\eqn{q})
#'model when estimating the variance factor via moving-average models (if
#'calculated; calculated for \code{alg = "OM"} and \code{alg = "NM"}).}
#'\item{Mcf}{the estimation method for the variance factor estimation; depends
#'on the chosen algorithm.}
#'\item{mu}{the smoothness parameter of the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{niterations}{the total number of iterations until convergence.}
#'\item{orig}{the original input series; input argument.}
#'\item{p.BIC}{the order p of the optimal AR(\eqn{p}) or ARMA(\eqn{p,q}) model
#'when estimating the variance factor via autoregressive or
#'autoregressive-moving average models (if calculated; calculated for
#'\code{alg = "OA"}, \code{alg = "NA"}, \code{alg = "OAM"} and
#'\code{alg = "NAM"}).}
#'\item{p}{the order of polynomial used in the IPI-algorithm; also used for the
#'final smoothing, if \code{method = "lpr"}; input argument.}
#'\item{q.BIC}{the order \eqn{q} of the optimal MA(\eqn{q}) or ARMA(\eqn{p,q})
#'model when estimating the variance factor via moving-average or
#'autoregressive-moving average models (if calculated; calculated for
#'\code{alg = "OM"},
#'\code{alg = "NM"}, \code{alg = "OAM"} and \code{alg = "NAM"}).}
#'\item{res}{the estimated residual series.}
#'\item{v}{the considered order of derivative of the trend; is always zero for
#'this function.}
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
#'\item{ye}{the nonparametric estimates of the trend.}
#'}
#'
#'@references
#' Beran, J. and Feng, Y. (2002). Local polynomial fitting with long-memory,
#' short-memory and antipersistent errors. Annals of the Institute of
#' Statistical Mathematics, 54(2), 291-311.
#'
#' Bühlmann, P. (1996). Locally adaptive lag-window spectral estimation.
#' Journal of Time Series Analysis, 17(3), 247-270.
#'
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
#'### Example 1: US-GDP ###
#'
#'# Logarithm of test data
#'# -> the logarithm of the data is assumed to follow the additive model
#'test_data <- gdpUS
#'y <- log(test_data$GDP)
#'
#'# Applied msmooth function for the trend
#'results <- msmooth(y, p = 1, mu = 1, bStart = 0.1, alg = "A", method = "lpr")
#'res <- results$res
#'ye <- results$ye
#'
#'# Plot of the results
#'t <- seq(from = 1947, to = 2019.25, by = 0.25)
#'matplot(t, cbind(y, ye), type = "ll", lty = c(3, 1), col = c(1, "red"),
#'  xlab = "Years", ylab = "Log-Quartlery US-GDP",
#'  main = "Log-Quarterly US-GDP vs. Trend, Q1 1947 - Q2 2019")
#'legend("bottomright", legend = c("Original series", "Estimated trend"),
#'  fill = c(1, "red"), cex = 0.7)
#'results
#'
#'\dontrun{
#'### Example 2: German Stock Index ###
#'
#'# The following procedure can be considered, if (log-)returns are assumed
#'# to follow a model from the general class of semiparametric GARCH-type
#'# models (including Semi-GARCH, Semi-Log-GARCH and Semi-APARCH models among
#'# others) with a slowly changing variance over time due to a deterministic,
#'# nonparametric scale function.
#'
#'# Obtain the logarithm of the squared returns
#'returns <- diff(log(dax$Close))   # (log-)returns
#'rt <- returns - mean(returns)     # demeaned (log-)returns
#'yt <- log(rt^2)                   # logarithm of the squared returns
#'
#'# Apply 'smoots' function to the log-data, because the logarithm of
#'# the squared returns follows an additive model with a nonparametric trend
#'# function, if the returns are assumed to follow a semiparametric GARCH-type
#'# model.
#'
#'# In this case, the setting 'alg = "A"' is used in combination with p = 3, as
#'# the resulting estimates appear to be more suitable than for 'alg = "B"'.
#'est <- msmooth(yt, p = 3, alg = "A")
#'m_xt <- est$ye                    # estimated trend values
#'
#'# Obtain the standardized returns 'eps' and the scale function 'scale.f'
#'res <- est$res                    # the detrended log-data
#'C <- -log(mean(exp(res)))         # an estimate of a constant value needed
#'                                   # for the retransformation
#'scale.f <- exp((m_xt - C) / 2)    # estimated values of the scale function in
#'                                   # the returns
#'eps <- rt / scale.f               # the estimated standardized returns
#'
#'# -> 'eps' can now be analyzed by any suitable GARCH-type model.
#'#    The total volatilities are then the product of the conditional
#'#    volatilities obtained from 'eps' and the scale function 'scale.f'.
#'}

msmooth <- function (y, p = c(1, 3), mu = c(0, 1, 2, 3), bStart = 0.15,
                     alg = c("A", "B", "N", "NA", "NAM", "NM", "O", "OA",
                       "OAM", "OM"),
                     method = c("lpr", "kr")) {


  if (length(y) <= 1 || !all(!is.na(y)) || !is.numeric(y)) {
    stop("The argument 'y' must be a numeric vector with length > 1 and ",
         "without NAs.")
  }

  if (!(length(p) %in% c(1, 2)) || !all(!is.na(p)) || !is.numeric(p)) {
    stop("The argument 'p' must be a single integer value (either 1 or 3).")
  }
  p <- floor(p)

  if (!(length(mu) %in% c(1, 4)) || !all(!is.na(mu)) || !is.numeric(mu)) {
    stop("The argument 'mu' must be a single integer value (0, 1, 2 or 3).")
  }
  mu <- floor(mu)

  if (length(bStart) != 1 || is.na(bStart) || !is.numeric(bStart) ||
      (bStart <= 0)) {
    stop("The argument 'bStart' must be a single non-NA double value with ",
         "bStart > 0.")
  }

  if (!(length(alg) %in% c(1, 10)) || !all(!is.na(alg)) ||
      !is.character(alg)) {
    stop("The argument 'alg' must be a single non-NA character value.")
  }

  if (!(length(method) %in% c(1, 2)) || !all(!is.na(method)) ||
      !is.character(method)) {
    stop("The argument 'method' must be a single non-NA character value.")
  }

  if (all(method == c("lpr", "kr"))) method <- "lpr"

  if (length(method) != 1 || !(method %in% c("lpr", "kr"))) {
    stop("Input of argument 'method' incorrect. Method not recognized.")
  }

  if (all(p == c(1, 3)) || method == "kr") p <- 1
  if (all(mu == c(0, 1, 2, 3))) mu <- 1

  if (length(p) != 1 || !(p %in% c(1, 3))) {
    stop("The argument 'p' must be a single integer value (either 1 or 3).")
  }

  if (all(alg == c("A", "B", "N", "NA", "NAM", "NM", "O", "OA",
                   "OAM", "OM"))) {
    if (p == 1) {
      alg <- "A"
    } else {
      alg <- "B"
    }
  }

  if (length(alg) != 1 || !(alg %in%
    c("A", "B", "N", "NA", "NAM", "NM", "O", "OA", "OAM", "OM"))) {
    stop("Input of argument 'alg' incorrect. Algorithm not recognized.")
  }

  if (length(mu) != 1 || !(mu %in% 0:3)) {
    stop("The argument 'mu' must be a single integer value (0, 1, 2 or 3).")
  }

  result <- msmoothCalc(y = y, p = p, mu = mu, bStart = bStart, alg = alg,
    method = method)

  class(result) <- "smoots"
  attr(result, "function") <- "msmooth"

  result
}
