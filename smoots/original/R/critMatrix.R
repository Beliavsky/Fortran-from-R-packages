#'ARMA Order Selection Matrix
#'
#'An information criterion is calculated for different orders of an
#'autoregressive-moving-average (ARMA) model.
#'
#'@param X a numeric vector that contains the observed time series ordered
#'from past to present; the series is assumed to follow an ARMA process.
#'@param p.max an integer value \eqn{>= 0} that defines the maximum
#'autoregressive order to calculate the criterion for; is set to \code{5}
#'by default; decimal numbers will be rounded off to integers.
#'@param q.max an integer value \eqn{>= 0} that defines the maximum
#'moving-average order to to calculate the criterion for; is set to \code{5}
#'by default; decimal numbers will be rounded off to integers.
#'@param criterion a character value that defines the information criterion
#'that will be calculated; the Bayesian Information Criterion (\code{"bic"})
#'and Akaike Information Criterion (\code{"aic"}) are the supported choices;
#'is set to \code{"bic"} by default.
#'@param include.mean a logical value; this argument regulates whether to
#'estimate the mean of the series (\code{TRUE}) or not (\code{FALSE}); is set
#'to \code{TRUE} by default.
#'
#'@export
#'
#'@importFrom stats arima
#'
#'@details
#'This function is part of the \code{smoots} package and was implemented under
#'version 1.1.0. The series passed to \code{X} is assumed to follow an
#'ARMA(\eqn{p,q}) model. A \code{p.max + 1} by \code{q.max + 1} matrix is
#'calculated for this series. More precisely, the criterion chosen via the
#'argument \code{criterion} is calculated for all combinations of orders
#'\eqn{p = 0, 1, ..., p_{max}}{p = 0, 1, ..., p_max} and
#'\eqn{q = 0, 1, ..., q_{max}}{q = 0, 1, ..., q_max}.
#'
#'Within the function, two information criteria are supported: the Bayesian
#'Information Criterion (BIC) and Akaike's Information Criterion (AIC). The AIC
#'is given by
#'\deqn{AIC_{p,q} := \ln(\hat{\sigma}_{p,q}^{2}) + \frac{2(p+q)}{n},}{AIC_[p,q]
#':= ln(hat[sigma]^[2]_[p,q]) + \{2(p + q)\}/n,}
#'where \eqn{\hat{sigma}_{p,q}^{2}}{hat[sigma]^[2]_[p,q]} is the estimated
#'innovation variance, \eqn{p} and \eqn{q} are the ARMA orders and \eqn{n} is
#'the number of observations.
#'
#'The BIC, on the other hand, is defined by
#'\deqn{BIC_{p,q} := k \ln(n) - 2\ln(\hat{L})}{BIC_[p,q] := k * ln(n) -
#'2ln(hat[L])}
#'with \eqn{k} being the number of estimated parameters and
#'\eqn{\hat{L}}{hat[L]} being the estimated Log-Likelihood. Since the parameter
#'\eqn{k} only differs with respect to the orders \eqn{p} and \eqn{q} for all
#'estimated models, the term \eqn{k \ln(n)}{k * ln(n)} is reduced to
#'\eqn{(p + q) \ln(n)}{(p + q) * ln(n)} within the function. Exemplarily,
#'if the mean of the series is estimated as well, it is usually considered
#'within the parameter \eqn{k} when calculating the BIC.
#'However, since the mean is estimated for all models, not considering this
#'estimated parameter within the calculation of the BIC will reduce all BIC
#'values by the same amount of \eqn{\ln(n)}{ln(n)}. Therefore, the selection
#'via this simplified criterion is still valid, if the number of the estimated
#'parameters only differs with respect to \eqn{p} and \eqn{q} between the
#'models that the BIC is obtained for.
#'
#'The optimal orders are considered to be the ones which minimize either the
#'BIC or the AIC. The use of the BIC is however recommended, because the BIC
#'is consistent, whereas the AIC is not.
#'
#'NOTE:
#'
#'Within this function, the \code{\link[stats]{arima}} function of the
#'\code{stats} package with its method \code{"CSS-ML"} is used throughout for
#'the estimation of ARMA models.
#'
#'@return
#'The function returns a \code{p.max + 1} by \code{q.max + 1} matrix, where the
#'rows represent the AR orders from \eqn{p = 0} to \eqn{p = p_{max}}{p = p_max}
#'and the columns represent the MA orders from \eqn{q = 0} to
#'\eqn{q = q_{max}}{q = q_max}. The values within the matrix are the values of
#'the previously selected information criterion for the different combinations
#'of \eqn{p} and \eqn{q}.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#'@examples
#'\dontrun{
#'# Simulate an ARMA(2,1) process
#'set.seed(23)
#'X.sim <- stats::arima.sim(model = list(ar = c(1.2, -0.71), ma = 0.46),
#'  n = 1000) + 13.1
#'# Application of the function
#'critMatrix(X.sim)
#'# Result: Via the BIC, the orders p.opt = 2 and q.opt = 1 are selected.
#'}

critMatrix <- function(X, p.max = 5, q.max = 5, criterion = c("bic", "aic"),
  include.mean = TRUE) {

  if (length(X) <= 1 || !all(!is.na(X)) || !is.numeric(X)) {
    stop("The argument 'X' must be a numberic vector with length ",
      "significantly larger than 1.")
  }
  if (length(p.max) != 1 || is.na(p.max) || !is.numeric(p.max) || p.max < 0) {
    stop("The argument 'p.max' must be a single integer >= 0.")
  }
  p.max <- floor(p.max)
  if (length(q.max) != 1 || is.na(q.max) || !is.numeric(q.max) || q.max < 0) {
    stop("The argument 'q.max' must be a single integer >= 0.")
  }
  q.max <- floor(q.max)

  if (!(length(criterion) %in% c(1, 2)) || !all(!is.na(criterion)) ||
      !is.character(criterion)) {
    stop("Input to argument 'criterion' not recognized.")
  }
  if (all(criterion == c("bic", "aic"))) criterion <- "bic"
  if (length(criterion) != 1 || !(criterion %in% c("bic", "aic"))) {
    stop("Input to argument 'criterion' not recognized.")
  }
  if (length(include.mean) != 1 || is.na(include.mean) ||
      !(include.mean %in% c(TRUE, FALSE))) {
    stop("The argument 'include.mean' must be a single logical value ",
         "(TRUE or FALSE).")
  }

  n <- length(X)
  matOut <- matrix(NA, ncol = q.max + 1, nrow = p.max + 1)
  for (i in 0:p.max) {
    for (j in 0:q.max) {
      suppressWarnings(
        est <- arima(X, order = c(i, 0, j), include.mean = include.mean)
      )
      if (criterion == "bic") {
        matOut[(i + 1), (j + 1)] <- -2 * est$loglik + log(n) * (i + j)
      } else {
        matOut[(i + 1), (j + 1)] <- est$aic
      }
    }
  }
  rownames(matOut) <- paste0("p=", 0:p.max)
  colnames(matOut) <- paste0("q=", 0:q.max)
  matOut
}
