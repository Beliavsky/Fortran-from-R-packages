#' Filtered historical simulation
#'
#' Calculates univariate Value at Risk and Expected Shortfall (Conditional
#' Value at Risk) by means of filtered historical simulation.
#' Volatility can be estimated with an exponentially weighted moving
#' average or a GARCH-type model.
#'
#' @param x a numeric vector of asset returns
#' @param p confidence level for VaR calculation; default is \code{0.975}
#' @param model model for estimating conditional volatility; options are \code{'EWMA'}
#' and \code{'GARCH'}; if \code{model = 'GARCH'}, additional arguments can be adjusted
#' via \code{...}; default is \code{'EWMA'}
#' @param lambda decay factor for the calculation of weights; default is \code{0.94}
#' @param nboot size of bootstrap sample; must be a single non-NA integer value
#' with \code{nboot > 0}; default is \code{NULL}
#' @param ... additional arguments of the \code{ugarchspec} function from the
#' \code{rugarch}-package; only applied if \code{model = 'GARCH'}; default
#' settings for the arguments \code{variance.model} and \code{mean.model} are:
#' \describe{
#' \item{\code{variance.model} = \code{list(model = 'sGARCH', garchOrder = c(1, 1))}}{}
#' \item{\code{mean.model} = \code{list(armaOrder = c(0, 0))}}{}
#' }
#'
#' @export
#'
#' @return Returns a list with the following elements:
#' \describe{
#' \item{VaR}{Calculated Value at Risk}
#' \item{ES}{Calculated Expected Shortfall (Conditional Value at Risk)}
#' \item{p}{Confidence level for VaR calculation}
#' \item{garchmod}{The model fit. Is the respective GARCH fit for
#' \code{model = "GARCH"} (see \code{rugarch} documentation) and  \code{'EWMA'} for
#' \code{model = "EWMA"}}
#' }
#' @examples
#' prices <- DAX$price_close
#' returns <- diff(log(prices))
#' # volatility weighting via EWMA
#' ewma <- fhs(x = returns, p = 0.975, model = "EWMA", lambda = 0.94,
#'             nboot = 10000)
#' ewma
#' # volatility weighting via GARCH
#' garch <- fhs(x = returns, p = 0.975, model = "GARCH", variance.model =
#' list(model = "sGARCH"), nboot = 10000)
#' garch

fhs <- function(x, p = 0.975, model = c("EWMA", "GARCH"), lambda = 0.94,
                nboot = NULL, ...) {
  if (length(x) <= 1 || any(is.na(x)) || !is.numeric(x)) {
    stop("A numeric vector of length > 1 and without NAs must be passed to",
         " 'x'.")
  }
  if (length(p) != 1 || is.na(p) || !is.numeric(p) || (p <= 0)) {
    stop("The argument 'p' must be a single non-NA double value with ",
         "0 < p < l.")
  }
  if (!length(model) %in% c(1, 2) || any(is.na(model)) ||
      !is.character(model) || !all(model %in% c("EWMA", "GARCH"))) {
    stop("A single character value must be passed to 'model'.",
         "Valid choices are 'EWMA' or 'GARCH'.")
  }
  if (length(lambda) != 1 || is.na(lambda) || !is.numeric(lambda) ||
      lambda < 0 || lambda >= 1) {
    stop("The argument 'lambda' must be a single non-NA double value with ",
         "0 < lambda < 1.")
  }
  if (length(nboot) != 1 || is.na(nboot) || !is.numeric(nboot) ||
      nboot <= 0) {
    stop("The argument 'nboot' must be a single non-NA integer value with ",
         "nboot > 0.")
  }
  if (all(model == c("EWMA", "GARCH")))
    model <- "EWMA"
  dots = list(...)
  if (!"mean.model" %in% names(dots)) {
    dots[["mean.model"]] <- list(armaOrder = c(0, 0))
  }
  n <- length(x)
  if (model == "GARCH") {
    spec <- do.call(what = rugarch::ugarchspec, args = dots)
    fit <- rugarch::ugarchfit(data = x, spec = spec)
    csig <- as.numeric(rugarch::sigma(fit))
    one.ahead.csig <-
      as.numeric(rugarch::sigma(rugarch::ugarchforecast(fit, n.ahead = 1)))
  }
  if (model == "EWMA") {
    cvar <- ewma(x, lambda = lambda)
    one.ahead.cvar <-  lambda * cvar[n] + (1 - lambda) * x[n]^2
    csig <- sqrt(cvar)
    one.ahead.csig <- sqrt(one.ahead.cvar)
    fit = "EWMA"
  }
  xz <- x / csig
  boot.xz <- sample(xz, size = nboot, replace = TRUE)
  boot.loss <- -(boot.xz * one.ahead.csig)
  VaR <- unname(stats::quantile(boot.loss, p))
  ES <- mean(boot.loss[boot.loss > VaR])
  results <- list(VaR = VaR, ES = ES, p = p, garchmod = fit)
  results
}
