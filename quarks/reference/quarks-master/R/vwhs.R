#' Volatility weighted historical simulation
#'
#' Calculates univariate Value at Risk and Expected Shortfall (Conditional
#' Value at Risk) by means of volatility weighted historical simulation.
#' Volatility can be estimated with an exponentially weighted moving average
#' or a GARCH-type model.
#'
#' @param x a numeric vector of asset returns
#' @param p confidence level for VaR calculation; default is \code{0.975}
#' @param model model for estimating conditional volatility; default is \code{'EWMA'}
#' @param lambda decay factor for the calculation of weights; default is \code{0.94}
#' @param ... additional arguments of the \code{ugarchspec} function from the
#' \code{rugarch}-package; the default settings for the arguments
#' \code{variance.model} and \code{mean.model} are \code{list(model = 'sGARCH',
#' garchOrder = c(1, 1))} and \code{list(armaOrder = c(0, 0))}, respectively
#'
#' @export
#'
#' @return Returns a list with the following elements:
#' \describe{
#' \item{VaR}{Calculated Value at Risk}
#' \item{ES}{Calculated Expected Shortfall (Conditional Value at Risk)}
#' \item{p}{Confidence level for VaR calculation}
#' \item{garchmod}{The model fit. Is the respective GARCH fit for
#' \code{model = 'GARCH'} (see \code{rugarch} documentation) and  \code{'EWMA'} for
#' \code{model = 'EWMA'}}
#' }
#' @examples
#' prices <- DAX$price_close
#' returns <- diff(log(prices))
#' # volatility weighting via EWMA
#' ewma <- vwhs(x = returns, p = 0.975, model = "EWMA", lambda = 0.94)
#' ewma
#' # volatility weighting via GARCH
#' garch <- vwhs(x = returns, p = 0.975, model = "GARCH", variance.model =
#' list(model = "sGARCH"))
#' garch

vwhs <- function(x, p = 0.975, model = c("EWMA", "GARCH"), lambda = 0.94, ...)
    {
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
    if (all(model == c("EWMA", "GARCH")))
        model <- "EWMA"
    dots <- list(...)
    if (!"mean.model" %in% names(dots)) {
        dots[["mean.model"]] <- list(armaOrder = c(0, 0))
    }
    if (model == "GARCH") {
        spec <- do.call(what = rugarch::ugarchspec, args = dots)
        fit <- rugarch::ugarchfit(data = x, spec = spec)
        csig <- as.numeric(rugarch::sigma(fit))
    }
    if (model == "EWMA") {
        cvar <- ewma(x, lambda = lambda)
        csig <- sqrt(cvar)
        fit <- "EWMA"
    }
    n <- length(x)
    xz <- x / csig
    loss <- -(xz * csig[n])
    VaR <- unname(stats::quantile(loss, p))
    ES <- mean(loss[loss > VaR])
    results <- list(VaR = VaR, ES = ES, p = p, garchmod = fit)
    results
}
