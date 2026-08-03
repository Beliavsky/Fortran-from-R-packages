#' Nonparametric calculation of univariate Value at Risk and Expected Shortfall
#'
#' Computes Value at Risk and Expected Shortfall (Conditional Value
#' at Risk) by means of plain and age-weighted historical simulation.
#'
#' @param x a numeric vector of asset returns
#' @param p confidence level for VaR calculation; default is \code{0.975}
#' @param method method to be used for calculation; default is \code{'plain'}
#' @param lambda decay factor for the calculation of weights; default is \code{0.98}
#'
#' @export
#'
#' @return Returns a list with the following elements:
#' \describe{
#' \item{VaR}{Calculated Value at Risk}
#' \item{ES}{Calculated Expected Shortfall (Conditional Value at Risk)}
#' \item{p}{Confidence level for VaR calculation}
#' }
#' @examples
#' prices <- DAX$price_close
#' returns <- diff(log(prices))
#' hs(x = returns, p = 0.975, method = 'plain')
#' hs(x = returns, p = 0.975, method = 'age', lambda = 0.98)

hs <- function(x, p = 0.975, method = c("age", "plain"), lambda = 0.98) {
    if (length(x) <= 1 || any(is.na(x)) || !is.numeric(x)) {
        stop("A numeric vector of length > 1 and without NAs must be passed to",
             " 'x'.")
    }
    if (!length(method) %in% c(1, 2) || any(is.na(method)) ||
        !is.character(method) || !all(method %in% c("age", "plain"))) {
        stop("A single character value must be passed to 'method'.",
             "Valid choices are 'age' or 'plain'.")
    }
    if (length(p) != 1 || is.na(p) || !is.numeric(p) || (p <= 0)) {
        stop("The argument 'p' must be a single non-NA double value with ",
             "0 < p < l.")
    }
    if (length(lambda) != 1 || is.na(lambda) || !is.numeric(lambda) ||
        lambda < 0 || lambda >= 1) {
        stop("The argument 'lambda' must be a single non-NA double value with ",
             "0 <= lambda < 1.")
    }

    if (all(method == c("age", "plain")))
        method <- "plain"

    if (method == "plain") {
        l <- sort(-x)
        VaR <- unname(stats::quantile(l, p))
        ES <- mean(l[l > VaR])
    }

    if (method == "age") {
        n <- length(x)
        w <- rev(lambda^((1:n) - 1) * (1 - lambda) / (1 - lambda^n))
        l <- sort(-x)
        l.ind <- sort(-x, index.return = TRUE)$ix
        pcum <- cumsum(w[l.ind])
        ind.high <- which(pcum > p)[1]
        ind.low <- ind.high - 1
        VaR.high <- l[ind.high]
        VaR.low <- l[ind.low]
        VaR <- VaR.low + (p - pcum[ind.low]) * (VaR.high - VaR.low) /
            (pcum[ind.high] - pcum[ind.low])
        ind.ES <- l.ind[which(l > VaR)]
        w.ES <- w[ind.ES]
        w.ES <- w.ES / sum(w.ES)
        l.ES <- l[l > VaR]
        ES <- sum(l.ES * w.ES)
    }
    results <- list(VaR = VaR, ES = ES, p = p)
    results
}
