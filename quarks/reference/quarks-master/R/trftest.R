#' Backtesting of Value-at-Risk via Traffic Light Test
#'
#' The Traffic Light Test, is applied to
#' previously calculated Value-at-Risk series.
#'
#' @param obj A list returned by the \code{rollcast} function, that contains
#' a Value-at-Risk series; any other list that follows the name conventions
#' of the \code{rollcast} function can be used as well.
#'
#' @importFrom stats 'pbinom'
#'
#' @export
#'
#' @details
#' This function uses an object returned by the \code{rollcast} function
#' of the \code{quarks} package as an input for the
#' function argument \code{obj}. A list with different elements, such as
#' the cumulative probabilities for the VaR series within \code{obj},
#' is returned. Instead of the list, only the traffic light backtesting results
#' are printed to the R console.
#'
#'
#' @return A list of class \code{quarks} is returned with the following elements.
#' \describe{
#' \item{model}{selected model for estimation}
#' \item{method}{selected method for estimation}
#' \item{p_VaR}{cumulative probability of observing the number of
#' breaches or fewer for (1 - \code{p})100\%-VaR}
#' \item{pot_VaR}{number of exceedances for (1 - \code{p})100\%-VaR}
#' \item{p}{coverage level for (1-\code{p})100\% VaR}
#' }
#'
#' @examples
#'
#' prices <- DAX$price_close
#' returns <- diff(log(prices))
#' n <- length(returns)
#' nout <- 250 # number of obs. for out-of-sample forecasting
#' nwin <- 500 # window size for rolling forecasts
#' results <- rollcast(x = returns, p = 0.975, method = 'age', nout = nout,
#'                      nwin = nwin)
#' trftest(results)
#'
#'
trftest <- function(obj) {
  if (!is.list(obj) && !is.data.frame(obj)) {
    stop("A list or data frame containing two vectors with equal
         length and without NAs must be passed to", " 'obj'.")
  }

  if (!inherits(obj, "quarks") && (length(obj[["loss"]]) <= 1 ||
                                   any(is.na(obj[["loss"]])) ||
                                   !is.numeric(obj[["loss"]]))) {
    stop("A numeric vector of length > 1 and without NAs must be passed to",
         " 'obj$loss'.")
  }

  if (!inherits(obj, "quarks") && (length(obj[["VaR"]]) <= 1 ||
                                   any(is.na(obj[["VaR"]])) ||
                                   !is.numeric(obj[["VaR"]]))) {
    stop("A numeric vector of length > 1 and without NAs must be passed to",
         " 'obj$VaR'.")
  }

  if (length(obj[["p"]]) != 1 || is.na(obj[["p"]]) || !is.numeric(obj[["p"]])) {
    stop("A single numeric value must be passed to"," 'p'")
  }

  if(inherits(obj, "quarks")) {
    loss <- -obj[["xout"]]
    model <- obj[["model"]]
    method <- obj[["method"]]
  }
  else {
    loss <- obj[["loss"]]
    if (is.null(model)) model <- NA
    if (is.null(method)) method <- NA
  }
  VaR <- obj[["VaR"]]
  p <- 1 - obj[["p"]]
  nout <- length(loss)

  pot_VaR <- length(VaR[VaR < loss])
  p_VaR <- pbinom(pot_VaR, nout, p)

  results <- list(model = model,
                  method = method,
                  p_VaR = p_VaR,
                  pot_VaR = pot_VaR,
                  p = p)
  class(results) <- "quarks"
  attr(results, "function") <- "trftest"
  results
}



