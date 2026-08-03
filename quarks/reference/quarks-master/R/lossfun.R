#' Loss Functions
#'
#' This functions allows for the calculation of loss functions in order to assess
#' the performance of models in regard to forecasting ES.
#'
#' @param obj a list that contains the following elements:
#' \describe{
#' \item{\code{loss}}{a numeric vector that contains the values of a loss series
#' ordered from past to present; is set to \code{NULL} by default}
#' \item{\code{ES}}{a numeric vector that contains the estimated values of the
#' ES for the same time points of the loss series \code{loss}; is set to
#' \code{NULL} by default}
#' }
#' Please note that a list returned by the \code{rollcast} function can be directly
#' passed to \code{lossfun}.
#' @param beta a single numeric value; a measure for the opportunity cost of
#' capital; default is \code{1e-04}.
#'
#' @export
#'
#' @details
#' Given a negative return series \code{obj$loss}, the corresponding Expected
#' Shortfall (ES) estimates \code{obj$ES} and a parameter \code{beta} that
#' defines the opportunity cost of capital, four different definitions of loss
#' functions are considered.
#'
#' @return an S3 class object, which is a list of
#' \describe{
#' \item{loss.fun1}{regulatory loss function}
#' \item{loss.fun2}{firm's loss function following Sarma et al. (2003)}
#' \item{loss.fun3}{loss function following Abad et al. (2015)}
#' \item{loss.fun4}{Feng's loss function; a compromise of regulatory and
#' firm's loss function}
#' }
#'
#' @references
#' Abad, P., Muela, S. B., & Martín, C. L. (2015). The role of the loss function
#' in value-at-risk comparisons. The Journal of Risk Model Validation, 9(1), 1-19.
#'
#' Sarma, M., Thomas, S., & Shah, A. (2003). Selection of Value-at-Risk models.
#' Journal of Forecasting, 22(4), 337-358.
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
#' loss <- -results$xout
#' ES <- results$ES
#' loss.data <- list(loss = loss, ES = ES)
#' lossfun(loss.data)
#'
#'# directly passing the output object of 'rollcast()' to 'lossfun()'
#' lossfun(results)
#'

lossfun <- function(obj = list(loss = NULL, ES = NULL), beta = 1e-04) {

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

  if (!inherits(obj, "quarks") && (length(obj[["ES"]]) <= 1 ||
                                   any(is.na(obj[["ES"]])) ||
                                   !is.numeric(obj[["ES"]]))) {
    stop("A numeric vector of length > 1 and without NAs must be passed to",
         " 'obj$ES'.")
  }

  if (length(beta) != 1 || is.na(beta) || !is.numeric(beta)) {
    stop("A single numeric value must be passed to"," 'beta'")
  }

  if(inherits(obj, "quarks")) {
    loss <- -obj[["xout"]]
  }
  else {
    loss <- obj[["loss"]]
  }
  ES <- obj[["ES"]]

  loss.ed <- (loss - ES)[loss > ES]
  loss.ex2 <- (loss - ES)[loss <= ES]
  loss.ex1 <- ES[loss <= ES]
  loss.ex23 <- (loss - ES)[loss <= ES & loss >= 0]
  loss.ex3 <- ES[loss <= ES & loss < 0]
  RLF <- sum(loss.ed^2)

  FLF2 <- sum(abs(loss.ex2)) * beta
  FLF23 <- sum(abs(loss.ex23)) * beta
  FLF3 <- sum(loss.ex3) * beta
  FLF1 <- sum(loss.ex1) * beta

  loss.fun1 <- RLF * 10000
  loss.fun2 <- (RLF + FLF1) * 10000
  loss.fun3 <- (RLF + FLF2) * 10000
  loss.fun4 <- (RLF + FLF23 + FLF3) * 10000
  result <- list(lossfun1 = loss.fun1,
                 lossfun2 = loss.fun2,
                 lossfun3 = loss.fun3,
                 lossfun4 = loss.fun4)
  message("\n", "Please note that the following results are multiplied with 10000.", "\n")
  return(result)
}
