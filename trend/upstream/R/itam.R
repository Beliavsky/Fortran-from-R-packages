#' @title Innovative Trend Analysis Method
#' @description
#' The function performs the innovative trend analysis method (Şen 2012).
#'
# @param x numeric vector or a time series object of class \dQuote{ts}
# @param conf.level numeric, the level of confidence
#'
#' @inheritParams sens.slope x conf.level
#'
#' @details
#' The magnitude of trend is calculated based on the difference between
#' the arithmeric means of the second half \eqn{\bar{y}} and the first half
#' \eqn{\bar{x}} of the time series:
#'
#' \deqn{
#'  b = 2 ~ \left( \bar{y} - \bar{x} \right) / n,
#' }
#'
#' with \eqn{n} the total number of observations in the full series.
#'
#' The corresponding standard deviation is
#' \deqn{
#' \sigma_{b} = \frac{2 ~ \sqrt(2)}{n \sqrt(n)} s \sqrt(1 - \rho_{xy}),
#' }
#'
#' with \eqn{s} the standard deviation of the full series and \eqn{\rho_{xy}}
#' the correlation coeffientient between the orderer first half \eqn{(X)}
#' and the  ordered second half \eqn{(Y)} sub-series.
#'
#' The test statistitic is:
#' \deqn{
#'  z = b / \sigma_{b}.
#' }
#'
#' The p-values are calculated from the standard normal distribution for the
#' two-sided case.
#'
#' @returns
#' A list of class \dQuote{htest} and \dQuote{itam}.
#'
#'  \item{estimates}{numeric, ITAM slope}
#'  \item{data.name}{character string that denotes the input data}
#'  \item{p.value}{the p-value}
#'  \item{statistic}{the z quantile of the standard normal distribution}
#'  \item{null.value}{the null hypothesis}
#'  \item{conf.in}{upper and lower confidence limit}
#'  \item{alternative}{the alternative hypothesis}
#'  \item{method}{character string that denotes the test}
#'  \item{df}{A data.frame with ordered x1 and ordered x2}
#'
#' @references
#' Şen, Z. (2012) Innovative Trend Analysis Methodology.
#' \emph{Journal of Hydrologic Engineering} \bold{17}, 1042--1046.
#' \doi{10.1061/(ASCE)HE.1943-5584.0000556}
#'
#' @examples
#' (out <- itam(Nile))
#' plot(out)
#'
#' @note Current Version is for complete observations only.
#' For odd \code{n} the mid-point of the series is omitted.
#'
#' @keywords htest ts nonparametric
#' @importFrom stats cor sd na.fail pnorm qnorm
#' @export
itam <- function(x, conf.level = 0.95) {
  stopifnot(is.numeric(x))
  na.fail(x)
  stopifnot(0 < conf.level & conf.level < 1)
  dname <- deparse(substitute(x))

  n <- length(x)
  if (n < 4) {
    stop("sample size must be greater than 3")
  }
  l <- ceiling(n / 2)
  # leave out the middle point of the series if
  # n is not even
  if ((as.integer(n) %% 2) != 0) {
    l <- l - 1
  }
  x1 <- sort(x[1:l])
  x2 <- sort(x[(n - l + 1):n])

  # estimate
  b <- 2 * (mean(x2) - mean(x1)) / n
  # sd of slope
  sigma <- 2 * sqrt(2) / (n * sqrt(n)) * sd(x) * sqrt(1 - cor(x1, x2))

  # z-value
  z = b / sigma
  # p-value
  pval <- 2 * min(0.5, pnorm(abs(z), lower.tail = FALSE))

  # confidence interval
  ci <- qnorm((1 - conf.level) / 2, lower.tail = FALSE) * sigma
  cint <- c(b - ci, b + ci)
  attr(cint, "conf.level") <- conf.level

  ans <- list(
    method = "Innovative Trend Analysis Method",
    data.name = dname,
    estimates = c(`Sen's ITA slope` = b, SD = sigma),
    statistic = c(z = z),
    p.value = pval,
    alternative = "two.sided",
    conf.int = cint,
    null.value = c(z = 0),
    df = data.frame(x1, x2)
  )
  class(ans) <- c("itam", "htest")
  return(ans)

}

#' @title Plotting itam-objects
#' @description
#' Plotting method for objects inheriting from class \dQuote{itam}
#'
#' @param x an object of class \dQuote{itam}
#' @param \ldots further arguments, currently ignored
#'
#' @seealso
#' \code{\link{itam}}
#'
#' @inherit itam references
#' @inherit itam examples
#'
#' @keywords hplot
#' @export
plot.itam <- function(x, ...) {
  df <- x$df
  xmin <- min(c(df$x1, df$x2))
  xmax <- max(c(df$x1, df$x2))

  plot(
    x2 ~ x1,
    data = df,
    type = "p",
    xlim = c(xmin, xmax),
    ylim = c(xmin, xmax),
    xlab = "1st half",
    ylab = "2nd half",
    main = x$method,
    asp = 1
  )

  #' @importFrom graphics abline
  abline(a = 0, b = 1, col = "black")
  invisible(x)
}
