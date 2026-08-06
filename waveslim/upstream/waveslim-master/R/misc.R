#' Autocovariance Functions via the Discrete Fourier Transform
#' 
#' Computes the autocovariance function (ACF) for a time series or the
#' cross-covariance function (CCF) between two time series.
#' 
#' The series is zero padded to twice its length before the discrete Fourier
#' transform is applied.  Only the values corresponding to nonnegative lags are
#' provided (for the ACF).
#' 
#' @usage my.acf(x)
#' @usage my.ccf(a, b)
#' @aliases my.acf my.ccf
#' @param x,a,b time series
#' @return The autocovariance function for all nonnegative lags or the
#' cross-covariance function for all lags.
#' @author B. Whitcher
#' @keywords ts
#' @examples
#' 
#' data(ibm)
#' ibm.returns <- diff(log(ibm))
#' plot(1:length(ibm.returns) - 1, my.acf(ibm.returns), type="h",
#'      xlab="lag", ylab="ACVS", main="Autocovariance Sequence for IBM Returns")
#' 
#' @export my.acf
my.acf <- function(x)
{
  n <- length(x)
  x <- c(x, rep(0, n))
  Re(fft(Mod(fft(x)) ^ 2, inverse = TRUE) / 2 / n ^ 2)[1:n]
}

my.ccf <- function(a, b) {
  n <- length(a)
  a <- c(a, rep(0, n))
  b <- c(b, rep(0, n))
  x <- Re(fft(fft(a) * Conj(fft(b)), inverse = TRUE)) / 2 / n ^ 2
  x[c((n + 2):(2 * n), 1:n)]
}
