#' Generate Stationary Gaussian Process Using Hosking's Method
#' 
#' Uses exact time-domain method from Hosking (1984) to generate a simulated
#' time series from a specified autocovariance sequence.
#' 
#' 
#' @param n Length of series.
#' @param acvs Autocovariance sequence of series with which to generate, must
#' be of length at least \code{n}.
#' @return Length \code{n} time series from true autocovariance sequence
#' \code{acvs}.
#' @author B. Whitcher
#' @references Hosking, J. R. M. (1984) Modeling persistence in hydrological
#' time series using fractional differencing, \emph{Water Resources Research},
#' \bold{20}, No. 12, 1898-1908.
#' 
#' Percival, D. B. (1992) Simulating Gaussian random processes with specified
#' spectra, \emph{Computing Science and Statistics}, \bold{22}, 534-538.
#' @keywords ts
#' @examples
#' 
#' dB <- function(x) 10 * log10(x)
#' per <- function (z) {
#'   n <- length(z)
#'   (Mod(fft(z))^2/(2 * pi * n))[1:(n%/%2 + 1)]
#' }
#' spp.sdf <- function(freq, delta, omega)
#'   abs(2 * (cos(2*pi*freq) - cos(2*pi*omega)))^(-2*delta)
#' data(acvs.andel8)
#' n <- 1024
#' \dontrun{
#' z <- hosking.sim(n, acvs.andel8[,2])
#' per.z <- 2 * pi * per(z)
#' par(mfrow=c(2,1), las=1)
#' plot.ts(z, ylab="", main="Realization of a Seasonal Long-Memory Process")
#' plot(0:(n/2)/n, dB(per.z), type="l", xlab="Frequency", ylab="dB",
#'      main="Periodogram")
#' lines(0:(n/2)/n, dB(spp.sdf(0:(n/2)/n, .4, 1/12)), col=2)
#' }
#' 
#' @export hosking.sim
hosking.sim <- function(n, acvs) {
  .C(C_hosking, tseries=rnorm(n), as.integer(n), as.double(acvs[1:n]))$tseries
}
