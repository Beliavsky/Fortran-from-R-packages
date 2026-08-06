#' Periodogram
#' 
#' Computation of the periodogram via the Fast Fourier Transform (FFT).
#' 
#' 
#' @param z time series
#' @author Author: Jan Beran; modified: Martin Maechler, Date: Sep 1995.
#' @keywords ts
per <- function(z) {
  n <- length(z)
  (Mod(fft(z)) ** 2 / (2 * pi * n)) [1:(n %/% 2 + 1)]
}
