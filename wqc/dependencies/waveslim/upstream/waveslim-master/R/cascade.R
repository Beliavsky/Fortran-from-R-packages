#' Higher-Order Wavelet Filters
#' 
#' Create a wavelet filter at arbitrary scale.
#' 
#' Uses \code{cascade} subroutine to compute higher-order wavelet coefficient
#' vector from a given filtering sequence.
#' 
#' @param wf.name Character string of wavelet filter.
#' @param filter.seq Character string of filter sequence.  \code{H} means
#' high-pass filtering and \code{L} means low-pass filtering.  Sequence is read
#' from right to left.
#' @param n Length of zero-padded filter.  Frequency resolution will be
#' \code{n}/2+1.
#' @return Vector of wavelet coefficients.
#' @author B. Whitcher
#' @seealso \code{\link{squared.gain}}, \code{\link{wave.filter}}.
#' @references Bruce, A. and H.-Y. Gao (1996).  \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' 
#' Doroslovacki, M. L. (1998) On the least asymmetric wavelets, \emph{IEEE
#' Transactions on Signal Processing}, \bold{46}, No. 4, 1125-1130.
#' 
#' Daubechies, I. (1992) \emph{Ten Lectures on Wavelets}, CBMS-NSF Regional
#' Conference Series in Applied Mathematics, SIAM: Philadelphia.
#' 
#' Morris and Peravali (1999) Minimum-bandwidth discrete-time wavelets,
#' \emph{Signal Processing}, \bold{76}, No. 2, 181-193.
#' 
#' Nielsen, M. (2001) On the Construction and Frequency Localization of Finite
#' Orthogonal Quadrature Filters, \emph{Journal of Approximation Theory},
#' \bold{108}, No. 1, 36-52.
#' @keywords ts
#' @examples
#' 
#' ## Figure 4.14 in Gencay, Selcuk and Whitcher (2001)
#' par(mfrow=c(3,1), mar=c(5-2,4,4-1,2))
#' f.seq <- "HLLLLL"
#' plot(c(rep(0,33), wavelet.filter("mb4", f.seq), rep(0,33)), type="l",
#'      xlab="", ylab="", main="D(4) in black, MB(4) in red")
#' lines(c(rep(0,33), wavelet.filter("d4", f.seq), rep(0,33)), col=2)
#' plot(c(rep(0,35), -wavelet.filter("mb8", f.seq), rep(0,35)), type="l",
#'      xlab="", ylab="", main="D(8) in black, -MB(8) in red")
#' lines(c(rep(0,35), wavelet.filter("d8", f.seq), rep(0,35)), col=2)
#' plot(c(rep(0,39), wavelet.filter("mb16", f.seq), rep(0,39)), type="l",
#'      xlab="", ylab="", main="D(16) in black, MB(16) in red")
#' lines(c(rep(0,39), wavelet.filter("d16", f.seq), rep(0,39)), col=2)
#' 
#' @export wavelet.filter
wavelet.filter <- function(wf.name, filter.seq = "L", n = 512)
{
  cascade <- function(f, x, j)
    {
      L <- length(f)
      N <- length(x)
      M <- (L - 1) * 2^j
      M1 <- M - L + 2
      M2 <- 2 * M - L + 2
      if(N > M1)
        stop("x is too long\n")
      else x <- c(x, rep(0, M1 - N))
      xj <- c(rep(0, M), x, rep(0, M))
      yj <- rep(0, M2)
      for(i in 1:L)
        yj <- yj + f[L - i + 1] * xj[1:M2 + (i - 1) * 2^j]
      yj
    }
  if(is.character(wf.name))
    wf <- wave.filter(wf.name)
  else
    wf <- wf.name
  J <- nchar(filter.seq)
  key <- rev(substring(filter.seq, 1:J, 1:J))
  f <- 1
  fl <- wf$lpf
  fh <- wf$hpf
  for(k in 1:J) {
    if(key[k] == "H")
      f <- cascade(fh, f, k - 1)
    else if(key[k] == "L")
      f <- cascade(fl, f, k - 1)
    else stop("Invalid filter.seq\n")
  }
  f
}



#' Squared Gain Function of a Filter
#' 
#' Produces the modulus squared of the Fourier transform for a given filtering
#' sequence.
#' 
#' Uses \code{cascade} subroutine to compute the squared gain function from a
#' given filtering sequence.
#' 
#' @param wf.name Character string of wavelet filter.
#' @param filter.seq Character string of filter sequence.  \code{H} means
#' high-pass filtering and \code{L} means low-pass filtering.  Sequence is read
#' from right to left.
#' @param n Length of zero-padded filter.  Frequency resolution will be
#' \code{n}/2+1.
#' @return Squared gain function.
#' @author B. Whitcher
#' @seealso \code{\link{wave.filter}}, \code{\link{wavelet.filter}}.
#' @keywords ts
#' @examples
#' 
#' par(mfrow=c(2,2))
#' f.seq <- "H"
#' plot(0:256/512, squared.gain("d4", f.seq), type="l", ylim=c(0,2),
#'      xlab="frequency", ylab="L = 4", main="Level 1")
#' lines(0:256/512, squared.gain("fk4", f.seq), col=2)
#' lines(0:256/512, squared.gain("mb4", f.seq), col=3)
#' abline(v=c(1,2)/4, lty=2)
#' legend(-.02, 2, c("Daubechies", "Fejer-Korovkin", "Minimum-Bandwidth"),
#'        lty=1, col=1:3, bty="n", cex=1)
#' f.seq <- "HL"
#' plot(0:256/512, squared.gain("d4", f.seq), type="l", ylim=c(0,4),
#'      xlab="frequency", ylab="", main="Level 2")
#' lines(0:256/512, squared.gain("fk4", f.seq), col=2)
#' lines(0:256/512, squared.gain("mb4", f.seq), col=3)
#' abline(v=c(1,2)/8, lty=2)
#' f.seq <- "H"
#' plot(0:256/512, squared.gain("d8", f.seq), type="l", ylim=c(0,2),
#'      xlab="frequency", ylab="L = 8", main="")
#' lines(0:256/512, squared.gain("fk8", f.seq), col=2)
#' lines(0:256/512, squared.gain("mb8", f.seq), col=3)
#' abline(v=c(1,2)/4, lty=2)
#' f.seq <- "HL"
#' plot(0:256/512, squared.gain("d8", f.seq), type="l", ylim=c(0,4),
#'      xlab="frequency", ylab="", main="")
#' lines(0:256/512, squared.gain("fk8", f.seq), col=2)
#' lines(0:256/512, squared.gain("mb8", f.seq), col=3)
#' abline(v=c(1,2)/8, lty=2)
#' 
#' @export squared.gain
squared.gain <- function(wf.name, filter.seq = "L", n = 512)
{
  cascade <- function(f, x, j)
    {
      L <- length(f)
      N <- length(x)
      M <- (L - 1) * 2^j
      M1 <- M - L + 2
      M2 <- 2 * M - L + 2
      if(N > M1)
        stop("x is too long\n")
      else x <- c(x, rep(0, M1 - N))
      xj <- c(rep(0, M), x, rep(0, M))
      yj <- rep(0, M2)
      for(i in 1:L)
        yj <- yj + f[L - i + 1] * xj[1:M2 + (i - 1) * 2^j]
      yj
    }
  if(is.character(wf.name))
    wf <- wave.filter(wf.name)
  else 
    wf <- wf.name
  J <- nchar(filter.seq)
  key <- rev(substring(filter.seq, 1:J, 1:J))
  f <- 1
  fl <- wf$lpf
  fh <- wf$hpf
  for(k in 1:J) {
    if(key[k] == "H")
      f <- cascade(fh, f, k - 1)
    else if(key[k] == "L")
      f <- cascade(fl, f, k - 1)
    else stop("Invalid filter.seq\n")
  }
  Mod(fft(c(f, rep(0, n - length(f))))[1:(n/2 + 1)])^2
}
