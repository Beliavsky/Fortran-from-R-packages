#' Kingsbury's Q-filters for the Dual-Tree Complex DWT
#' 
#' Kingsbury's Q-filters for the dual-tree complex DWT.
#' 
#' These cofficients are rounded to 8 decimal places.
#' 
#' @aliases dualfilt1 AntonB
#' @return \item{af}{List (i=1,2) - analysis filters for tree i} \item{sf}{List
#' (i=1,2) - synthesis filters for tree i} Note: \code{af[[2]]} is the reverse
#' of \code{af[[1]]}.
#' @author Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher
#' @seealso \code{\link{dualtree}}
#' @references Kingsbury, N.G. (2000).  A dual-tree complex wavelet transform
#' with improved orthogonality and symmetry properties, \emph{Proceedings of
#' the IEEE Int. Conf. on Image Proc.} (ICIP).
#' @keywords ts
#' @export dualfilt1
dualfilt1 <- function() {

  af1 <- c(0.03516384000000, 0,
           0, 0,
           -0.08832942000000, -0.11430184000000,
           0.23389032000000, 0,
           0.76027237000000, 0.58751830000000,
           0.58751830000000, -0.76027237000000,
           0, 0.23389032000000,
           -0.11430184000000, 0.08832942000000,
           0, 0,
           0, -0.03516384000000)
  af1 <- matrix(af1, ncol=2, byrow=TRUE)
  af2 <- c(0, -0.03516384000000,
           0, 0,
           -0.11430184000000, 0.08832942000000,
           0, 0.23389032000000,
           0.58751830000000, -0.76027237000000,
           0.76027237000000, 0.58751830000000,
           0.23389032000000, 0,
           -0.08832942000000, -0.11430184000000,
           0, 0,
           0.03516384000000, 0)
  af2 <- matrix(af2, ncol=2, byrow=TRUE)
  sf1 <- af1[nrow(af1):1, ]
  sf2 <- af2[nrow(af2):1, ]
  list(af = list(af1, af2), sf = list(sf1, sf2))
}



#' Farras nearly symmetric filters
#' 
#' Farras nearly symmetric filters for orthogonal 2-channel perfect
#' reconstruction filter bank and Farras filters organized for the dual-tree
#' complex DWT.
#' 
#' 
#' @aliases farras FSfarras
#' @return \item{af}{List (i=1,2) - analysis filters for tree i} \item{sf}{List
#' (i=1,2) - synthesis filters for tree i}
#' @author Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher
#' @seealso \code{\link{afb}}, \code{\link{dualtree}}, \code{\link{dualfilt1}}.
#' @references A. F. Abdelnour and I. W. Selesnick. \dQuote{Nearly symmetric
#' orthogonal wavelet bases}, Proc. IEEE Int. Conf. Acoust., Speech, Signal
#' Processing (ICASSP), May 2001.
#' @keywords ts
FSfarras <- function() {

  af1 <- c(0, 0,
           -0.08838834764832,  -0.01122679215254,
           0.08838834764832,   0.01122679215254,
           0.69587998903400,   0.08838834764832,
           0.69587998903400,   0.08838834764832,
           0.08838834764832,  -0.69587998903400,
           -0.08838834764832,   0.69587998903400,
           0.01122679215254,  -0.08838834764832,
           0.01122679215254,  -0.08838834764832,
           0, 0)
  af1 <- matrix(af1, ncol=2, byrow=TRUE)
  sf1 <- af1[nrow(af1):1, ]
  af2 <- c(0.01122679215254, 0,
           0.01122679215254, 0,
           -0.08838834764832, -0.08838834764832,
           0.08838834764832, -0.08838834764832,
           0.69587998903400, 0.69587998903400,
           0.69587998903400, -0.69587998903400,
           0.08838834764832, 0.08838834764832,
           -0.08838834764832, 0.08838834764832,
           0, 0.01122679215254,
           0, -0.01122679215254)
  af2 <- matrix(af2, ncol=2, byrow=TRUE)
  sf2 <- af2[nrow(af2):1, ]
  list(af = list(af1, af2), sf = list(sf1, sf2))
}

farras <- function() {

  af <- c(0, -0.01122679215254, 0, 0.01122679215254,
          -0.08838834764832, 0.08838834764832,
          0.08838834764832, 0.08838834764832,
          0.69587998903400, -0.69587998903400,
          0.69587998903400, 0.69587998903400,
          0.08838834764832, -0.08838834764832,
          -0.08838834764832, -0.08838834764832,
          0.01122679215254, 0, 0.01122679215254, 0)
  af <- matrix(af, nrow=10, byrow=TRUE)
  sf <- af[nrow(af):1, ]
  list(af = af, sf = sf)
}



#' Miscellaneous Functions for Dual-Tree Wavelet Software
#' 
#' Miscellaneous functions for dual-tree wavelet software.
#' 
#' 
#' @usage cshift(x, m)
#' @usage cshift2D(x, m)
#' @usage pm(a, b)
#' @aliases cshift cshift2D pm
#' @param x N-point vector
#' @param m amount of shift
#' @param a,b input parameters
#' @return \item{y}{vector \code{x} will be shifed by \code{m} samples to the
#' left or matrix \code{x} will be shifed by \code{m} samples down.}
#' \item{u}{(a + b) / sqrt(2)} \item{v}{(a - b) / sqrt(2)}
#' @author Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher
#' @keywords ts
cshift <- function(x, m) {

  N <- length(x)
  n <- 0:(N-1)
  n <- (n-m) %% N
  y <- x[n+1]
  y
}



#' Filter Banks for Dual-Tree Wavelet Transforms
#' 
#' Analysis and synthesis filter banks used in dual-tree wavelet algorithms.
#' 
#' The functions \code{afb2D.A} and \code{sfb2D.A} implement the convolutions,
#' either for analysis or synthesis, in one dimension only.  Thus, they are the
#' workhorses of \code{afb2D} and \code{sfb2D}.  The output for the analysis
#' filter bank along one dimension (\code{afb2D.A}) is a list with two elements
#' \describe{ \item{lo}{low-pass subband} \item{hi}{high-pass subband} } where
#' the dimension of analysis will be half its original length.  The output for
#' the synthesis filter bank along one dimension (\code{sfb2D.A}) will be the
#' output array, where the dimension of synthesis will be twice its original
#' length.
#' 
#' @usage afb(x, af)
#' @usage afb2D(x, af1, af2 = NULL)
#' @usage afb2D.A(x, af, d)
#' @usage sfb(lo, hi, sf)
#' @usage sfb2D(lo, hi, sf1, sf2 = NULL)
#' @usage sfb2D.A(lo, hi, sf, d)
#' @aliases afb afb2D afb2D.A sfb sfb2D sfb2D.A
#' @param x vector or matrix of observations
#' @param af analysis filters.  First element of the list is the low-pass
#' filter, second element is the high-pass filter.
#' @param af1,af2 analysis filters for the first and second dimension of a 2D
#' array.
#' @param sf synthesis filters.  First element of the list is the low-pass
#' filter, second element is the high-pass filter.
#' @param sf1,sf2 synthesis filters for the first and second dimension of a 2D
#' array.
#' @param d dimension of filtering (d = 1 or 2)
#' @param lo low-frequecy coefficients
#' @param hi high-frequency coefficients
#' @return In one dimension the output for the analysis filter bank
#' (\code{afb}) is a list with two elements \item{lo}{Low frequecy output}
#' \item{hi}{High frequency output} and the output for the synthesis filter
#' bank (\code{sfb}) is the output signal.
#' 
#' In two dimensions the output for the analysis filter bank (\code{afb2D}) is
#' a list with four elements \item{lo}{low-pass subband} \item{hi[[1]]}{'lohi'
#' subband} \item{hi[[2]]}{'hilo' subband} \item{hi[[3]]}{'hihi' subband} and
#' the output for the synthesis filter bank (\code{sfb2D}) is the output array.
#' @author Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher
#' @keywords ts
#' @examples
#' 
#' ## EXAMPLE: afb, sfb
#' af = farras()$af
#' sf = farras()$sf
#' x = rnorm(64)
#' x.afb = afb(x, af)
#' lo = x.afb$lo
#' hi = x.afb$hi
#' y = sfb(lo, hi, sf)
#' err = x - y
#' max(abs(err))
#' 
#' ## EXAMPLE: afb2D, sfb2D
#' x = matrix(rnorm(32*64), 32, 64)
#' af = farras()$af
#' sf = farras()$sf
#' x.afb2D = afb2D(x, af, af)
#' lo = x.afb2D$lo
#' hi = x.afb2D$hi
#' y = sfb2D(lo, hi, sf, sf)
#' err = x - y
#' max(abs(err))
#' 
#' ## Example: afb2D.A, sfb2D.A
#' x = matrix(rnorm(32*64), 32, 64)
#' af = farras()$af
#' sf = farras()$sf
#' x.afb2D.A = afb2D.A(x, af, 1)
#' lo = x.afb2D.A$lo
#' hi = x.afb2D.A$hi
#' y = sfb2D.A(lo, hi, sf, 1)
#' err = x - y
#' max(abs(err))
afb <- function(x, af) {

  N <- length(x)
  L <- nrow(af)/2
  x <- cshift(x,-L)
  
  ## lowpass filter
  lo <- convolve(x, af[,1], conj=FALSE, type="open")
  lo <- cshift(lo,-(2*L-1))
  lo <- lo[seq(1, length(lo), by=2)]
  lo[1:L] <- lo[N/2+(1:L)] + lo[1:L]
  lo <- lo[1:(N/2)]
  
  ## highpass filter
  hi <- convolve(x, af[,2], conj=FALSE, type="open")
  hi <- cshift(hi,-(2*L-1))
  hi <- hi[seq(1, length(hi), by=2)]
  hi[1:L] <- hi[N/2+(1:L)] + hi[1:L]
  hi <- hi[1:(N/2)]

  list(lo = lo, hi = hi)
}



#' Dual-tree Complex Discrete Wavelet Transform
#' 
#' One- and two-dimensional dual-tree complex discrete wavelet transforms
#' developed by Kingsbury and Selesnick \emph{et al.}
#' 
#' In one dimension \eqn{N} is divisible by \eqn{2^J} and
#' \eqn{N\ge2^{J-1}\cdot\mbox{length}(\mbox{\code{af}})}.
#' 
#' In two dimensions, these two conditions must hold for both \eqn{M} and
#' \eqn{N}.
#' 
#' @usage dualtree(x, J, Faf, af)
#' @usage idualtree(w, J, Fsf, sf)
#' @usage dualtree2D(x, J, Faf, af)
#' @usage idualtree2D(w, J, Fsf, sf)
#' @aliases dualtree idualtree dualtree2D idualtree2D
#' @param x N-point vector or MxN matrix.
#' @param w DWT coefficients.
#' @param J number of stages.
#' @param Faf analysis filters for the first stage.
#' @param af analysis filters for the remaining stages.
#' @param Fsf synthesis filters for the last stage.
#' @param sf synthesis filters for the preceeding stages.
#' @return For the analysis of \code{x}, the output is \item{w}{DWT
#' coefficients.  Each wavelet scale is a list containing the real and
#' imaginary parts.  The final scale (J+1) contains the low-pass filter
#' coefficients.} For the synthesis of \code{w}, the output is \item{y}{output
#' signal}
#' @author Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher
#' @seealso \code{\link{FSfarras}}, \code{\link{farras}},
#' \code{\link{convolve}}, \code{\link{cshift}}, \code{\link{afb}},
#' \code{\link{sfb}}.
#' @keywords ts
#' @examples
#' 
#' ## EXAMPLE: dualtree
#' x = rnorm(512)
#' J = 4
#' Faf = FSfarras()$af
#' Fsf = FSfarras()$sf
#' af = dualfilt1()$af
#' sf = dualfilt1()$sf
#' w = dualtree(x, J, Faf, af)
#' y = idualtree(w, J, Fsf, sf)
#' err = x - y
#' max(abs(err))
#' 
#' ## Example: dualtree2D
#' x = matrix(rnorm(64*64), 64, 64)
#' J = 3
#' Faf = FSfarras()$af
#' Fsf = FSfarras()$sf
#' af = dualfilt1()$af
#' sf = dualfilt1()$sf
#' w = dualtree2D(x, J, Faf, af)
#' y = idualtree2D(w, J, Fsf, sf)
#' err = x - y
#' max(abs(err))
#' 
#' ## Display 2D wavelets of dualtree2D.m
#' 
#' J <- 4
#' L <- 3 * 2^(J+1)
#' N <- L / 2^J
#' Faf <- FSfarras()$af
#' Fsf <- FSfarras()$sf
#' af <- dualfilt1()$af
#' sf <- dualfilt1()$sf
#' x <- matrix(0, 2*L, 3*L)
#' w <- dualtree2D(x, J, Faf, af)
#' w[[J]][[1]][[1]][N/2, N/2+0*N] <- 1
#' w[[J]][[1]][[2]][N/2, N/2+1*N] <- 1
#' w[[J]][[1]][[3]][N/2, N/2+2*N] <- 1
#' w[[J]][[2]][[1]][N/2+N, N/2+0*N] <- 1
#' w[[J]][[2]][[2]][N/2+N, N/2+1*N] <- 1
#' w[[J]][[2]][[3]][N/2+N, N/2+2*N] <- 1
#' y <- idualtree2D(w, J, Fsf, sf)
#' image(t(y), col=grey(0:64/64), axes=FALSE)
#' 
dualtree <- function(x, J, Faf, af) {

  ## normalization
  x <- x/sqrt(2)
  w <- vector("list", J+1)

  ## Tree 1
  w[[1]] <- vector("list", 2)
  temp <- afb(x, Faf[[1]])
  x1 <- temp$lo
  w[[1]][[1]] <- temp$hi
  if(J > 1) {
    for(j in 2:J) {
      w[[j]] <- vector("list", 2)
      temp <- afb(x1, af[[1]])
      x1 <- temp$lo
      w[[j]][[1]] <- temp$hi
    }
  }
  w[[J+1]] <- vector("list", 2)
  w[[J+1]][[1]] <- x1

  ## Tree 2
  temp <- afb(x, Faf[[2]])
  x2 <- temp$lo
  w[[1]][[2]] <- temp$hi
  if(J > 1) {
    for(j in 2:J) {
      temp <- afb(x2, af[[2]])
      x2 <- temp$lo
      w[[j]][[2]] <- temp$hi
    }
  }
  w[[J+1]][[2]] <- x2
  w
}

sfb <- function(lo, hi, sf) {

  N <- 2*length(lo)
  L <- nrow(sf)
  ## lo <- upfirdn(lo, sf[,1], 2, 1)
  lo <- c(matrix(c(rep(0, N/2), lo), nrow=2, byrow=TRUE))
  lo <- convolve(lo, sf[,1], conj=FALSE, type="open")
  lo <- cshift(lo, -L)
  ## hi <- upfirdn(hi, sf[,2], 2, 1)
  hi <- c(matrix(c(rep(0, N/2), hi), nrow=2, byrow=TRUE))
  hi <- convolve(hi, sf[,2], conj=FALSE, type="open")
  hi <- cshift(hi, -L)
  
  y <- lo + hi
  y[1:(L-2)] <- y[1:(L-2)] + y[N+1:(L-2)]
  y <- y[1:N]
  ## y = cshift(y, 1-L/2);
  y <- cshift(y, 1-L/2)
  y
}

idualtree <- function(w, J, Fsf, sf) {

  ## Tree 1
  y1 <- w[[J+1]][[1]]
  if(J > 1) {
    for(j in J:2) {
      y1 <- sfb(y1, w[[j]][[1]], sf[[1]])
    }
  }
  y1 <- sfb(y1, w[[1]][[1]], Fsf[[1]])
  
  ## Tree 2
  y2 <- w[[J+1]][[2]]
  if(J > 1) {
    for(j in J:2) {
      y2 <- sfb(y2, w[[j]][[2]], sf[[2]])
    }
  }
  y2 <- sfb(y2, w[[1]][[2]], Fsf[[2]])

  ## normalization
  y <- (y1 + y2)/sqrt(2)
  y
}

