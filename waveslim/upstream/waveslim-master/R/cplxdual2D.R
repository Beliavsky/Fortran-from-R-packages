#' Dual-tree Complex 2D Discrete Wavelet Transform
#' 
#' Dual-tree complex 2D discrete wavelet transform (DWT).
#' 
#' 
#' @usage cplxdual2D(x, J, Faf, af)
#' @usage icplxdual2D(w, J, Fsf, sf)
#' @aliases cplxdual2D icplxdual2D
#' @param x 2D array.
#' @param w wavelet coefficients.
#' @param J number of stages.
#' @param Faf first stage analysis filters for tree i.
#' @param af analysis filters for the remaining stages on tree i.
#' @param Fsf last stage synthesis filters for tree i.
#' @param sf synthesis filters for the preceeding stages.
#' @return For the analysis of \code{x}, the output is \item{w}{wavelet
#' coefficients indexed by \code{[[j]][[i]][[d1]][[d2]]}, where
#' \eqn{j=1,\ldots,J} (scale), \eqn{i=1} (real part) or \eqn{i=2} (imag part),
#' \eqn{d1=1,2} and \eqn{d2=1,2,3} (orientations).} For the synthesis of
#' \code{w}, the output is \item{y}{output signal.}
#' @author Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher
#' @seealso \code{\link{FSfarras}}, \code{\link{farras}}, \code{\link{afb2D}},
#' \code{\link{sfb2D}}.
#' @keywords ts
#' @examples
#' 
#' \dontrun{
#' ## EXAMPLE: cplxdual2D
#' x = matrix(rnorm(32*32), 32, 32)
#' J = 5
#' Faf = FSfarras()$af
#' Fsf = FSfarras()$sf
#' af = dualfilt1()$af
#' sf = dualfilt1()$sf
#' w = cplxdual2D(x, J, Faf, af)
#' y = icplxdual2D(w, J, Fsf, sf)
#' err = x - y
#' max(abs(err))
#' }
#' 
cplxdual2D <- function(x, J, Faf, af) {

  ## Dual-Tree Complex 2D Discrete Wavelet Transform
  ##
  ## USAGE:
  ##   w = cplxdual2D(x, J, Faf, af)
  ## INPUT:
  ##   x - 2-D array
  ##   J - number of stages
  ##   Faf{i}: first stage filters for tree i
  ##   af{i}:  filters for remaining stages on tree i
  ## OUTPUT:
  ##   w{j}{i}{d1}{d2} - wavelet coefficients
  ##       j = 1..J (scale)
  ##       i = 1 (real part); i = 2 (imag part)
  ##       d1 = 1,2; d2 = 1,2,3 (orientations)
  ##   w{J+1}{m}{n} - lowpass coefficients
  ##       d1 = 1,2; d2 = 1,2 
  ## EXAMPLE:
  ##   x = rand(256);
  ##   J = 5;
  ##   [Faf, Fsf] = FSfarras;
  ##   [af, sf] = dualfilt1;
  ##   w = cplxdual2D(x, J, Faf, af);
  ##   y = icplxdual2D(w, J, Fsf, sf);
  ##   err = x - y;
  ##   max(max(abs(err)))
  ##
  ## WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
  ## http://eeweb.poly.edu/iselesni/WaveletSoftware/

  ## normalization
  x <- x/2
  w <- vector("list", J+1)
  
  for (m in 1:2) {
    w[[1]][[m]] <- vector("list", 2)
    for (n in 1:2) {
      w[[1]][[m]][[n]] <- vector("list", 2)
      temp <- afb2D(x, Faf[[m]], Faf[[n]])
      lo <- temp$lo
      w[[1]][[m]][[n]] <- temp$hi
      if (J > 1) {
        for (j in 2:J) {
          temp <- afb2D(lo, af[[m]], af[[n]])
          lo <- temp$lo
          w[[j]][[m]][[n]] <- temp$hi
        }
        w[[J+1]][[m]][[n]] <- lo
      }
    }
  }

  for (j in 1:J) {
    for (m in 1:3) {
      w[[j]][[1]][[1]][[m]] <- pm(w[[j]][[1]][[1]][[m]])
      w[[j]][[2]][[2]][[m]] <- pm(w[[j]][[2]][[2]][[m]])
      w[[j]][[1]][[2]][[m]] <- pm(w[[j]][[1]][[2]][[m]])
      w[[j]][[2]][[1]][[m]] <- pm(w[[j]][[2]][[1]][[m]])
    }
  }
  return(w)
}

icplxdual2D <- function(w, J, Fsf, sf) {

  ## Inverse Dual-Tree Complex 2D Discrete Wavelet Transform
  ## 
  ## USAGE:
  ##   y = icplxdual2D(w, J, Fsf, sf)
  ## INPUT:
  ##   w - wavelet coefficients
  ##   J - number of stages
  ##   Fsf - synthesis filters for final stage
  ##   sf - synthesis filters for preceeding stages
  ## OUTPUT:
  ##   y - output array
  ## See cplxdual2D
  ##
  ## WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
  ## http://eeweb.poly.edu/iselesni/WaveletSoftware/

  for (j in 1:J) {
    for (m in 1:3) {
      w[[j]][[1]][[1]][[m]] <- pm(w[[j]][[1]][[1]][[m]])
      w[[j]][[2]][[2]][[m]] <- pm(w[[j]][[2]][[2]][[m]])
      w[[j]][[1]][[2]][[m]] <- pm(w[[j]][[1]][[2]][[m]])
      w[[j]][[2]][[1]][[m]] <- pm(w[[j]][[2]][[1]][[m]])
    }
  }

  y <- matrix(0, 2*nrow(w[[1]][[1]][[1]][[1]]), 2*ncol(w[[1]][[1]][[1]][[1]]))
  for (m in 1:2) {
    for (n in 1:2) {
      lo <- w[[J+1]][[m]][[n]]
      if (J > 1) {
        for (j in J:2) {
          lo <- sfb2D(lo, w[[j]][[m]][[n]], sf[[m]], sf[[n]])
        }
        lo <- sfb2D(lo, w[[1]][[m]][[n]], Fsf[[m]], Fsf[[n]])
        y <- y + lo
      }
    }
  }

  ## normalization
  return(y/2)
}

