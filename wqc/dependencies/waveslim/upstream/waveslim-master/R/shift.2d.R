#' Circularly Shift Matrices from a 2D MODWT
#' 
#' Compute phase shifts for wavelet sub-matrices based on the ``center of
#' energy'' argument of Hess-Nielsen and Wickerhauser (1996).
#' 
#' The "center of energy" technique of Wickerhauser and Hess-Nielsen (1996) is
#' employed to find circular shifts for the wavelet sub-matrices such that the
#' coefficients are aligned with the original series.  This corresponds to
#' applying a (near) linear-phase filtering operation.
#' 
#' @param z Two-dimensional MODWT object
#' @param inverse Boolean value on whether to perform the forward or inverse
#' operation.
#' @return Two-dimensional MODWT object with circularly shifted coefficients.
#' @author B. Whitcher
#' @seealso \code{\link{phase.shift}}, \code{\link{modwt.2d}}.
#' @references Hess-Nielsen, N. and M. V. Wickerhauser (1996) Wavelets and
#' time-frequency analysis, \emph{Proceedings of the IEEE}, \bold{84}, No. 4,
#' 523-540.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @examples
#' 
#' n <- 512
#' G1 <- G2 <- dnorm(seq(-n/4, n/4, length=n))
#' G <- 100 * zapsmall(outer(G1, G2))
#' G <- modwt.2d(G, wf="la8", J=6)
#' k <- 50
#' xr <- yr <- trunc(n/2) + (-k:k)
#' par(mfrow=c(3,3), mar=c(1,1,2,1), pty="s")
#' for (j in names(G)[1:9]) {
#'   image(G[[j]][xr,yr], col=rainbow(64), axes=FALSE, main=j)
#' }
#' Gs <- shift.2d(G)
#' for (j in names(G)[1:9]) {
#'   image(Gs[[j]][xr,yr], col=rainbow(64), axes=FALSE, main=j)
#' }
#' 
#' @export shift.2d
shift.2d <- function(z, inverse=FALSE) {
  ## "Center of Energy"
  coe <- function(g) { 
    sum(0:(length(g)-1) * g^2) / sum(g^2)
  }

  wf <- attributes(z)$wavelet
  h <- wave.filter(wf)$hpf
  g <- wave.filter(wf)$lpf
  
  J <- (length(z) - 1) / 3
  m <- nrow(z[[1]])
  n <- ncol(z[[1]])
  
  nu.H <- round(2^(1:J-1) * (coe(g) + coe(h)) - coe(g), 0)
  nu.Hm <- ifelse(nu.H/m < 1, nu.H, nu.H - trunc(nu.H/m) * m)
  nu.Hn <- ifelse(nu.H/n < 1, nu.H, nu.H - trunc(nu.H/n) * n)
  nu.G <- round((2^(1:J) - 1) * coe(g), 0)
  nu.Gm <- ifelse(nu.G/m < 1, nu.G, nu.G - trunc(nu.G/m) * m)
  nu.Gn <- ifelse(nu.G/n < 1, nu.G, nu.G - trunc(nu.G/n) * n)
  
  if (!inverse) {
    ## Apply the phase shifts
    for (j in 0:(J-1)) {
      Hm.order <- c((nu.H[j+1]+1):m, 1:nu.H[j+1])
      Hn.order <- c((nu.H[j+1]+1):n, 1:nu.H[j+1])
      Gm.order <- c((nu.G[j+1]+1):m, 1:nu.G[j+1])
      Gn.order <- c((nu.G[j+1]+1):n, 1:nu.G[j+1])
      z[[3*j+1]] <- z[[3*j+1]][Gm.order, Hn.order]
      z[[3*j+2]] <- z[[3*j+2]][Hm.order, Gn.order]
      z[[3*j+3]] <- z[[3*j+3]][Hm.order, Hn.order]
    } 
    z[[3*J+1]] <- z[[3*J+1]][Gm.order, Gn.order]
  } else {
    ## Apply the phase shifts "reversed"
    for (j in 0:(J-1)) {
      Hm.order <- c((m-nu.H[j+1]+1):m, 1:(m-nu.H[j+1]))
      Hn.order <- c((n-nu.H[j+1]+1):n, 1:(n-nu.H[j+1]))
      Gm.order <- c((m-nu.G[j+1]+1):m, 1:(m-nu.G[j+1]))
      Gn.order <- c((n-nu.G[j+1]+1):n, 1:(n-nu.G[j+1]))
      z[[3*j+1]] <- z[[3*j+1]][Gm.order, Hn.order]
      z[[3*j+2]] <- z[[3*j+2]][Hm.order, Gn.order]
      z[[3*j+3]] <- z[[3*j+3]][Hm.order, Hn.order]
    }
    z[[3*J+1]] <- z[[3*J+1]][Gm.order, Gn.order]
  }
  return(z)
}

