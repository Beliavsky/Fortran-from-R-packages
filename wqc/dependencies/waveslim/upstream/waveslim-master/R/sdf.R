#' Spectral Density Functions for Long-Memory Processes
#' 
#' Draws the spectral density functions (SDFs) for standard long-memory
#' processes including fractional difference (FD), seasonal persistent (SP),
#' and seasonal fractional difference (SFD) processes.
#' 
#' 
#' @usage fdp.sdf(freq, d, sigma2 = 1)
#' @usage spp.sdf(freq, d, fG, sigma2 = 1)
#' @usage spp2.sdf(freq, d1, f1, d2, f2, sigma2 = 1)
#' @usage sfd.sdf(freq, s, d, sigma2 = 1)
#' @aliases fdp.sdf spp.sdf spp2.sdf sfd.sdf
#' @param freq vector of frequencies, normally from 0 to 0.5
#' @param d,d1,d2 fractional difference parameter
#' @param fG,f1,f2 Gegenbauer frequency
#' @param s seasonal parameter
#' @param sigma2 innovations variance
#' @return The power spectrum from an FD, SP or SFD process.
#' @author B. Whitcher
#' @seealso \code{\link{fdp.mle}}, \code{\link{spp.mle}}.
#' @keywords ts
#' @examples
#' 
#' dB <- function(x) 10 * log10(x)
#' 
#' fdp.main <- expression(paste("FD", group("(",d==0.4,")")))
#' sfd.main <- expression(paste("SFD", group("(",list(s==12, d==0.4),")")))
#' spp.main <- expression(paste("SPP",
#'     group("(",list(delta==0.4, f[G]==1/12),")")))
#' 
#' freq <- 0:512/1024
#' 
#' par(mfrow=c(2,2), mar=c(5-1,4,4-1,2), col.main="darkred")
#' plot(freq, dB(fdp.sdf(freq, .4)), type="l", xlab="frequency",
#'      ylab="spectrum (dB)", main=fdp.main)
#' plot(freq, dB(spp.sdf(freq, .4, 1/12)), type="l", xlab="frequency",
#'      ylab="spectrum (dB)", font.main=1, main=spp.main)
#' plot(freq, dB(sfd.sdf(freq, 12, .4)), type="l", xlab="frequency",
#'      ylab="spectrum (dB)", main=sfd.main)
fdp.sdf <- function(freq, d, sigma2 = 1)
  sigma2 / ((2 * sin(pi * freq)) ^ 2) ^ d

bandpass.fdp <- function(a, b, d)
  2 * integrate(fdp.sdf, lower = a, upper = b, d = d)$value

spp.sdf <- function(freq, d, fG, sigma2 = 1)
  sigma2 * abs(2 * (cos(2 * pi * freq) - cos(2 * pi * fG))) ^ (-2 * d)

spp2.sdf <- function(freq, d1, f1, d2, f2, sigma2 = 1) {
  sigma2 * abs(2 * (cos(2 * pi * freq) - cos(2 * pi * f1))) ^ (-2 * d1) * 
    abs(2 * (cos(2 * pi * freq) - cos(2 * pi * f2))) ^ (-2 * d2)
}

sfd.sdf <- function(freq, s, d, sigma2=1)
  sigma2 / (2 * (1 - cos(s * 2 * pi * freq))) ^ d

bandpass.spp <- function(a, b, d, fG) {
  if (fG > a && fG < b) {
    result1 <- integrate(spp.sdf, lower=a, upper=fG, d=d, fG=fG)$value
    result2 <- integrate(spp.sdf, lower=fG, upper=b, d=d, fG=fG)$value
  }
  else {
    result1 <- integrate(spp.sdf, lower=a, upper=b, d=d, fG=fG)$value
    result2 <- 0
  }
  return(2*(result1 + result2))
}

bandpass.spp2 <- function(a, b, d1, f1, d2, f2) {
  a1 <- a
  b1 <- b
  if(a1 < f1 && b1 > f2) {
    a2 <- f1
    b2 <- f2
    result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
    result2 <- integrate(spp2.sdf, a1, b2, d1=d1, f1=f1, d2=d2, f2=f2)$value
    result3 <- integrate(spp2.sdf, b2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
  }
  else {
    if (a1 < f1 && b1 < f2) {
      a2 <- f1
      result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
      result2 <- integrate(spp2.sdf, a2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
      result3 <- 0
    }
    else {
      if (a1 < f1 && b1 > f1 && b1 < f2) {
        a2 <- f1
        result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result2 <- integrate(spp2.sdf, a2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result3 <- 0
      }
      else {
        if (a1 > f1 && a1 < f2 && b1 > f2) {
          a2 <- f2
        result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result2 <- integrate(spp2.sdf, a2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result3 <- 0
        }
        else {
          result1 <- integrate(spp2.sdf, a1, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
          result2 <- 0
          result3 <- 0
        }
      }
    }
  }
  return(2 * (result1 + result2 + result3))
}

#' Variance of a Seasonal Persistent Process
#' 
#' Computes the variance of a seasonal persistent (SP) process using a
#' hypergeometric series expansion.
#' 
#' See Lapsa (1997).  The subroutine to compute a hypergeometric series was
#' taken from \emph{Numerical Recipes in C}.
#' 
#' @usage spp.var(d, fG, sigma2 = 1)
#' @usage Hypergeometric(a, b, c, z)
#' @aliases spp.var Hypergeometric
#' @param d Fractional difference parameter.
#' @param fG Gegenbauer frequency.
#' @param sigma2 Innovations variance.
#' @param a,b,c,z Parameters for the hypergeometric series.
#' @return The variance of an SP process.
#' @author B. Whitcher
#' @references Lapsa, P.M. (1997) Determination of Gegenbauer-type random
#' process models.  \emph{Signal Processing} \bold{63}, 73-90.
#' 
#' Press, W.H., S.A. Teukolsky, W.T. Vetterling and B.P. Flannery (1992)
#' \emph{Numerical Recipes in C}, 2nd edition, Cambridge University Press.
#' @keywords ts
#' @export spp.var
spp.var <- function(d, fG, sigma2 = 1) {
  ## Hypergeometric series representation of the variance taken from
  ## Lapsa (1997)
  omega <- 2 * pi * fG
  A <- sigma2 / 2 / sqrt(pi) * gamma(1 - 2 * d) / gamma(3 / 2 - 2 * d) * sin(omega) ^(1 - 4 * d)
  P1 <- Hypergeometric(1 - 2 * d, 1 - 2 * d, 3 / 2 - 2 * d, sin(omega / 2) ^ 2)
  P2 <- Hypergeometric(1 - 2 * d, 1 - 2 * d, 3 / 2 - 2 * d, cos(omega / 2) ^ 2)
  return(A * (P1 + P2))
}



Hypergeometric <- function(a, b, c, z) {
  ## Recursive implementation taken from Numerical Recipes in C (6.12)
  ## Press, Teukolsky, Vetterling and Flannery (1992)
  fac <- 1
  temp <- fac
  aa <- a
  bb <- b
  cc <- c
  for (n in 1:1000) {
    fac <- fac * (aa * bb) / cc
    fac <- fac * z / n
    series <- temp + fac
    if (series == temp)
      return(series)
    temp <- series
    aa <- aa + 1
    bb <- bb + 1
    cc <- cc + 1
  }
  stop("convergence failure in Hypergeometric")
}
