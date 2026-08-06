#' Simulate Seasonal Persistent Processes Using the DWPT
#' 
#' A seasonal persistent process may be characterized by a spectral density
#' function with an asymptote occuring at a particular frequency in
#' \eqn{[0,\frac{1}{2})}{[0,1/2)}.  It's time domain representation was first
#' noted in passing by Hosking (1981).  Although an exact time-domain approach
#' to simulation is possible, this function utilizes the discrete wavelet
#' packet transform (DWPT).
#' 
#' Two subroutines are used, the first selects an adaptive orthonormal basis
#' for the true spectral density function (SDF) while the second computes the
#' bandpass variances associated with the chosen orthonormal basis and SDF.
#' Finally, when \eqn{M>N}{\code{M} > \code{N}} a uniform random variable is
#' generated in order to select a random piece of the simulated time series.
#' For more details see Whitcher (2001).
#' 
#' @param N Length of time series to be generated.
#' @param wf Character string for the wavelet filter.
#' @param delta Long-memory parameter for the seasonal persistent process.
#' @param fG Gegenbauer frequency.
#' @param M Actual length of simulated time series.
#' @param adaptive Logical; if \code{TRUE} the orthonormal basis used in the
#' DWPT is adapted to the ideal spectrum, otherwise the orthonormal basis is
#' performed to a maximum depth.
#' @param epsilon Threshold for adaptive basis selection.
#' @return Time series of length \code{N}.
#' @author B. Whitcher
#' @seealso \code{\link{hosking.sim}} for an exact time-domain method and
#' \code{\link{wave.filter}} for a list of available wavelet filters.
#' @references Hosking, J. R. M. (1981) Fractional Differencing,
#' \emph{Biometrika}, \bold{68}, No. 1, 165-176.
#' 
#' Whitcher, B. (2001) Simulating Gaussian Stationary Time Series with
#' Unbounded Spectra, \emph{Journal of Computational and Graphical Statistics},
#' \bold{10}, No. 1, 112-134.
#' @keywords ts
#' @examples
#' 
#' ## Generate monthly time series with annual oscillation
#' ## library(ts) is required in order to access acf()
#' x <- dwpt.sim(256, "mb16", .4, 1/12, M=4, epsilon=.001)
#' par(mfrow=c(2,1))
#' plot(x, type="l", xlab="Time")
#' acf(x, lag.max=128, ylim=c(-.6,1))
#' data(acvs.andel8)
#' lines(acvs.andel8$lag[1:128], acvs.andel8$acf[1:128], col=2)
#' 
#' @export dwpt.sim
dwpt.sim <- function(N, wf, delta, fG, M=2, adaptive=TRUE, epsilon=0.05) {
  M <- M*N
  J <- log(M, 2)
  jn <- rep(1:J, 2^(1:J))
  jl <- length(jn)

  if( adaptive ) {
    Basis <- find.adaptive.basis(wf, J, fG, epsilon) 
  } else {
    Basis <- numeric(jl)
    a <- min((1:jl)[jn == J])
    b <- max((1:jl)[jn == J])
    Basis[a:b] <- 1
  }

  Index <- (1:jl)[as.logical(Basis)]
  Length <- 2^jn

  variance <- bandpass.var.spp(delta, fG, J, Basis, Length)
  z <- vector("list", jl)
  class(z) <- "dwpt"
  attr(z, "wavelet") <- wf

  for(i in Index)
    z[[i]] <- rnorm(M/Length[i], sd=sqrt(Length[i]*variance[i]))

  x <- idwpt(z, Basis)
  xi <- trunc(runif(1, 1, M-N))
  return(x[xi:(xi+N-1)])
}



#' Determine an Orthonormal Basis for the Discrete Wavelet Packet Transform
#' 
#' Subroutine for use in simulating seasonal persistent processes using the
#' discrete wavelet packet transform.
#' 
#' The squared gain functions for a Daubechies (extremal phase or least
#' asymmetric) wavelet family are used in a filter cascade to compute the value
#' of the squared gain function for the wavelet packet filter at the
#' Gengenbauer frequency.  This is done for all nodes of the wavelet packet
#' table.
#' 
#' The idea behind this subroutine is to approximate the relationship between
#' the discrete wavelet transform and long-memory processes, where the squared
#' gain function is zero at frequency zero for all levels of the DWT.
#' 
#' @param wf Character string; name of the wavelet filter.
#' @param J Depth of the discrete wavelet packet transform.
#' @param fG Gegenbauer frequency.
#' @param eps Threshold for the squared gain function.
#' @return Boolean vector describing the orthonormal basis for the DWPT.
#' @author B. Whitcher
#' @seealso Used in \code{\link{dwpt.sim}}.
#' @keywords ts
#' @export find.adaptive.basis
find.adaptive.basis <- function(wf, J, fG, eps) {
  H <- function(f, L) {
    H <- 0
    for(l in 0:(L/2-1))
      H <- H + choose(L/2+l-1,l) * cos(pi*f)^(2*l)
    H <- 2 * sin(pi*f)^L * H
    return(H)
  }

  G <- function(f, L) {
    G <- 0
    for(l in 0:(L/2-1))
      G <- G + choose(L/2+l-1,l) * sin(pi*f)^(2*l)
    G <- 2 * cos(pi*f)^L * G
    return(G)
  }
  
  L <- wave.filter(wf)$length
  jn <- rep(1:J, 2^(1:J))
  jl <- length(jn)
  U <- numeric(jl)
  U[1] <- G(fG, L)
  U[2] <- H(fG, L)
  for(j in 2:J) {
    jj <- min((1:jl)[jn == j])
    jp <- (1:jl)[jn == j-1]
    for(n in 0:(2^j/2-1)) {
      if (n%%2 == 0) {
        U[jj + 2 * n + 1] <- U[jp[n+1]] * H(2^(j-1)*fG, L)
        U[jj + 2 * n] <- U[jp[n+1]] * G(2^(j-1)*fG, L)
      } else {
        U[jj + 2 * n] <- U[jp[n+1]] * H(2^(j-1)*fG, L)
        U[jj + 2 * n + 1] <- U[jp[n+1]] * G(2^(j-1)*fG, L)
      }
    }
  }
  return(ortho.basis(U < eps))
}

#' Bandpass Variance for Long-Memory Processes
#' 
#' Computes the band-pass variance for fractional difference (FD) or seasonal
#' persistent (SP) processes using numeric integration of their spectral
#' density function.
#' 
#' See references.
#' 
#' @usage bandpass.fdp(a, b, d)
#' @usage bandpass.spp(a, b, d, fG)
#' @usage bandpass.spp2(a, b, d1, f1, d2, f2)
#' @usage bandpass.var.spp(delta, fG, J, Basis, Length)
#' @aliases bandpass.fdp bandpass.spp bandpass.spp2 bandpass.var.spp
#' @param a Left-hand boundary for the definite integral.
#' @param b Right-hand boundary for the definite integral.
#' @param d,delta,d1,d2 Fractional difference parameter.
#' @param fG,f1,f2 Gegenbauer frequency.
#' @param J Depth of the wavelet transform.
#' @param Basis Logical vector representing the adaptive basis.
#' @param Length Number of elements in Basis.
#' @return Band-pass variance for the FD or SP process between \eqn{a} and
#' \eqn{b}.
#' @author B. Whitcher
#' @references McCoy, E. J., and A. T. Walden (1996) Wavelet analysis and
#' synthesis of stationary long-memory processes, \emph{Journal for
#' Computational and Graphical Statistics}, \bold{5}, No. 1, 26-56.
#' 
#' Whitcher, B. (2001) Simulating Gaussian stationary processes with unbounded
#' spectra, \emph{Journal for Computational and Graphical Statistics},
#' \bold{10}, No. 1, 112-134.
#' @keywords ts
bandpass.var.spp <- function(delta, fG, J, Basis, Length) {
  a <- unlist(sapply(2^(1:J)-1, seq, from=0, by=1)) / (2*Length)
  b <- unlist(sapply(2^(1:J), seq, from=1, by=1)) / (2*Length)
  bp.var <- rep(0, length(Basis))
  for(jn in (1:length(Basis))[as.logical(Basis)]) {
    if(fG < a[jn] | fG > b[jn])
      bp.var[jn] <- 2*integrate(spp.sdf, a[jn], b[jn], d=delta, fG=fG)$value
    else {
      result1 <- 2*integrate(spp.sdf, a[jn], fG, d=delta, fG=fG)$value
      result2 <- 2*integrate(spp.sdf, fG, b[jn], d=delta, fG=fG)$value
      bp.var[jn] <- result1 + result2
    }
  }
  return(bp.var)
}

