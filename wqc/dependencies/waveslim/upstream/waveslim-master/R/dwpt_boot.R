#' Bootstrap Time Series Using the DWPT
#' 
#' An adaptive orthonormal basis is selected in order to perform the naive
#' bootstrap within nodes of the wavelet packet tree.  A bootstrap realization
#' of the time series is produce by applying the inverse DWPT.
#' 
#' A subroutines is used to select an adaptive orthonormal basis for the
#' piecewise-constant approximation to the underlying spectral density function
#' (SDF).  Once selected, sampling with replacement is performed within each
#' wavelet packet coefficient vector and the new collection of wavelet packet
#' coefficients are reconstructed into a bootstrap realization of the original
#' time series.
#' 
#' @param y Not necessarily dyadic length time series.
#' @param wf Name of the wavelet filter to use in the decomposition.  See
#' \code{\link{wave.filter}} for those wavelet filters available.
#' @param J Depth of the discrete wavelet packet transform.
#' @param p Level of significance for the white noise testing procedure.
#' @param frac Fraction of the time series that should be used in constructing
#' the likelihood function.
#' @return Time series of length $N$, where $N$ is the length of \code{y}.
#' @author B. Whitcher
#' @seealso \code{\link{dwpt.sim}}, \code{\link{spp.mle}}
#' @references Percival, D.B., S. Sardy, A. Davision (2000) Wavestrapping Time
#' Series: Adaptive Wavelet-Based Bootstrapping, in B.J. Fitzgerald, R.L.
#' Smith, A.T. Walden, P.C. Young (Eds.)  \emph{Nonlinear and Nonstationary
#' Signal Processing}, pp. 442-471.
#' 
#' Whitcher, B. (2001) Simulating Gaussian Stationary Time Series with
#' Unbounded Spectra, \emph{Journal of Computational and Graphical Statistics},
#' \bold{10}, No. 1, 112-134.
#' 
#' Whitcher, B. (2004) Wavelet-Based Estimation for Seasonal Long-Memory
#' Processes, \emph{Technometrics}, \bold{46}, No. 2, 225-238.
#' @keywords ts
#' @export dwpt.boot
dwpt.boot <- function(y, wf, J=log(length(y),2)-1, p=1e-04, frac=1) {

  N <- length(y)
  if(N/2^J != trunc(N/2^J)) 
    stop("Sample size is not divisible by 2^J")
  
  ## Perform discrete wavelet packet transform (DWPT) on Y
  y.dwpt <- dwpt(y, wf, n.levels=J)
  n <- length(y)
  if(frac < 1) {
    for(i in 1:length(y.dwpt)) {
      vec <- y.dwpt[[i]]
      ni <- length(vec)
      j <- rep(1:J, 2^(1:J))[i]
      vec[trunc(frac * n/2^j):ni] <- NA
      y.dwpt[[i]] <- vec
    }
  }
  y.basis <- as.logical(ortho.basis(portmanteau.test(y.dwpt, p, type="other")))

  ## Taken from my 2D bootstrapping methodology
  resample.dwpt <- y.dwpt
  for(i in 1:length(y.basis)) {
    m <- length(y.dwpt[[i]])
    if(y.basis[i])
      resample.dwpt[[i]] <- sample(y.dwpt[[i]], replace=TRUE)
    else
      resample.dwpt[[i]] <- rep(NA, m)
  }
  idwpt(resample.dwpt, y.basis)
}

