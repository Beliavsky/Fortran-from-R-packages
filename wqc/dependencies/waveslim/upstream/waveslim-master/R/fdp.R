#' Wavelet-based Maximum Likelihood Estimation for a Fractional Difference
#' Process
#' 
#' Parameter estimation for a fractional difference (long-memory, self-similar)
#' process is performed via maximum likelihood on the wavelet coefficients.
#' 
#' The variance-covariance matrix of the original time series is approximated
#' by its wavelet-based equivalent.  A Whittle-type likelihood is then
#' constructed where the sums of squared wavelet coefficients are compared to
#' bandpass filtered version of the true spectrum.  Minimization occurs only
#' for the fractional difference parameter \eqn{d}, while variance is estimated
#' afterwards.
#' 
#' @param y Dyadic length time series.
#' @param wf Name of the wavelet filter to use in the decomposition.  See
#' \code{\link{wave.filter}} for those wavelet filters available.
#' @param J Depth of the discrete wavelet transform.
#' @return List containing the maximum likelihood estimates (MLEs) of \eqn{d}
#' and \eqn{\sigma^2}, along with the value of the likelihood for those
#' estimates.
#' @author B. Whitcher
#' @references M. J. Jensen (2000) An alternative maximum likelihood estimator
#' of long-memory processes using compactly supported wavelets, \emph{Journal
#' of Economic Dynamics and Control}, \bold{24}, No. 3, 361-387.
#' 
#' McCoy, E. J., and A. T. Walden (1996) Wavelet analysis and synthesis of
#' stationary long-memory processes, \emph{Journal for Computational and
#' Graphical Statistics}, \bold{5}, No. 1, 26-56.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @examples
#' 
#' ## Figure 5.5 in Gencay, Selcuk and Whitcher (2001)
#' fdp.sdf <- function(freq, d, sigma2=1)
#'   sigma2 / ((2*sin(pi * freq))^2)^d
#' dB <- function(x) 10 * log10(x)
#' per <- function(z) {
#'   n <- length(z)
#'   (Mod(fft(z))**2/(2*pi*n))[1:(n %/% 2 + 1)]
#' }
#' data(ibm)     
#' ibm.returns <- diff(log(ibm))
#' ibm.volatility <- abs(ibm.returns)
#' ibm.vol.mle <- fdp.mle(ibm.volatility, "d4", 4)
#' freq <- 0:184/368
#' ibm.vol.per <- 2 * pi * per(ibm.volatility)
#' ibm.vol.resid <- ibm.vol.per/ fdp.sdf(freq, ibm.vol.mle$parameters[1])
#' par(mfrow=c(1,1), las=0, pty="m")
#' plot(freq, dB(ibm.vol.per), type="l", xlab="Frequency", ylab="Spectrum")
#' lines(freq, dB(fdp.sdf(freq, ibm.vol.mle$parameters[1],
#'                        ibm.vol.mle$parameters[2]/2)), col=2)
#' 
#' @export fdp.mle
fdp.mle <- function(y, wf, J=log(length(y),2))
{
  fdpML <- function(d, y) {
    
    y.dwt <- y[[1]]
    n <- y[[2]]
    J <- y[[3]]

    ## Establish the limits of integration for the band-pass variances
    a <- c(1/2^c(1:J+1), 0)
    b <- 1/2^c(0:J+1)

    ## Define some useful parameters for computing the likelihood
    length.j <- n / c(2^(1:J), 2^J)
    scale.j <- c(2^(1:J+1), 2^(J+1))

    ## Initialize various parameters for computing the approximate ML
    bp.var <- numeric(J+1)

    ## Compute the band-pass variances according to d
    omega.diag <- NULL
    for(j in 1:(J+1)) {
      bp.var[j] <- integrate(fdp.sdf, a[j], b[j], d=d)$value
      omega.diag <- c(omega.diag, scale.j[j] * rep(bp.var[j], length.j[j]))
    }
    
    ## Compute approximate maximum likelihood 
    n * log(sum(y.dwt^2 / omega.diag) / n) +
      sum(length.j * log(scale.j * bp.var)) 
  }

  n <- length(y)
  y.dwt <- as.vector(unlist(dwt(y, wf, n.levels=J)))

  ## Compute MLE of d (limited to stationary region)
  result <- optimize(fdpML, interval=c(-0.5,0.5), maximum=FALSE,
                     y=list(y.dwt, n, J))

  ## Compute MLE of sigma_epsilon^2
  a <- c(1/2^c(1:J+1), 0)
  b <- 1/2^c(0:J+1)
  length.j <- n / c(2^(1:J), 2^J)
  scale.j <- c(2^(1:J+1), 2^(J+1))
  bp.var <- numeric(J+1)
  omega.diag <- NULL
  for(j in 1:(J+1)) {
    bp.var[j] <- integrate(fdp.sdf, a[j], b[j], d=result$minimum)$value
    omega.diag <- c(omega.diag, scale.j[j] * rep(bp.var[j], length.j[j]))
  }
  sigma2 <- sum(y.dwt^2 / omega.diag) / n

  list(parameters=c(result$minimum, sigma2), objective=result$objective)
}

