#' Calculating Thomson's Spectral Multitapers by Inverse Iteration
#' 
#' This is now a wrapper to the function multitaper::dpss().
#' 
#' @param n length of data taper(s)
#' @param k number of data tapers; 1, 2, 3, ... (do not use 0!)
#' @param nw product of length and half-bandwidth parameter (w)
#' @return \item{v}{matrix of data tapers (cols = tapers)}
#' \item{eigen}{eigenvalue associated with each data taper, discarded} 
#' @author B. Whitcher
#' @seealso \code{\link{sine.taper}}.
#' @references Percival, D. B. and A. T. Walden (1993) \emph{Spectral Estimation for
#' Physical Applications: Multitaper and Conventional Univariate Techniques},
#' Cambridge University Press.
#' @keywords ts
#' @export dpss.taper
dpss.taper <- function(n, k, nw = 4) {
  out <- multitaper::dpss(n, k, nw)
  return(out$v)
}

#' Computing Sinusoidal Data Tapers
#' 
#' Computes sinusoidal data tapers directly from equations.
#' 
#' See reference.
#' 
#' @param n length of data taper(s)
#' @param k number of data tapers
#' @return A vector or matrix of data tapers (cols = tapers).
#' @author B. Whitcher
#' @references Riedel, K. S. and A. Sidorenko (1995) Minimum bias multiple
#' taper spectral estimation, \emph{IEEE Transactions on Signal Processing},
#' \bold{43}, 188-195.
#' @keywords ts
#' @export sine.taper
sine.taper <- function(n, k) {
  tapers <- NULL
  for(i in 1:k)
    tapers <- cbind(tapers, sqrt(2/(n+1)) * sin((pi*i*1:n)/(n+1)))
  return(tapers)
}
