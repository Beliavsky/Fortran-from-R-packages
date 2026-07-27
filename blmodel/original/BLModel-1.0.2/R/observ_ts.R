#'Example of distribution of views -- Student t-distribution
#'
#'@description Function observ_ts computes density of Student t-distribution of views using the formula \cr
#'\eqn{f(x) = c_k*(1 +(x-q)^{T}*\Sigma^{-1}*(x-q)/df)^{(-(df+k)/2)}}, \cr
#'where \eqn{c_k} is a normalization constant (depends on the dimension of \eqn{x} and \eqn{q}) and \eqn{\Sigma} is the dispersion matrix. 
#'
#'@param x Data points matrix which collects in rows coordinates of points in which distribution density is  computed.
#'@param q Vector of investor's views.
#'@param covmat Covariance matrix of the distribution; dispersion matrix \eqn{\Sigma} is computed from \code{covmat}.
#'@param df Number of degrees of freedom of Students t-distribution.
#'
#'@return function returns a vector of observation distribution densities in data points x.
#'
#'@examples
#' k =3
#'observ_ts (x = matrix(c(rep(0.5,k),rep(0.2,k)),k,2), q = matrix(0,k,1), covmat = diag(k), 
#'          df=5)
#'
#'@references Kotz, S.,  Nadarajah, S., Multivariate t Distributions and Their Applications. Cambridge University Press,  2004.
#'
#'@export

observ_ts  <- function (x, q, covmat, df = 5)
{
  # for Student t-distribution of observations
  
  dfp = df
  k = ncol(covmat)
  Omega = (dfp-2)/dfp * covmat # dispersion matrix
  
  n = ncol(x)
  
  q = matrix(q,k,1) 
  aux1 =  t(q %*% matrix(1,1,n)- x)
  Omega_inv = solve(Omega)
  tpm = rowSums(((aux1) %*% Omega_inv) * (aux1))
  ck = gamma((dfp+k)/2)/(gamma(dfp/2)*sqrt(dfp^k*pi^k*det(Omega))) 
  odf = ck * (1 + tpm/dfp)^(-(dfp+k)/2)
  
  odf = (as.vector(odf, mode="numeric"))  
  odf = cbind((odf))
  dimnames ( odf) = NULL 
  return (odf)
}
