#'Example of distribution of views -- normal distribution
#'
#'@description Function observ_normal computes density of normal distribution of views using the formula \cr
#'\eqn{f(x) = c_k*\exp(-((x-q)^{T}*covmat^{-1}*(x-q))/2)},\cr
#'where \eqn{c_k} is a normalization constant (depends on the dimension of \eqn{x} and \eqn{q}).
#'
#'@param x Data points matrix which collects in rows coordinates of points in which distribution density is  computed.
#'@param q Vector of investor's views.
#'@param covmat Covariance matrix of the distribution.
#'
#'@return function returns a vector of distribution densities in data points x.
#'
#'@examples
#' k =3
#' observ_normal (x = matrix(c(rep(0.5,k),rep(0.2,k)),k,2), q = matrix(0,k,1), 
#'                covmat = diag(k)) 
#'
#'@references Palczewski, J., Palczewski, A., Black-Litterman Model for Continuous Distributions (2016). Available at SSRN: https://ssrn.com/abstract=2744621.
#'
#'@export

observ_normal  <- function (x, q, covmat )
{
  
  # for normal distributions of observations
  k = ncol(covmat)
  
  Omega = covmat # dispersion matrix for normal distribution
  n = ncol(x)
  q = matrix(q,k,1) 
  aux1 =  t(q %*% matrix(1,1,n)- x)
     
  Omega_inv = solve(Omega)
  ck = 1/sqrt((2*pi)^k*det(Omega))
  tpm = rowSums(((aux1) %*% Omega_inv) * (aux1))
  odf = ck * exp(-tpm/2)
  
  odf = (as.vector(odf, mode="numeric"))  
  odf = cbind((odf)) 
  dimnames ( odf) = NULL 
  return (odf)
  
}
