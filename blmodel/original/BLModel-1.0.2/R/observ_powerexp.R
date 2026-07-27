#'Example of distribution of views -- power exponential distribution
#'
#'@description Function observ_powerexp computes density of power exponential distribution of views   using the formula\cr
#'\eqn{f(x) = c_k*\exp(- ((x-q)^{T}*\Sigma^{-1}*(x-q))^{\beta}/2)},\cr
#' where \eqn{c_k} is a normalization constant (depends on the dimension of \eqn{x} and \eqn{q}) and \eqn{\Sigma} is the dispersion matrix.  
#'
#'@param x Data points matrix which collects in rows coordinates of points in which distribution density is  computed.
#'@param q Vector of investor's views.
#'@param covmat Covariance matrix of the distribution; dispersion matrix \eqn{\Sigma} is computed from \code{covmat}. 
#'@param beta Shape parameter of the power exponential distribution.
#'
#'@return function returns a vector of distribution densities in data points x.
#'
#'@examples
#' k =3
#'observ_powerexp (x = matrix(c(rep(0.5,k),rep(0.2,k)),k,2), q = matrix(0,k,1),
#'                covmat = diag(k), beta = 0.6)
#'
#'@references Gomez, E., Gomez-Villegas, M., Marin, J.,  A multivariate generalization of the power exponential family of distributions. Commun. Statist. Theory Methods, 27 (1998), 589--600.
#'DOI: 10.1080/03610929808832115
#'@export

observ_powerexp  <- function (x, q,  covmat, beta = 0.6)
{
  # for power-exponential distribution of observations
  
  betal = beta
  k = ncol(covmat)
  Omega = k*gamma(k/(2*betal))/(2^(1/betal)*gamma((k+2)/(2*betal)))*covmat # dispersion matrix
  
  n = ncol(x)
  
  q = matrix(q,k,1) 
  aux1 =  t(q %*% matrix(1,1,n)- x)
  Omega_inv = solve(Omega)
  tpm = rowSums(((aux1) %*% Omega_inv) * (aux1))
  ck = k*gamma(k/2)/(pi^(k/2)* gamma(1 + k/(2*betal)) * 2^(1 + k/(2*betal)) *sqrt(det(Omega))) 
  odf = ck * exp(-tpm^betal/2 )
  
  odf = (as.vector(odf, mode="numeric"))  
  odf = cbind((odf))
  dimnames ( odf) = NULL 
  return (odf)
}
