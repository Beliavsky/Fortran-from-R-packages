#'for a given k-dimentional discrete distribution the function computes expected value vector and coviariance matrix of this distribution
#'
#'@param returns_coef g
#'@param returns g
#'
#'

.discrete_variance  <- function (returns_coef, returns )
{

  k = ncol(returns)-1
  n = nrow(returns)
  pk = returns[,k+1]
  ra = as.matrix(returns[,1:k])
  if( is.null(rownames(ra))) {
    clab <- as.character(1:k)
  } else {
    clab <- colnames(ra)
  }

  mu = matrix(0,1,k)
  for (i in 1:n ){
    mu = mu + ra[i, ] * pk[i]
  }

  rr = (ra -  matrix(1,n,1)%*%mu)
  dimnames (rr) = NULL

  for (i in 1:n ){
    rr[i, ] =  rr[i, ] * sqrt(pk[i])
  }
  mu = mu * returns_coef
  colnames (mu) = clab
  rownames (mu) = as.character("assets excess returns")


  cov = (t(rr) %*% rr) * returns_coef
  colnames (cov) = clab
  rownames (cov) = clab

  return (list (mu_disc = mu, variance = cov))

}
