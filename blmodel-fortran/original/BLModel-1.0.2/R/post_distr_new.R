#'Computes posterior distribution for discrete prior
#'
#'@param returns g
#'@param mu_m g
#'@param q g
#'@param P g
#'@param covmat g
#'@param tau g
#'@param FUN g
#'@param PARAM g
#'
#'

.post_distr_new <- function (dat, mu_m, q, P, covmat, tau, FUN, PARAM)
{
  # center the returns and shift to new mean mu_m
  k = ncol(dat)-1
  n = nrow(dat)
  pk = dat[,k+1]
  ra = as.matrix(dat[,1:k])


  mu = matrix(0,1,k)
  for (i in 1:n ){
    mu = mu + ra[i, ] * pk[i]
  }
  dimnames (mu) = NULL


  rr = ra -  matrix(1,n,1) %*% (mu - mu_m)
  colnames (rr) = colnames(ra)

  view_cov_type = switch(PARAM,
                         diag = TRUE,
                         full = FALSE,
                         stop( "wrong 'cov_matrix' type"))
  aux =  P %*% (covmat %*% t(P))
  if (view_cov_type)
  {
    Omega =  .make_diag (.diag_of(aux / tau))
  }
  else
  {
    Omega = aux/tau
  }

  data_points = P %*% t(rr)
  new_prob = pk * FUN(data_points, q, Omega)
  new_prob = cbind( new_prob/sum(new_prob))
  postdf = cbind(rr, new_prob)

  return (postdf)
}
