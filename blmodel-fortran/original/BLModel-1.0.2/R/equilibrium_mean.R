#'Solves the inverse optimization to mean-risk standard optimization  problem  to find equilibrium returns.
#' The function is invoked by BL_post_distr and arguments are supplemented by BL_post_distr.
#'
#'@description The function computes the vector of equilibrium returns implied by a market portfolio.
#'The vector of means for the mean-risk optimization problem is found by inverse optimization. \cr
#'The optimization problem is:\cr
#'\eqn{\min F(w_m^{T} r)}\cr
#'subject to\cr
#'\eqn{w_m^{T} E(r) \ge RM},\cr
#'where \cr
#'\eqn{F} is a risk measure -- one from the list  c("CVAR", "DCVAR", "LSAD", "MAD"),\cr
#'\eqn{r}  is a time series of market returns,\cr
#'\eqn{w_m}  is  market portfolio,\cr
#'\eqn{RM}  is  market expected return.
#'
#'@param dat Time series of returns data; dat = cbind(rr, pk), where \eqn{rr} is an array (time series) of market asset returns,
#' for \eqn{n} returns and \eqn{k} assets it is an array with \eqn{\dim(rr) = (n, k)},
#' \eqn{pk} is a vector of length \eqn{n} containing probabilities of returns.
#'@param w_m Market portfolio.
#'@param RM Market_expected_return.
#'@param risk A risk measure, one from the list  c("CVAR", "DCVAR", "LSAD", "MAD").
#'@param alpha Value of alpha quantile in the definition of risk measures CVAR and DCVAR. Can be any number when risk measure is parameter free.
#'@return
#' \tabular{llll}{
#'\code{market_returns}  \tab a vector of market returns obtain by inverse optimization; this is vector \eqn{E(r)}\cr
#'
#'\code{ }  \tab  from the description of this function.
#'}
#'
#'@examples
#' 
#'# In normal usage all data are supplemented by function BL_post_distr.
#'library(mvtnorm)
#'k = 3 
#'num =100
#'dat <-  cbind(rmvnorm (n=num, mean = rep(0,k), sigma=diag(k)), matrix(1/num,num,1)) 
#'# a data sample with num rows and (k+1) columns for k assets;
#'w_m <- rep(1/k,k) # market portfolio.
#'RM = 0.05 # market expected return.
#'equilibrium_mean (dat, w_m, RM, risk = "CVAR", alpha = 0.95) 
#'
#'@references Palczewski, J., Palczewski, A., Black-Litterman Model for Continuous Distributions (2016). Available at SSRN: https://ssrn.com/abstract=2744621.
#'@export


equilibrium_mean  <- function (dat, w_m, RM, risk = c("CVAR", "DCVAR", "LSAD", "MAD"), alpha=0.95)
{

  k = ncol(dat)-1
  n = nrow(dat)
  x_m = w_m
  RM = RM
  if (length(x_m) != k){
    stop( paste( "Length of a vector of assets weights  must be the same",
                 "as a number of assets for which returns data are provided.\n" ) )
  }
  if (RM <= 0){
    stop( paste("Market portfolio return must be positive.\n" ) )
  }

  risk = toupper(risk)
  risk <- match.arg(risk)
  cvarind = switch(risk,
                   CVAR = TRUE,
                   DCVAR = TRUE,
                   LSAD = FALSE,
                   MAD = FALSE)

  # center the returns
  # changing sign we pass from returns to losses
  ra = -as.matrix(dat[,1:k])
  pk = dat[,k+1]

  mu = matrix(0,1,k)
  for (i in 1:n ){
    mu = mu + ra[i, ] * pk[i]
  }
  dimnames (mu) = NULL
  rr = ra -  matrix(1,n,1)%*%mu
  dimnames (rr) = NULL

  ## Labels
  if( is.null(colnames(ra))) {
    ralab <- as.character(1:k)
  } else {
    ralab <- colnames(ra)
  }

  returns_m = rr %*% x_m
  o = order (returns_m)
  sorted_returns = as.matrix(returns_m [o])
  weight = as.matrix( pk[o])
  if (cvarind) {
    index = sum(cumsum(weight) < alpha) +1
    mu_m =  t( t(rr [o [(index + 1):length(o)],])%*% (weight[(index + 1):length(o)])) + (sum(weight[1:index]) - alpha) * rr [o[index],]
  }
  else {
    index = sum(sorted_returns  <= 0 )
    mu_m =  t( t(rr [o [1:(index)],])%*% (weight[1:(index)]))
  }

  u0 = (mu_m %*% x_m)/RM

  mu_m = mu_m / u0[1,1]
  mu_m = t(mu_m)
  rownames (mu_m) = ralab
  colnames (mu_m) = as.character("assets excess returns")
  return (list (market_returns = mu_m))
}
