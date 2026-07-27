#'@title Portfolio optimization which finds an optimal portfolio with the smallest distance to a benchmark.
#'
#'@description PortfolioOptimProjection is a linear program for financial portfolio optimization. The function finds an optimal portfolio 
#'which has the smallest distance to a benchmark portfolio given by \code{bvec}.
#'Solution is by the algorithm due to Zhao and Li modified to account for the fact that the benchmark portfolio \code{bvec} has the dimension of portfolio weights 
#'and the solved linear program  has a much higher dimension since the solution vector to the LP problem consists of a set of primal variables: financial portfolio weights, 
#'auxiliary variables coming from the reduction of the mean-risk problem to a linear program and also a set of dual variables  depending  
#'on the number of constrains in the primal problem (see Palczewski).  
#'
#'@usage PortfolioOptimProjection (dat, portfolio_return,
#'risk=c("CVAR","DCVAR","LSAD","MAD"), alpha=0.95, bvec,
#'Aconstr=NULL, bconstr=NULL, LB=NULL, UB=NULL, maxiter=500, tol=1e-7)
#'
#'@param dat Time series of returns data; dat = cbind(rr, pk), where \eqn{rr} is an array (time series) of asset returns,  
#' for \eqn{n} returns and \eqn{k} assets it is an array with \eqn{\dim(rr) = (n, k)},    
#' \eqn{pk} is a vector of length \eqn{n} containing probabilities of returns.
#'@param portfolio_return Target portfolio return.
#'@param risk Risk measure chosen for optimization; one of "CVAR", "DCVAR", "LSAD", "MAD", where
#' "CVAR" -- denotes Conditional Value-at-Risk (CVaR),
#' "DCVAR" -- denotes deviation CVaR,
#' "LSAD" -- denotes Lower Semi Absolute Deviation,
#' "MAD" -- denotes Mean Absolute Deviation.
#'
#'@param alpha Value of alpha quantile used to compute portfolio VaR and CVaR; used also as quantile value for risk measures CVAR and DCVAR.
#'@param bvec Benchmark portfolio, a vector of length k; function \code{PortfolioOptimProjection} finds an optimal portfolio with the smallest distance to \code{bvec}.
#'@param Aconstr Matrix defining additional constraints,   \eqn{\dim (Aconstr) = (m,k)}, where
#' \eqn{k} -- number of assets, \eqn{m} -- number of constraints.
#'@param bconstr Vector defining additional constraints, length (\eqn{bconstr}) \eqn{ = m}.
#'@param LB Vector of length k, lower bounds of portfolio weights \eqn{\theta}; warning: condition LB = NULL is equivalent to LB = rep(0, k) (lower bound zero).
#'@param UB Vector of length k, upper bounds for portfolio weights \eqn{\theta}.
#'@param maxiter Maximal number of iterations. 
#'@param tol Accuracy of computations, stopping rule.
#'
#'
#'@return PortfolioOptimProjection returns a list with items:
#' \tabular{llll}{
#'\code{return_mean} \tab vector of asset returns  mean values. \cr
#'
#'\code{mu} \tab  realized portfolio return.\cr
#'
#'\code{theta} \tab  portfolio weights.\cr
#'
#'\code{CVaR} \tab  portfolio CVaR.\cr
#'
#'\code{VaR} \tab  portfolio VaR.\cr
#'
#'\code{MAD} \tab  portfolio MAD.\cr
#'
#'\code{risk} \tab  portfolio risk measured by the risk measure chosen for optimization.\cr
#'
#'\code{new_portfolio_return} \tab  modified target portfolio return; when the original target portfolio return \cr
#'
#'\code{ } \tab is to high for the problem, the optimization problem is solved for \cr
#'
#'\code{ } \tab new_portfolio_return as the target return. \cr

#'}
#'
#'@examples
#' 
#'library(mvtnorm)
#'k = 3 
#'num =100
#'dat <-  cbind(rmvnorm (n=num, mean = rep(0,k), sigma=diag(k)), matrix(1/num,num,1)) 
#'# a data sample with num rows and (k+1) columns for k assets;  
#'w_m <- rep(1/k,k) # benchmark portfolio, a vector of length k, 
#'port_ret = 0.05 # portfolio target return
#'alpha_optim = 0.95
#'
#'# minimal constraints set: \sum theta_i = 1
#'# has to be in two inequalities: 1 - \epsilon <= \sum theta_i <= 1 +\epsilon
#'a0 <- rep(1,k)
#'Aconstr <- rbind(a0,-a0)
#'bconstr <- c(1+1e-8, -1+1e-8)
#'
#'LB <- rep(0,k) 
#'UB <- rep(1,k)  
#'
#'res <- PortfolioOptimProjection(dat, port_ret, risk="MAD",  
#'alpha=alpha_optim, w_m, Aconstr, bconstr, LB, UB, maxiter=200, tol=1e-7)
#'
#'cat ( c("Projection optimal portfolio:\n\n"))
#'cat(c("weights \n"))
#'print(res$theta)
#'
#'
#'cat (c ("\n mean = ", res$mu, " risk  = ", res$risk,    "\n CVaR = ", res$CVaR, " VaR = ",
#' res$VaR, "\n MAD = ", res$MAD,  "\n\n"))
#'
#'
#'@references Palczewski, A., LP Algorithms for Portfolio Optimization: The PortfolioOptim Package, R Journal, 10(1) (2018), 308--327. DOI:10.32614/RJ-2018-028.
#' 
#'
#'Zhao, Y-B., Li, D., Locating the least 2-norm solution of linear programs via a path-following method, SIAM Journal on Optimization, 12 (2002), 893--912. DOI:10.1137/S1052623401386368.
#'
#'
#'@export


PortfolioOptimProjection <- function (dat, portfolio_return, risk = c("CVAR", "DCVAR", "LSAD", "MAD"), alpha = 0.95, bvec, 
							Aconstr =NULL, bconstr= NULL, LB = NULL, UB = NULL, maxiter = 500, tol = 1e-7)
{

  ep = tol
  rdigits <- -round( log10( tol ) )

  k = ncol(dat)-1
  n = nrow(dat)

  if (!is.null(Aconstr)){
    nVar <- k  # number of variables
    nCon <- length(bconstr)  # number of constraints

    if( !all.equal( dim( Aconstr ), c( nCon, nVar ) ) == TRUE ) {
      stop( paste( "Matrix A must have as many rows as constraints (=elements of vector b)",
                   " and as many columns as variables (=assets).\n" ) )
    }
  }

  if (!is.null(LB)){
    if (length(LB) != k){
      stop( paste( "Length of vector LB (=lower limits for portfolio weights) must be the same",
                   "as a number of assets.\n" ) )
    }
  }

  if (!is.null(UB)){
    if (length(UB) != k){
      stop( paste( "Length of vector UB (=upper limits for portfolio weights) must be the same",
                   "as a number of assets.\n" ) )
    }
  }

  if (length(bvec) != k){
    stop( paste( "Length of vector bvec \n",
                 "(selected optimal portfolio  is a projection of bvec on a set of optimal portfolios ) \n",
                 "must be the same as a number of assets.\n" ) )
  }

  pk = dat[,k+1]
  ra = as.matrix(dat[,1:k])

  ## Labels
  if( is.null(colnames(ra))) {
    ralab <- as.character(1:k)
  } else {
    ralab <- colnames(ra)
  }

  mu = matrix(0,1,k)
  for (i in 1:n ){
    mu = mu + ra[i, ] * pk[i]
  }
  colnames (mu) = ralab
  rownames (mu) = as.character("assets excess returns")

  lowb = LB[mu<0] * mu[mu<0]
  upb = UB[mu>0] * mu[mu>0]

  if (!is.null(UB) && ((sum(upb) + sum(lowb)) <= portfolio_return)) {
    cat ( c("Expected portfolio return is to high for the problem\n"))
    cat ( c("portfolio_return is set to highest accepted value \n"))
    portfolio_return =  (sum(upb) + sum(lowb)) - portfolio_return/20
    cat (c ("new portfolio_return = ", portfolio_return, "\n"))
  }

  rr = ra -  matrix(1,n,1)%*%mu
  dimnames (rr) = NULL

  risk = toupper(risk)

  risk <- match.arg(risk)

  cvarind = switch(risk,
                   CVAR = TRUE,
                   DCVAR = TRUE,
                   LSAD = FALSE,
                   MAD = FALSE)


  #############################################
  #
  #  solving optimization problem 
  #  where Lagrangian is given by
  #  L_\rho (x,y,z) = c^T x + y^T(z - A x + b) - \rho(\sum_1^n \log x_i + \sum_1^n \log z_i + 1/2 \rho^p(|B(x - x_eq)|^2 - |y|^2)
  #  this correspond to the solution of problem
  #  x*s = \rho e
  #  y *z = \rho e
  #  s +A^T y - c = \rho^p B^T(B(x - x_eq))
  #  z - A x + b = \rho^p y
  #
  #  with constraints
  #
  #  w^T*rr_k+\xi +eta_k  >= - \sum(LB * rr_k)
  #  eta_k >= 0
  #  \mu^T*w >= portfolio_return - \sum(LB * \mu)
  #  w <= UB - LB
  #
  #  x = (w[1:k], \xi, eta[1:n]) or x = (w[1:k], eta[1:n])
  #
  ###############################################

  if (cvarind) {
    clin = c(rep(0,k), 1, (1/(1-alpha))* t(pk))

  }
  else {
    clin = c(rep(0,k), t(pk))

  }

  if (cvarind) {
    initvec = c(bvec, rep(0,n+1))
    Bproj = .make_diag(c(rep(1,k),rep(0.001,n+1)))
  }
  else {
    initvec = c(bvec, rep(0,n))
    Bproj = .make_diag(c(rep(1,k),rep(0.001,n)))
  }

  if (is.null(Aconstr) ){
    Acon1 = NULL
  }
  else {
    rcon = nrow(Aconstr)
    if (cvarind) {
      a0 = matrix(0,rcon,n+1)}
    else {
      a0 = matrix(0,rcon,n)
    }
    Acon1 = cbind(Aconstr, a0)
  }

  if (cvarind) {
    a1 = c(mu,rep(0, n+1))}
  else {
    a1 = c(mu,rep(0, n))
  }

  if (cvarind) {
    a2 = matrix(1,n,1)}
  else {
    a2 = NULL
  }
  a3 = .make_diag(rep(1,n))
  amat1 = cbind(rr, a2, a3 )   # rr_k^T w + \xi + eta_k >= -sum(LB * r_k)

  a4 = .make_diag(rep(-1,k))
  if (cvarind) {
    a5 = matrix(0, k, n+1)}
  else {
    a5 = matrix(0, k, n)
  }

  amat2 = cbind(a4, a5)    #  w <= UB - LB
  if (is.null(UB) || is.null(LB)) amat2 = NULL

  Amat = rbind(Acon1, a1, amat1, amat2)
  dimnames (Amat) = NULL


  if (is.null(bconstr)) {
    bcon1 = NULL
  }
  else {
    if (is.null(LB)){
      bcon1 = bconstr
    }
    else {
      bcon1 = bconstr - LB %*% t(Aconstr)
    }
  }

  if (is.null(LB)){
    prcor = 0
    bcor1 = rep(0,n)
    bcor2 = NULL
    LB = 0
  }
  else {
    prcor = LB %*% t(mu)
    bcor1 = LB %*% t(rr)
    if ( is.null(UB)){
      bcor2 = NULL
    }
    else{
      bcor2 = LB - UB
    }
  }


  bmat = c(bcon1, portfolio_return - prcor, -bcor1, bcor2 )
  bmat = as.vector(bmat, mode="numeric")



  res = .ZI_projection(clin, Amat, bmat, initvec, Bproj, maxiter, tol, k)

  weights = res$xproj


  port_weights =  (weights[1:k]+LB )
  weight = round(c(port_weights), digits = rdigits)
  weight = matrix(weight, k,1)
  rownames (weight) = ralab
  colnames (weight) = as.character("portfolio weights")
  ret_mu = mu %*% port_weights

  loss = -ra %*% port_weights
  risk_st = .RISK_post (loss, pk, alpha)

  RM = switch(risk,
              CVAR = risk_st$cvar,
              DCVAR = risk_st$cvar+ret_mu,
              LSAD = risk_st$mad/2,
              MAD = risk_st$mad)
  if (ret_mu > (portfolio_return +ep)) {
    cat ( c("Target portfolio return is to small for the problem.\n"))
    cat ( c("Problem has been solved for target portfolio return = ", round(ret_mu, digits = rdigits), "\n"))
  }

  return (list (return_mean = (mu), mu = round(ret_mu, digits = rdigits), theta  = weight, CVaR = round(risk_st$cvar, digits = rdigits), 
		VaR = round(risk_st$var, digits = rdigits), MAD = round(risk_st$mad, digits = rdigits), risk = round(RM, digits = rdigits), 
		new_portfolio_return = portfolio_return ))


}
