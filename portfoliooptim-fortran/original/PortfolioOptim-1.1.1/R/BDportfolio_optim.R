#'@title Portfolio Optimization by Benders decomposition 
#' 
#'@description BDportfolio_optim is a linear program  for financial portfolio optimization. 
#' Portfolio risk is measured by one of the risk measures from the list  c("CVAR", "DCVAR", "LSAD", "MAD"). 
#' Benders decomposition method is explored to enable optimization for very large returns samples (\eqn{\sim 10^6}). 
#' 
#' The optimization problem is:\cr  
#' \eqn{\min F({\theta^{T}} r)} \cr  
#' over \cr 
#' \eqn{\theta^{T} E(r)} \eqn{\ge} \eqn{portfolio\_return}, \cr 
#' \eqn{LB} \eqn{\le \theta \le} \eqn{UB}, \cr 
#' \eqn{Aconstr} \eqn{\theta \le} \eqn{bconstr}, \cr 
#' where \cr 
#' \eqn{F} is a measure of risk; \cr 
#' \eqn{r} is a time series of returns of assets; \cr 
#' \eqn{\theta} is a vector of portfolio weights. \cr 
#' 
#'@usage BDportfolio_optim(dat, portfolio_return,  
#'risk=c("CVAR", "DCVAR","LSAD","MAD"), alpha=0.95,  
#'Aconstr=NULL, bconstr=NULL, LB=NULL, UB=NULL, maxiter=500,tol=1e-8) 
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
#'@param alpha Value of alpha quantile used to compute portfolio VaR and CVaR; used also as quantile value for risk measures CVAR and DCVAR. 
#'@param Aconstr Matrix defining additional constraints,   \eqn{\dim (Aconstr) = (m,k)}, where 
#' \eqn{k} -- number of assets, \eqn{m} -- number of constraints. 
#'@param bconstr Vector defining additional constraints, length (\eqn{bconstr}) \eqn{ = m}. 
#'@param LB Vector of length k, lower bounds of portfolio weights \eqn{\theta}; warning: condition LB = NULL is equivalent to LB = rep(0, k) (lower bound zero). 
#'@param UB Vector of length k, upper bounds for portfolio weights \eqn{\theta}. 
#'@param maxiter Maximal number of iterations.  
#'@param tol Accuracy of computations, stopping rule. 
#' 
#' 
#'@return BDportfolio_optim returns a list with items: 
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
#'library (Rsymphony)  
#'library(Rglpk) 
#'library(mvtnorm)
#'k = 3 
#'num =100
#'dat <-  cbind(rmvnorm (n=num, mean = rep(0,k), sigma=diag(k)), matrix(1/num,num,1)) 
#'# a data sample with num rows and (k+1) columns for k assets; 
#'port_ret = 0.05 # target portfolio return 
#'alpha_optim = 0.95 
#' 
#'# minimal constraints set: \eqn{\sum \theta_{i} = 1} 
#'# has to be in two inequalities: \eqn{1 - \epsilon <= \sum \theta_{i} <= 1 + \epsilon} 
#'a0 <- rep(1,k) 
#'Aconstr <- rbind(a0,-a0) 
#'bconstr <- c(1+1e-8, -1+1e-8) 
#' 
#'LB <- rep(0,k) 
#'UB <- rep(1,k) 
#' 
#'res <- BDportfolio_optim(dat, port_ret, "CVAR", alpha_optim, 
#' Aconstr, bconstr, LB, UB, maxiter=200, tol=1e-8) 
#' 
#'cat ( c("Benders decomposition portfolio:\n\n")) 
#'cat(c("weights \n")) 
#'print(res$theta) 
#' 
#'cat(c("\n mean = ", res$mu, " risk  = ", res$risk, 
#'"\n CVaR = ", res$CVaR, " VaR = ", res$VaR, "\n MAD = ", res$MAD, "\n\n")) 
#'  
#' 
#'@references Benders, J.F.,   Partitioning procedures for solving mixed-variables programming problems. Number. Math., 4 (1962),  238--252, reprinted in 
#' Computational Management Science 2 (2005), 3--19. DOI: 10.1007/s10287-004-0020-y. 
#' 
#'Konno, H., Piecewise linear risk function and portfolio optimization, Journal of the Operations Research Society of Japan, 33 (1990), 139--156. 
#' 
#'Konno, H., Yamazaki, H.,  Mean-absolute deviation portfolio optimization model and its application to Tokyo stock market. Management Science, 37 (1991), 519--531. 
#' 
#'Konno, H., Waki, H.,   Yuuki, A., Portfolio optimization under lower partial risk measures, Asia-Pacific Financial Markets, 9 (2002), 127--140. DOI: 10.1023/A:1022238119491. 
#' 
#'Kunzi-Bay, A.,  Mayer, J., Computational aspects of minimizing conditional value at risk. Computational Management Science, 3 (2006), 3--27. DOI: 10.1007/s10287-005-0042-0. 
#' 
#'Rockafellar, R.T., Uryasev, S.,  Optimization of conditional value-at-risk. Journal of Risk, 2 (2000), 21--41. DOI: 10.21314/JOR.2000.038.  
#' 
#'Rockafellar, R. T., Uryasev, S.,  Zabarankin, M.,  Generalized deviations in risk analysis. Finance and Stochastics, 10 (2006), 51--74. DOI: 10.1007/s00780-005-0165-8. 
#' 
#' 
#'@export 


BDportfolio_optim <- function (dat, portfolio_return, risk = c("CVAR", "DCVAR", "LSAD", "MAD"),  alpha = 0.95,  
					Aconstr =NULL, bconstr= NULL, LB = NULL, UB = NULL, maxiter = 500, tol = 1e-8 ) 
{ 

 ep = tol # tolerance for optimal solution 
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



 clin = c(rep(0,k), 1) 
 if (cvarind)  clin = c(clin, 1/(1-alpha)) 


 if (is.null(Aconstr) ){ 
   Acon1 = NULL 
 } 
 else { 
   rcon = nrow(Aconstr) 
   if (cvarind) { 
     a0 = matrix(0,rcon,2)} 
   else { 
     a0 = matrix(0,rcon,1) 
   } 
   Acon1 = cbind(Aconstr, a0) 
 } 

 if (cvarind) { 
   a1 = c(mu,0,0)} 
 else { 
   a1 = c(mu,0) 
 } 

 a4 = .make_diag(rep(1,k)) 
 if (cvarind) { 
   a5 = matrix(0, k, 2)} 
 else { 
   a5 = matrix(0, k, 1) 
 } 
 amat2 = cbind(a4, a5) 
 if (is.null(UB) || is.null(LB)) amat2 = NULL 

 Amat = rbind(Acon1, -a1, amat2) 

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
   bcor1 = 0 
   bcor2 = NULL 
   LB = 0 
 } 
 else { 
   prcor = LB %*% t(mu) 
   if ( is.null(UB)){ 
     bcor2 = NULL 
   } 
   else{ 
     bcor2 = UB - LB 
   } 
 } 


 bmat = c(bcon1,-portfolio_return + prcor, bcor2) 


 K_set <- rep(1,n) 

 loop = 0 

 repeat 
 { 
   loop = loop+1 
   eta = t(pk * K_set) %*% rr 
   p_i = sum( pk * K_set ) 
   if (cvarind) { 
     a3 = c(eta,p_i,1)} 
   else { 
     a3 = c(eta,1) 
   } 
   Amat = rbind(Amat, -a3) 
   dimnames (Amat) = NULL 

   if (length(LB) > 1 ) bcor1 = LB %*% t(eta) 


   bmat = c(bmat, bcor1) 
   bmat = as.vector(bmat, mode="numeric") 

   res = Rsymphony::Rsymphony_solve_LP(clin,  Amat, dir = rep( "<=", length(bmat)), bmat,  max = FALSE, node_limit = 3000)

   sol = res$solution 

   if (cvarind) { 
     K_tem = rr %*% (sol[1:k]+ LB) + sol[k+1]*matrix(1,n,1)} 
   else { 
     K_tem = rr %*% (sol[1:k]+ LB) 
   } 
   K_set <- K_tem < 0 

   w_plus = -t(pk * K_set) %*% K_tem 
   if (cvarind) { 
     w_star = -(eta %*% (sol[1:k]+ LB) + p_i *sol[k+1])} 
   else { 
     w_star = -(eta %*% (sol[1:k]+ LB)) 
   } 
   if (cvarind) { 
     F_bar = sol[k+1] + w_plus/ (1- alpha) 
     F_down = sol[k+1] + w_star/ (1- alpha) 
   } 
   else { 
     F_bar =   w_plus 
     F_down =   w_star 
   } 


   if (loop == 1) F_min = F_bar 

   F_min = min(F_bar, F_min) 

   if (F_min - F_down <= ep)   break 
   if (loop > maxiter){ 
			cat(c("Maximal number of iteration has been exceeded. \n")) 
			break 
	} 

 } 
 port_weights =  (sol[1:k]+LB ) 
 weight = round(c(port_weights), digits = rdigits) 
 weight = matrix(weight, k,1) 
 rownames (weight) = ralab 
 colnames (weight) = as.character("portfolio weights") 
 ret_mu = mu %*% port_weights 

 RM = F_down 

 RM = switch(risk, 
             CVAR = RM - ret_mu, 
             DCVAR = RM, 
             LSAD = RM, 
             MAD = 2* RM) 

 if (ret_mu > (portfolio_return +ep)) { 
   cat ( c("Target portfolio return is to small for the problem.\n")) 
   cat ( c("Problem has been solved for target portfolio return = ", round(ret_mu, digits = rdigits), "\n")) 
 } 

 loss = -ra %*% port_weights 

 risk_st = .RISK_post (loss, pk, alpha) 

 return (list (return_mean = mu, mu = round(ret_mu, digits = rdigits), theta  = weight, CVaR = round(risk_st$cvar, digits = rdigits),  
		VaR = round(risk_st$var, digits = rdigits), MAD = round(risk_st$mad, digits = rdigits), risk = round(RM, digits = rdigits),  
		new_portfolio_return = portfolio_return )) 


} 
