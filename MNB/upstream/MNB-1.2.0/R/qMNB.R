f.Ymas <- function(par,ymas,mu.mas){

  const <- exp(lgamma(ymas+par[1])-lfactorial(ymas)-lgamma(par[1]))
  q <- par[1]/(par[1]+mu.mas)
  func <- const*(q^(par[1]))*((1-q)^(ymas))
  return(func)
}

#' Randomized quantile residual
#'
#' @description randomized quantile residual is available to assess possible departures
#' from the multivariate negative binomial model for fitting correlated data with overdispersion.
#' @param par the maximum likelihood estimates.
#' @param formula The structure matrix of covariates of dimension n x p (in models that include an intercept x
#' should contain a column of ones).
#' @param dataSet data
#' @details The randomized quantile residual (Dunn and Smyth, 1996), which
#' follow a standard normal distribution is used to assess departures from the multivariate negative binomial model.
#' @return Randomized quantile Residuals
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @references
#' \itemize{
#' \item Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical
#' Statistics, 5, 236-244.
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
#' \item Fabio, L. C., Villegas, C., Mamun, A. S., and Carrasco, J. M. F. (2025).
#' Residual analysis for discrete correlated data in the multivariate approach.
#' Brazilian Journal of Biometrics, 43, e43728.
#' }
#' @examples
#'
#' data(seizures)
#' head(seizures)
#'
#' star <-list(phi=1, beta0=1, beta1=1, beta2=1, beta3=1)
#' mod <- fit.MNB(formula=Y ~ trt + period +
#' trt:period + offset(log(weeks)),star=star,dataSet=seizures,tab=FALSE)
#' par <- mod$par
#' names(par)<-c()
#'
#' res.q <- qMNB(par=par,formula=Y ~ trt + period + trt:period +
#' offset(log(weeks)),dataSet=seizures)
#'
#' plot(res.q,ylim=c(-3,4.5),ylab="Randomized quantile residual",
#' xlab="Index",pch=15,cex.lab = 1.5, cex = 0.6, bg = 5)
#' abline(h=c(-2,0,2),lty=3)
#' #identify(res.q)
#'
#'
#' data(alzheimer)
#' head(alzheimer)
#'
#' star <- list(phi=10,beta1=2, beta2=0.2)
#' mod <- fit.MNB(formula = Y ~ trat, star = star, dataSet = alzheimer,tab=FALSE)
#'
#' par<- mod$par
#' names(par) <- c()
#' re.q <- qMNB(par=par,formula = Y ~ trat, dataSet = alzheimer)
#' head(re.q)
#'
#' @export
#' @import stats


qMNB <- function(par,formula,dataSet){
  Y <-  stats::model.response(data=stats::model.frame(formula, dataSet))
  X <-  stats::model.matrix(formula, dataSet)
  off <- stats::model.extract(stats::model.frame(formula,dataSet),"offset")

  dataSet.ind <- split(dataSet, f= dataSet$ind)
  n <-  length(dataSet.ind)
  p <-  dim(X)[2]
  mi <- dim(dataSet.ind[[1]])[1]
  N <-  n*mi
  YI <- apply(matrix(Y,mi,n), 2,sum)

  phi <- par[1]
  beta <- par[2:(p+1)]

  if(methods::is(off,"NULL")){
    eta <- X %*% beta} else{eta <- X %*% beta + off}

  mu.i <- apply(matrix(exp(eta),mi,n), 2,sum)

  rq <- numeric(n)

  for(i in 1:n){
    Fy <- 0
    for(j in 0:YI[i]){
      temp1 <- f.Ymas(par=par,ymas=j,mu.mas=mu.i[i])
      Fy <- Fy + temp1
    }

    Fy.menos <- 0
    for(k in 0:(YI[i]-1)){
      temp2 <- f.Ymas(par=par,ymas=k,mu.mas=mu.i[i])
      Fy.menos <- Fy.menos + temp2
    }

    u <- stats::runif(1,Fy.menos,Fy)
    rq[i] <- stats::qnorm(u)

  }
  return(rq)
}
