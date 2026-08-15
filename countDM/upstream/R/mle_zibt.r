#' @export
#' @import maxLik
#' @import stats
#' @import numbers
#' @importFrom miscTools stdEr

mle_zibellt <- function(x, lambda, theta, pi) {
n<- length(x)
n
S0 <- sum( x == 0)
x1 <- x[x > 0]
x1
S <- sum(x1)

com <- (n - S0)

loglik<-function(p){
     pi <- p[1]
    lambda <- p[2]
 	  theta <- p[3]

loglik<-S0*log(p[1]+(1-p[1])*exp(p[3]*(1-exp(p[2]))))+com*log(1-p[1])+com*p[3]*(1-exp(p[2]))+S*log(p[2])-sum( lgamma(x1 + 1) )+sum(log(sapply(x1,theta=p[3],TP)))

}

res <- suppressWarnings(maxLik(loglik,start=c(lambda, theta, pi)))
  aux <- cbind(coef(res), miscTools::stdEr(res))
  colnames(aux) = c("MLE", "SE")

  aux1 = cbind(AIC(res))
  colnames(aux1) = c("AIC")
  list(Estimates = aux, AIC= aux1)
}
