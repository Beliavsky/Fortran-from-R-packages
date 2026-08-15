#' @export
#' @import maxLik
#' @import stats
#' @import numbers
#' @importFrom miscTools stdEr
#'
mle_zoibell <- function(x, alpha, beta, theta) {
n<- length(x)
n
S0 <- sum( x == 0)
S1 <- sum( x == 1 )
x1 <- x[x > 1]
x1
S <- sum(x1)
com <- ( n - S0 - S1 )

  loglik<-function(p){

    alpha <- p[1]
	beta <- p[2]
	theta <- p[3]
    loglik <-  S0*log(p[1])+S1*log(p[2])+(com)*log(1-p[1]-p[2])+S*log(p[3])+com*(1-exp(p[3]))-com*log(1-exp(1-exp(p[3]))-p[3]*exp(1-exp(p[3])))-sum( lgamma(x1 + 1) )+sum(log(sapply(x1,bell)))

  }
  res <- suppressWarnings(maxLik(loglik,start=c(alpha, beta, theta)))
  aux <- cbind(coef(res), miscTools::stdEr(res))
  colnames(aux) = c("MLE", "SE")

  aux1 = cbind(AIC(res))
  colnames(aux1) = c("AIC")
  list(Estimates = aux, AIC= aux1)
}
