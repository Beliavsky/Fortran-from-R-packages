#' @export
#' @import maxLik
#' @import stats
#' @import numbers
#' @importFrom miscTools stdEr
#'
mle_zip <- function(x, alpha, theta) {

 no <- sum(x == 0)
n <- length(x)
n1 <- n - no
x1 <- x[x > 0]

  loglik<-function(p){
    alpha <- p[1]
    theta <- p[2]
    loglik <- no * log(p[1] + (1 - p[1]) * exp(-p[2])) + n1 * log(1-p[1]) + sum(dpois(x1, p[2], log = TRUE))

  }
  res <- suppressWarnings(maxLik(loglik,start=c(alpha,theta)))
  aux <- cbind(coef(res), miscTools::stdEr(res))
  colnames(aux) = c("MLE", "SE")

  aux1 = cbind(AIC(res))
  colnames(aux1) = c("AIC")
  list(Estimates = aux, AIC= aux1)
}
