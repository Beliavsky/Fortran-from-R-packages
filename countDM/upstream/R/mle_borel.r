#' @export
#' @import maxLik
#' @import stats
#' @importFrom miscTools stdEr
mle_borel <- function(x, alpha) {

  if (any(alpha < 0) | any(alpha > 1))
    warning("alpha must be in 0, 1")
n <- length(x)

  loglik<-function(p){

    alpha <- p[1]
    loglik <- -p[1]*sum(x) + sum((x - 1) * log(p[1] * x)) - sum(lgamma(x +1))

  }
  res <- suppressWarnings(maxLik(loglik,start=c(alpha)))
  aux <- cbind(coef(res), miscTools::stdEr(res))
  colnames(aux) = c("MLE", "SE")

  aux1 = cbind(AIC(res))
  colnames(aux1) = c("AIC")
  list(Estimates = aux, AIC= aux1)
}

