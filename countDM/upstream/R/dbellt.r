#' @export
#' @import stats

dbellt<-function(x, lambda, theta, log = FALSE) {
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(theta <= 0))
    stop(paste("theta must be greater than 0"))

    TPx <- c()
    for (i in 1:length(x)) {
        TPx[i] <- TP(x[i], theta=theta)
    }

	lf <- x * log(lambda)+ theta * ( -exp(lambda) + 1 ) + log(TPx) - lgamma(x + 1)
    if (log == TRUE) {
        return(lf)
    }
    else {
        return(exp(lf))
    }
}

