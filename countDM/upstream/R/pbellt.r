#' @export
#' @import stats
pbellt<-function (q, lambda, theta, lower.tail = TRUE, log.p = FALSE) {
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(theta <= 0))
    stop(paste("theta must be greater than 0"))
  if (any(q < 0))
    stop(paste("q must be 0 or greater than 0"))
  prob <- c()
  for ( i in 1:length(q) )  prob[i] <- sum( dbellt(0:q[i], lambda=lambda, theta=theta, log = FALSE) )
  if ( lower.tail == FALSE ) {
    if ( log.p == FALSE ) {
      prob <- 1 - prob
    } else  prob <- log(1 - prob)
  }  else if ( log.p == TRUE )  prob <- log(prob)
  prob
}
