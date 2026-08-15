#' @export
#' @import stats
dzibellt<-function (x, lambda, theta, pi, log = FALSE) {
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(theta <= 0))
    stop(paste("theta must be greater than 0"))
  if (any(pi < 0) | any(pi > 1))
    warning("pi must be in 0, 1")
  rval <- log(1 - pi) + dbellt(x, lambda, theta, log = TRUE)
  if (any(x0 <- (x == 0L)))
    rval[x0] <- log(exp(rval) + pi)[x0]
  if (log)
    rval
  else exp(rval)
}
