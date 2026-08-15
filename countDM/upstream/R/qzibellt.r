#' @export
#' @import stats

qzibellt <- function (p, lambda, theta, pi, lower.tail = TRUE, log.p = FALSE) {
  if (any(lambda <= 0))
    stop(paste("lambda must be greater than 0"))
  if (any(theta <= 0))
    stop(paste("theta must be greater than 0"))
  if ( any(pi < 0) | any(pi > 1) )  warning("pi must be in [0, 1]")
  if (any(p < 0) | any(p > 1))
    warning("pi must be in 0, 1")
  if (log.p)   p <- exp(p)
  if (!lower.tail)  p <- 1 - p
  p <- pmax(0, (p - pi) / (1 - pi))
  qbellt(p, lambda = lambda, theta = theta, lower.tail = TRUE, log.p = FALSE)

}
