#' Calculate expected number of classes in a sample of size n given theta
#'
#' @description
#' The expected number of classes from the Ewens distribution is given by \eqn{\theta \sum_{j=1}^{n} \frac{1}{\theta + j - 1}}. This is often more convenient than integrating across the PMF given by `dewens_k`
#' 
#' @param n The sample size
#' @param theta The non-negative parameter governing expected sample diversity
#' 
ewens_k_exact <- function(n, theta) {
    theta * sum(1 / (theta + seq(0, n - 1)))
}
