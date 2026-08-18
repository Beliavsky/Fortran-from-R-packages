#' Probability mass function for the number of classes from a Ewens distribution
#'
#' @details
#' The number of classes from a Ewens distribution with parameter \eqn{\theta} is given by the expression \deqn{Pr(K = k) = \lvert{} S^k_n \rvert{} \frac{\theta^k}{\theta (\theta + 1) ... (\theta + n - 1)}}, where \eqn{\lvert{}S^k_n \rvert{}} is the absolute value of a Stirling number of the first kind.
#' 
#' @param k An integer number of classes at which to evaluate the PMF
#' @param n A sample size not less than k
#' @param theta A non-negative parameter governing the expected sample diversity.
#' @returns The probability of observing k classes
#' @examples
#' x <- sample(LETTERS, 120, replace = TRUE)
#' dewens_k(1, 20, theta = 1) ## Pretty unlikely we just see one class
#'
#' @export
dewens_k <- function(k, n, theta) {
    abs(copula::Stirling1(n, k)) * (theta ^ k) /
        exp(sum(log(theta + seq_len(n) - 1)))
}
