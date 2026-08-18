#' Draw from the Ewens distribution
#'
#' Returns a vector with class membership
#'
#' @details
#' Although this command is described as sampling from the Ewens distribution, it is easier to think of it as a particular instantiation of the Chinese Restaurant Process, run for n "customers". The $j$th customer
#'  - sits at a new table with probability \eqn{\frac{\theta}{j - 1 + \theta}}, or
#'  - sits at an occupied table with probability \eqn{\frac{c}{j - 1 + \theta}}
#' where $c$ is the number of customers already at each table.
#' 
#' @param n The sample size.
#' @param theta A non-negative parameter governing the expected sample diversity.
#' @returns A vector of length n consisting of numeric class labels.
#' @examples
#' rewens(100, 1)
#' rewens(120, 0.5)
#' rewens(10, 0) ## equal to rep(1, 10)
#' @export
rewens <- function(n, theta = 1) {
    .Call("rewens_c", as.integer(n), as.numeric(theta), PACKAGE = "ewens")
}
