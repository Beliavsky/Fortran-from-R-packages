#' Draw from a generalized Chinese Restaurant Process
#' 
#' @param n The sample size.
#' @param theta A non-negative parameter governing the expected sample diversity.
#' @param alpha A parameter between zero and one inclusive governing the expected sample diversity
#' @returns A vector of length n consisting of numeric class labels.
#' @examples
#' rewens(100, 1)
#' rewens(120, 0.5)
#' rewens(10, 0)
#' @export
gcrp <- function(n, alpha = 0, theta = 1) {
    .Call("gcrp_c", as.integer(n), as.numeric(alpha),
          as.numeric(theta), PACKAGE = "ewens")
}
