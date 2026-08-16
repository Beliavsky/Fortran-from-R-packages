#' The simulate_hp
#'
#' @description
#' Auxiliar function to generate a single observation for HYPERPO.
#'
#' @param mu a single value of the mu parameter.
#' @param sigma a single value of the sigma parameter.
#'
#' @return
#' a single value for the HYPERPO distribution.
#'
#' @export
#' @useDynLib DiscreteDists
#' @importFrom Rcpp sourceCpp
simulate_hp <- function(sigma, mu) {
  pochammer <- function(a, r) if (r == 0) 1 else prod(a:(a + r - 1))
  u <- stats::runif(1)
  y <- 0
  p <- 0
  value <- f11_cpp(gamma=sigma, lambda=mu)
  while (p < u) {
    p <- p + mu ^ y / (value * pochammer(sigma, y))
    y <- y + 1
  }
  y - 1
}
