#' @title Ziggurat Pseudo Random Number Generators
#'
#' @description Pseudo random numbers drawn from one of three
#' distributions are offered, as proposed by Marsaglia and Tsang
#' (2000, <doi:10.18637/jss.v005.i08>) and improved by Leong et al
#' (2005, <doi:10.18637/jss.v012.i07>.
#'
#' @param n The size of the returned vector
#' @param s The chosen seed
#'
#' @return A vector of pseudo random-number generated draws from either the Normal,
#' Exponential or Uniform distribution, as selected.
#'
#' @references
#' <doi:10.18637/jss.v005.i08>
#' <doi:10.18637/jss.v012.i07>
#'
#' @examples
#' zsetseed(123)
#' zrnorm(3)
#' round(summary(zrunif(1e6)), 4)
zrnorm <- function(n) {
    .Call(cpp_zrnorm, n)
}

#' @rdname zrnorm
zrexp <- function(n) {
    .Call(cpp_zrexp, n)
}

#' @rdname zrnorm
zrunif <- function(n) {
    .Call(cpp_zrunif, n)
}

#' @rdname zrnorm
zsetseed <- function(s) {
    invisible(.Call(cpp_zsetseed, s))
}
