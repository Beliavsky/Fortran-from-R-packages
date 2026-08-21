rbp <- function(n, lambda) {
  x <- cbind( rpois(n, lambda[1]), rpois(n, lambda[2]) )
  x + rpois(n, lambda[3])
}
