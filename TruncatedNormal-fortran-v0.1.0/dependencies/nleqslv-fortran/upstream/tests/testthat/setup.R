expectedMessage1 <- "Function criterion near zero"
expectedMessage2 <- "x-values within tolerance 'xtol'"

common_test_f <- function(x) {
  y <- numeric(3)
  y[1] <- x[1] + x[2] - x[1]*x[2] - 2
  y[2] <- x[1] + x[3] - x[1]*x[3] - 3
  y[3] <- x[2] + x[3] - 4
  return(y)
}

common_test_jac <- function(x) {
  J <- matrix(0,nrow=3,ncol=3)
  J[,1] <- c(1-x[2], 1-x[3], 0)
  J[,2] <- c(1-x[1], 0, 1)
  J[,3] <- c(0, 1-x[1], 1)
  J
}

common_test_xsol <- c(-.5, 5/3 , 7/3)
