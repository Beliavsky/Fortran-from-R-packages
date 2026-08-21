#'
#' check monotonicity, etc
#'

library(goftest)
xx1 <- seq(0,20, length=1024)
xx2 <- seq(20, 100, by=0.1)
x <- c(xx1, xx2[-1])
for(n in c(1:10, 20, 50, 100, 200, 500, 1000, Inf)) {
  Fx <- pCvM(x, n)
  mindy <- min(diff(Fx))
  if(mindy < 0)
    stop(paste0("pCvM(x, ", n, ") is not monotone"))
}
