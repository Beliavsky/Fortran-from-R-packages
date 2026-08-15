#' @export
TP <- function(x,theta) {
if (x == 0) return(1)
Stirling2 <- function(n,m){
if (0 > m || m > n) stop("m must be in 0..n")
k <- 0:m
sig <- rep(c(1,-1)*(-1)^m, length= m+1)
ga <- gamma(k+1)
round(sum(sig * k^n /(ga * rev(ga))))
}
vBell <- numeric()
for(j in 1:x) vBell[j] <- Stirling2(x,j)*theta^j
return(sum(vBell))

}
