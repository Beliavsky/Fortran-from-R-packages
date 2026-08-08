if (!require(mvtnorm)) {
    stop("package mvtnorm required")
}
if (!require(devtools)) {
    stop("package devtools required")
}

set.seed(123)

N <- 200
k <- 2
T <- 100

x.mean <- rep(0,k)
x.cov <- diag(k)
mu <- rnorm(k, 0, 1)
Omega <- diag(k)

X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
XB <- colSums(X * B)
log.p <- XB - log1p(exp(XB))
Y <- sapply(log.p, function(q) return(rbinom(1,T,exp(q))))

binary <- list(Y=Y, X=X, T=T)

usethis::use_data(binary, overwrite=TRUE)
