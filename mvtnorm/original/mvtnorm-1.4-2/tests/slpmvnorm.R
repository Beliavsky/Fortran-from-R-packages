
library("mvtnorm")

set.seed(290875)

if (require("numDeriv", quietly = TRUE)) {


chk <- function(...) stopifnot(isTRUE(all.equal(...)))

J <- 5

chks <- function(dg, tolerance = .Machine$double.eps^(1 / 4)) {

    prm <- runif(J * (J + c(-1, 1)[dg + 1L]) / 2)
    L <- ltMatrices(prm, diag = dg)
    a <- matrix(-2, nrow = J, ncol = 1)
    b <- matrix( 2, nrow = J, ncol = 1)
    M <- 10000L

    l <- function(x) {
        x <- ltMatrices(x, diag = dg)
        lpmvnorm(a, b, chol = x, M = M, seed = 29)
    }

    s <- function(x)
        slpmvnorm(a, b, chol = x, M = M, seed = 29)

    rl <- l(L)
    rs <- s(L)
    chk(rl, rs$logLik)
    if (dg) {
        chk(rs$chol, ltMatrices(grad(l, unclass(L)), diag = dg), 
            tolerance = tolerance)
    } else {
        chk(c(Lower_tri(rs$chol)), grad(l, unclass(L)), 
            tolerance = tolerance)
    }

    l <- function(x) {
        x <- ltMatrices(x, diag = dg)
        lpmvnorm(a, b, invchol = x, M = M, seed = 29)
    }

    s <- function(x)
        slpmvnorm(a, b, invchol = x, M = M, seed = 29)

    rl <- l(L)
    rs <- s(L)
    chk(rl, rs$logLik)
    if (dg) {
        chk(rs$invchol, ltMatrices(grad(l, unclass(L)), diag = dg), 
            tolerance = tolerance)
    } else {
        chk(c(Lower_tri(rs$invchol)), grad(l, unclass(L)), 
            tolerance = tolerance)    }

    l <- function(x)
        lpmvnorm(a, b, mean = x, chol = L, M = M, seed = 29)
    s <- function(x)
        slpmvnorm(a, b, mean = x, chol = L, M = M, seed = 29)

    x <- numeric(J)
    rl <- l(x)
    rs <- s(x)
    chk(rl, rs$logLik)
    chk(grad(l, x), c(rs$mean), tolerance = tolerance)
    x <- 1:J
    rl <- l(x)
    rs <- s(x)
    chk(rl, rs$logLik)
    chk(grad(l, x), c(rs$mean), tolerance = tolerance)

    l <- function(x)
        lpmvnorm(a, b, invcholmean = x, chol = L, M = M, seed = 29)
    s <- function(x)
        slpmvnorm(a, b, invcholmean = x, chol = L, M = M, seed = 29)

    x <- numeric(J)
    rl <- l(x)
    rs <- s(x)
    chk(rl, rs$logLik)
    chk(grad(l, x), c(rs$invcholmean), tolerance = tolerance)
    x <- 1:J
    rl <- l(x)
    rs <- s(x)
    chk(rl, rs$logLik)
    chk(grad(l, x), c(rs$invcholmean), tolerance = tolerance)
}

chks(TRUE)
chks(FALSE)

### check scores for conditional distributions

.cmvnorm <- function(invchol, which_given, given) {
    L <- invchol
    J <- dim(L)[2L]
    tmp <- matrix(0, ncol = ncol(given), nrow = J - length(which_given))
    center <- Mult(L, rbind(given, tmp))
    center <- center[-which_given,,drop = FALSE]
    L <- L[,-which_given]
    return(list(center = center, invchol = L))
}

J <- (cJ <- 4) + (dJ <- 4)
N <- 3
M <- 30
ltM <- function(x) ltMatrices(x, diag = FALSE, byrow = TRUE)
ltD <- function(x) ltMatrices(x, diag = TRUE, byrow = TRUE)
prm <- matrix(runif(J * (J - 1) / 2 * N), ncol = N)
L <- ltM(prm)

obs <- matrix(rnorm(J * N), ncol = N)
lwr <- -abs(obs)
upr <- abs(obs)

w <- matrix(runif((dJ - 1) * M), ncol = M)

### without scaling, diag(L) == 1
## score wrt L
j <- 1:cJ
ll <- function(x) {
    L <- ltM(x)
    cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- ltM(matrix(grad(ll, unclass(L)), ncol = N))
diagonals(a) <- 0

cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
             invchol = cd$invchol, w = w)

chk(as.array(a[,-j]), as.array(s$invchol), check.attributes = FALSE)

## score wrt obs
ll <- function(x) {
    cd <- .cmvnorm(invchol = L, which_given = j, given = x)
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, obs[1:cJ,,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

tmp0 <- solve(cd$invchol, s$mean, transpose = TRUE)
aL <- as.array(L)[-(1:cJ), 1:cJ,]
lst <- tmp0[rep(1:dJ, cJ),,drop = FALSE]
dobs <- -margin.table(aL * array(lst, dim = dim(aL)), 2:3)

chk(a, dobs, check.attributes = FALSE)

## score wrt lower
ll <- function(x) {
    cd <- .cmvnorm(invchol = L, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(x, upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, lwr[-(1:cJ),,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$lower, check.attributes = FALSE)

## score wrt upper
ll <- function(x) {
    cd <- .cmvnorm(invchol = L, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], x, center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, upr[-(1:cJ),,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$upper, check.attributes = FALSE)

### after scaling
LD <- invcholD(L)

## score wrt LD (!)
# use center
ll <- function(x) {
    LD <- ltD(x)
    cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w, logLik = TRUE)
}

a1 <- ltD(matrix(grad(ll, unclass(LD)), ncol = N))
# use mean
ll <- function(x) {
    LD <- ltD(x)
    cd <- cond_mvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], mean = cd$mean, 
            invchol = cd$invchol, w = w, logLik = FALSE)
}

a2 <- ltMatrices(matrix(grad(function(...) sum(ll(...)), unclass(LD)), ncol = N), 
                diag = TRUE, byrow = TRUE)

chk(a1, a2)

cd <- cond_mvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s1 <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], mean = cd$mean, 
             invchol = cd$invchol, w = w)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s2 <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
             invchol = cd$invchol, w = w)

### needs to be FALSE
# all.equal(s1, s2)

chk(a1[,-j], s2$invchol, check.attributes = FALSE)

# score wrt obs
ll <- function(x) {
    cd <- .cmvnorm(invchol = LD, which_given = j, given = x)
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, obs[1:cJ,,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, 
        w = w)

tmp0 <- solve(cd$invchol, s$mean, transpose = TRUE)
aL <- as.array(LD)[-(1:cJ), 1:cJ,]
lst <- tmp0[rep(1:dJ, cJ),,drop = FALSE]
dobs <- -margin.table(aL * array(lst, dim = dim(aL)), 2:3)

chk(a, dobs, check.attributes = FALSE)

## score wrt lower
ll <- function(x) {
    cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(x, upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, lwr[-(1:cJ),,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$lower, check.attributes = FALSE)

## score wrt upper
ll <- function(x) {
    cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], x, center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, upr[-(1:cJ),, drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$upper, check.attributes = FALSE)

### one-dimensional conditional distribution

J <- (cJ <- 4) + (dJ <- 1)
prm <- matrix(runif(J * (J - 1) / 2 * N), ncol = N)
L <- ltM(prm)

obs <- matrix(rnorm(J * N), ncol = N)
lwr <- -abs(obs)
upr <- abs(obs)

w <- matrix(runif((dJ - 1) * M), ncol = M)

### without scaling, diag(L) == 1
## score wrt L not needed (a constant)
j <- 1:cJ

## score wrt obs
ll <- function(x) {
    cd <- .cmvnorm(invchol = L, which_given = j, given = x)
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, obs[1:cJ,,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

tmp0 <- solve(cd$invchol, s$mean, transpose = TRUE)
aL <- as.array(L)[-(1:cJ), 1:cJ,]
lst <- tmp0[rep(1:dJ, cJ),,drop = FALSE]
dobs <- -aL * array(lst, dim = dim(aL))

chk(a, dobs, check.attributes = FALSE)

## score wrt lower
ll <- function(x) {
    cd <- .cmvnorm(invchol = L, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(x, upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, lwr[-(1:cJ),,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$lower, check.attributes = FALSE)

## score wrt upper
ll <- function(x) {
    cd <- .cmvnorm(invchol = L, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], x, center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, upr[-(1:cJ),,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = L, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$upper, check.attributes = FALSE)

### after scaling
LD <- invcholD(L)

## score wrt LD (!)
# use center
ll <- function(x) {
    LD <- ltD(x)
    cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w, logLik = TRUE)
}

a1 <- ltD(matrix(grad(ll, unclass(LD)), ncol = N))
# use mean
ll <- function(x) {
    LD <- ltD(x)
    cd <- cond_mvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], mean = cd$mean, 
            invchol = cd$invchol, w = w, logLik = FALSE)
}

a2 <- ltMatrices(matrix(grad(function(...) sum(ll(...)), unclass(LD)), ncol = N), 
                diag = TRUE, byrow = TRUE)

chk(a1, a2)

cd <- cond_mvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s1 <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], mean = cd$mean, 
             invchol = cd$invchol, w = w)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s2 <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
             invchol = cd$invchol, w = w)

### needs to be FALSE
# all.equal(s1, s2)

chk(a1[,-j], s2$invchol, check.attributes = FALSE)

# score wrt obs
ll <- function(x) {
    cd <- .cmvnorm(invchol = LD, which_given = j, given = x)
    lpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, obs[1:cJ,,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, 
        w = w)

tmp0 <- solve(cd$invchol, s$mean, transpose = TRUE)
aL <- as.array(LD)[-(1:cJ), 1:cJ,]
lst <- tmp0[rep(1:dJ, cJ),,drop = FALSE]
dobs <- -(aL * array(lst, dim = dim(aL)))

chk(a, dobs, check.attributes = FALSE)

## score wrt lower
ll <- function(x) {
    cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(x, upr[-j,,drop = FALSE], center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, lwr[-(1:cJ),,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$lower, check.attributes = FALSE)

## score wrt upper
ll <- function(x) {
    cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[1:cJ,,drop = FALSE])
    lpmvnorm(lwr[-j,,drop = FALSE], x, center = cd$center, 
            invchol = cd$invchol, w = w)
}
a <- matrix(grad(ll, upr[-(1:cJ),,drop = FALSE]), ncol = N)

cd <- .cmvnorm(invchol = LD, which_given = j, given = obs[j,,drop = FALSE])
s <- slpmvnorm(lwr[-j,,drop = FALSE], upr[-j,,drop = FALSE], center = cd$center, invchol = cd$invchol, w = w)

chk(a, s$upper, check.attributes = FALSE)

### check scores for extremely small likelihoods
### use independence model as ground truth
J <- 10

yl <- matrix(-as.double(rep(1, J) / 2000), ncol = 1)
yr <- abs(yl)

w <- matrix(.5, nrow = J - 1, ncol = 1)

## L = diag(J)
L <- ltMatrices(rep(1, J * (J + 1) / 2), diag = TRUE)
L[] <- 0
diagonals(L) <- 1

### constant b/c prob < .Machine$double.eps
lpmvnorm(lower = yl, upper = yr, invchol = L, w = w)

s <- slpmvnorm(lower = yl, upper = yr, invchol = L, w = w)[c("logLik", "invchol")]

chk(0, unique(c(diagonals(s$invchol))))
}
