## -*- truncate-lines: t; -*-

## library("tinytest")
## library("neighbours")

steps <- if (Sys.getenv("ES_PACKAGE_TESTING_73179826243954") != "")
             50000 else 1000


## ------ logical: constant k
N <- neighbourfun(type = "logical",
                  kmin = 3,
                  kmax = 3)

x <- logical(10)
x[1:3] <- TRUE

for (i in 1:steps)
    expect_equal(3, sum(N(x)))



## ------ logical: switch 1 element
N <- neighbourfun(type = "logical")
x <- logical(10)
for (i in 1:steps) {
    xn <- N(x)
    expect_equal(sum(xn != x), 1)
    x <- xn
}



## ------ logical: switch 3 elements
N <- neighbourfun(type = "logical", stepsize = 3)
x <- logical(10)
for (i in 1:steps) {
    xn <- N(x)
    expect_equal(sum(xn != x), 3)
    x <- xn
}



## ------ logical: switch 3 elements
N <- neighbourfun(type = "logical", stepsize = 3, n = 10)
x <- logical(10)
for (i in 1:steps) {
    xn <- N(x)
    expect_equal(sum(xn != x), 3)
    x <- xn
}

## ------ logical: between 0 and 5 elements
N <- neighbourfun(type = "logical", stepsize = 1, kmin = 0, kmax = 5)
x <- logical(10)
for (i in 1:steps) {
    xn <- N(x)
    expect_true(sum(xn) <= 5)
    x <- xn
}





## ------ numeric: change 1 element, fixed stepsize
x.min <- -0.05
x.max <-  0.05
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.01,
                  random = FALSE)
x <- rep(1/25, 25)
for (i in 1:steps) {
    xn <- N(x)
    expect_true(all(xn <= x.max))
    expect_true(all(xn >= x.min))
    x <- xn
}


## ------ numeric: change 1 element
x.min <- -0.05
x.max <-  0.05
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.015)
x <- rep(1/25, 25)
for (i in 1:steps) {
    xn <- N(x)
    expect_true(all(xn <= x.max))
    expect_true(all(xn >= x.min))
    x <- xn
}



## ------ numeric: 'active' with min/max vectors
n <- 20
x.min <- rep(-0.1, n)
x.max <- seq(0.1, 0.2, length.out = n)
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.1,
                  active = 1:5)
x <- rep(1/n, n)
for (i in 1:steps) {
    xn <- N(x)
    expect_true(all(xn <= x.max + sqrt(.Machine$double.eps)))
    expect_true(all(xn + sqrt(.Machine$double.eps) >= x.min))
    expect_true(all(x[6:20] == xn[6:20]))
    x <- xn
}




## ------ numeric: updating
x.min <- -0.05
x.max <-  0.05
A <- array(rnorm(100*25), dim = c(100,25))
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.015,
                  update = "Ax",
                  A = A)
x <- rep(1/25, 25)
attr(x, "Ax") <- A %*% x

for (i in 1:steps) {
    xn <- N(x, A)
    expect_true(all(xn <= x.max))
    expect_true(all(xn >= x.min))
    x <- xn
}

x <- rep(1/25, 25)
attr(x, "Ax") <- A %*% x
for (i in 1:steps) {
    x <- N(x, A)
}
expect_equivalent(A %*%x, attr(x, "Ax"))




## ------ numeric: updating, active
x.min <- -0.1
x.max <-  0.1
active <- 1:5
A <- array(rnorm(100*25), dim = c(100,25))
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.015,
                  update = "Ax",
                  A = A,
                  active = 1:5)
x <- rep(1/25, 25)
attr(x, "Ax") <- A %*% x

for (i in 1:steps) {
    xn <- N(x, A)
    expect_true(all(xn <= x.max))
    expect_true(all(xn >= x.min))
    x <- xn
}

x0 <- x <- rep(1/25, 25)
attr(x, "Ax") <- A %*% x
for (i in 1:steps) {
    x <- N(x, A)
}
for (i in 1:steps) {
    xn <- N(x)
    expect_true(all(xn <= x.max + sqrt(.Machine$double.eps)))
    expect_true(all(xn + sqrt(.Machine$double.eps) >= x.min))
    expect_true(all(x[6:25] == xn[6:25]))
    x <- xn
}
expect_true(all(x[6:25] == x0[6:25]))





## ------ numeric: updating, with a constant added
x.min <- -0.01
x.max <-  0.01
tol <- 1e-14
A <- array(rnorm(100*25), dim = c(100, 25))
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.015,
                  update = "Ax",
                  A = A)
x0 <- rep(1/25, 25)
x <- rep(0, 25)
attr(x, "Ax") <- A %*% (x + x0)

for (i in 1:steps) {
    xn <- N(x, A)
    expect_true(all(xn - tol <= x.max))
    expect_true(all(xn + tol >= x.min))
    x <- xn
}
expect_equivalent(A %*% (x + x0),   attr(x, "Ax"))



## ------ numeric: change 1 element, fixed stepsize, min/max budget
x.min <- -0.20
x.max <-  0.20
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.01,
                  random = FALSE,
                  sum = c(0.8, 1.2))
x <- rep(1/25, 25)
for (i in 1:steps) {
    xn <- N(x)
    expect_true(all(round(xn - x, 8) %in% c(0, -0.01, 0.01)))
    expect_true(all(xn <= x.max))
    expect_true(all(xn >= x.min))
    expect_true(sum(xn) <= 1.2 + sqrt(.Machine$double.eps))
    expect_true(sum(xn) >= 0.8 - sqrt(.Machine$double.eps))
    expect_true(sum(xn != x) <= 1L)
    x <- xn
}


## ------ numeric: change 1 element, random stepsize, min/max budget
x.min <- -0.20
x.max <-  0.20
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.01,
                  random = TRUE,
                  sum = c(0.8, 1.2))
x <- rep(1/25, 25)
for (i in 1:steps) {
    xn <- N(x)
    expect_true(all(xn <= x.max))
    expect_true(all(xn >= x.min))
    expect_true(sum(xn) <= 1.2 + sqrt(.Machine$double.eps))
    expect_true(sum(xn) >= 0.8 - sqrt(.Machine$double.eps))
    expect_true(sum(xn != x) <= 1L)
    x <- xn
}


## ------ numeric: change 1 element, random stepsize, no budget
x.min <- -0.20
x.max <-  0.20
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.01,
                  random = TRUE,
                  sum = FALSE)
x <- rep(1/25, 25)
for (i in 1:steps) {
    xn <- N(x)
    expect_true(sum(xn) != sum(x))
    expect_true(all(xn >= x.min - sqrt(.Machine$double.eps)))
    expect_true(all(xn <= x.max + sqrt(.Machine$double.eps)))
    expect_true(sum(xn != x) <= 1L)
    x <- xn
}





## ------ numeric: change 1 element, random stepsize, sum == 0
x.min <- -0.05
x.max <-  0.05
N <- neighbourfun(type = "numeric",
                  min = x.min,
                  max = x.max,
                  stepsize = 0.01)
x <- runif(20, max = x.max)
x <- x - mean(x)

for (i in 1:steps) {
    xn <- N(x)
    expect_equivalent(sum(xn), 0)
    expect_true(all(xn >= x.min))
    expect_true(all(xn <= x.max))
    x <- xn
}




## ------ permute: with list
N <- neighbourfun(type = "permute", stepsize = 2)
x <- as.list(letters[1:5])
xx <- N(x)
expect_equal(mode(xx), "list")
expect_equal(sum(unlist(x) == unlist(xx)), 3)
expect_equal(sum(unlist(x) == sort(unlist(xx))), 5)








## TODO
##   create example of find
## LSopt <- function(OF, algo = list(), ...) {
##     xc  <- algo$x0
##     xcF <- OF(xc, ...)
##     for (s in seq_len(algo$nI)) {
##         xn <- algo$neighbour(xc, ...)
##         xnF <- OF(xn, ...)
##         if (xnF <= xcF) {
##             xc  <- xn
##             xcF <- xnF
##         }
##     }
##     list(xbest = xc, OFvalue = xcF)
## }

## target <- runif(20)



## deviation <- function(x, target) {
##     xx <- x - target
##     crossprod(xx)
## }

## sol <- LSopt(deviation,
##              list(neighbour = neighbourfun(type = "numeric",
##                                            ## length = 20,
##                                            stepsize = 0.03),
##                   x0 = runif(length(target)),
##                   nI = 500000,
##                   printBar = FALSE),
##              target = target)

## sol$xbest - target
