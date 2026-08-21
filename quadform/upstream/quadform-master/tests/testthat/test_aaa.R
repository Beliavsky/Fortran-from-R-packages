## This file follows the structure of aaa.R in the free group package.

## Define a checker function, and call it at the end.  

test_that("Test suite aaa.R",{


# First we need a helper function to create random complex matrices
# (NB: we cannot use the `cmvnorm` package because that depends on the
# `emulator` package):

tester <- function(a,b,SMALL = 1e-6){expect_true(all(abs(a-b) < SMALL))}

checker <- function(M,x,y,x1,y1){
    tester(ht(x)                , Conj(t(x)))
    tester(ht(x)                , t(Conj(x)))
    tester(cprod(x)             , ht(x) %*% x)
    tester(cprod(x,x1)          , ht(x) %*% x1)
    tester(tcprod(x)            , x %*% ht(x))
    tester(tcprod(x,x1)         , x %*% ht(x1))
    tester(quad.form(M,x)       , ht(x) %*% M %*% x)
    tester(quad.form(solve(M),x), ht(x) %*% solve(M) %*% x)
    tester(quad3.form(M,x,x1)   , ht(x) %*% M %*% x1)
    tester(quad3.form.inv(M,x,x1)   , ht(x) %*% solve(M) %*% x1)
    tester(quad3.tform(M,y,y1)  , y %*% M %*% ht(y1))
    tester(quad.tform(M,y)      , y %*% M %*% ht(y))
    tester(quad.tform.inv(M,y)  , y %*% solve(M) %*% ht(y))
    tester(quad.diag(M,x)       , diag(ht(x) %*% M %*% x))
    tester(quad.diag(M,x)       , diag(quad.form(M,x)))
    tester(quad.tdiag(M,y)      , diag(y %*% M %*% ht(y)))
    tester(quad.tdiag(M,y)      , diag(quad.tform(M,y)))
    tester(quad3.diag(M,x,x1)   , diag(ht(x) %*% M %*% x1))
    tester(quad3.diag(M,x,x1)   , diag(quad3.form(M,x,x1)))
    tester(quad3.tdiag(M,y,y1)  , diag(y %*% M %*% ht(y1)))
    tester(quad3.tdiag(M,y,y1)  , diag(quad3.tform(M,y,y1)))
    tester(quad.trace(M,x)      , sum(quad.diag(M,x)))
    tester(quad.ttrace(M,y)     , sum(quad.tdiag(M,y)))
    tester(quad.trace(M,x)      , sum(diag(quad.form(M,x))))
    tester(quad.ttrace(M,y)     , sum(diag(quad.tform(M,y))))

    if(is.numeric(M)){# should be "is.real"
        M <- cprod(M) # to ensure positive-definiteness
        tester(quad.form.chol(ht(chol(M)),x), quad.form(M,x))
    }
}

mat_r <- function(row,col){ matrix(rnorm(row*col),row,col)}
mat_c <- function(row,col){ matrix(rnorm(row*col)+1i*rnorm(row*col),row,col)}


a <- 2
b <- 3  # fails if b=1
M <- mat_r(a,a)
x <- mat_r(a,b)
y <- mat_r(b,a)
x1 <- mat_r(a,b)
y1 <- mat_r(b,a)
checker(M,x,y,x1,y1)

M <- mat_c(a,a)
x <- mat_c(a,b)
y <- mat_c(b,a)
x1 <- mat_c(a,b)
y1 <- mat_c(b,a)
checker(M,x,y,x1,y1)



})

