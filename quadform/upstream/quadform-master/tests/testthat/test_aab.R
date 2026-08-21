## This file follows the structure of aaa.R in the free group package.

## Define a checker function, and call it at the end.  

test_that("Test suite aab.R, quad3.form()",{


# First we need a helper function to create random complex matrices
# (NB: we cannot use the `cmvnorm` package because that depends on the
# `emulator` package):

tester <- function(a,b,SMALL = 1e-6){expect_true(all(abs(a-b) < SMALL))}

checker_mat <- function(M,l,r){
    tester(quad3.form(M,l,r)     , ht(l) %*% M %*% r)
    tester(quad3.form_ab(M,l,r)  , quad3.form(M,l,r))
    tester(quad3.form_bc(M,l,r)  , quad3.form(M,l,r))
}

mat_r <- function(row,col){ matrix(rnorm(row*col),row,col)}
mat_c <- function(row,col){ matrix(rnorm(row*col)+1i*rnorm(row*col),row,col)}

checker <- function(a,b,c){
    M <- mat_r(a,a)
    l <- mat_r(a,b)
    r <- mat_r(a,c)
    checker_mat(M,l,r)

    M <- mat_c(a,a)
    l <- mat_c(a,b)
    r <- mat_c(a,c)
    checker_mat(M,l,r)

}

checker(20,2,5)
checker(20,5,2)



})

