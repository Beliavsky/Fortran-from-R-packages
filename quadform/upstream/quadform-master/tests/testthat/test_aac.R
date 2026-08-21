## This file follows the structure of aaa.R in the free group package.

## Define a checker function, and call it at the end.  

test_that("Test suite aac.R, quad3.tform()",{


# First we need a helper function to create random complex matrices
# (NB: we cannot use the `cmvnorm` package because that depends on the
# `emulator` package):

tester <- function(a,b,SMALL = 1e-6){expect_true(all(abs(a-b) < SMALL))}

checker_mat <- function(M,l,r){
    tester(quad3.tform(M,l,r)    , l %*% M %*% ht(r))
    tester(quad3.tform_ab(M,l,r) , quad3.tform(M,l,r))
    tester(quad3.tform_bc(M,l,r) , quad3.tform(M,l,r))
}

mat_r <- function(row,col){ matrix(rnorm(row*col),row,col)}
mat_c <- function(row,col){ matrix(rnorm(row*col)+1i*rnorm(row*col),row,col)}

checker <- function(a,b,c){
    M <- mat_r(a,a)
    l <- mat_r(b,a)
    r <- mat_r(c,a)
    checker_mat(M,l,r)

    M <- mat_c(a,a)
    l <- mat_c(b,a)
    r <- mat_c(c,a)
    checker_mat(M,l,r)

}

checker(20,2,5)
checker(20,5,2)



})

