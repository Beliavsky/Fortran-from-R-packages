### Does maxControl stuff behave?
### 
### do not run unless 'NOT_CRAN' explicitly defined
### (Suggested by Sebastian Meyer and others)
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    message("skipping slow optimizer tests")
    q("no")
}
### test for:
### 1. create maxControl object
### 2. SGA_batchSize NULL
### 3. negative batch size
### 4. more than 1 batch size
### SG_clip: NULL, negative, more than one
###
### printing:
### * #of cols, rows

library(maxLik)
set.seed(3)

### ---------- create maxControl object
maxControl(tol=1e-4, lambdatol=1e-5, qrtol=1e-6, qac="marquardt",
           marquardt_lambda0=0.1, marquardt_lambdaStep=3, marquardt_maxLambda=1e10,
           nm_alpha=2, nm_beta=1, nm_gamma=4,
           sann_temp=5, sann_tmax=100, sann_randomSeed=1,
           SGA_momentum=0.9, Adam_momentum1=0.5, Adam_momentum2=0.55,
           SG_learningRate=0.5, SG_batchSize=10, SG_clip=1000, 
           SG_patience=7, SG_patienceStep=10,
           iterlim=10, printLevel=3)

### ---------- SG_batchSize
expect_silent(maxControl(SG_batchSize=NULL))
expect_error(maxControl(SG_batchSize=-1))  # should fail
expect_error(maxControl(SG_batchSize=2:3))  # should fail

expect_silent(maxControl(SG_clip=NULL))
expect_error(maxControl(SG_clip=-1))  # fails
expect_error(maxControl(SG_clip=2:3))  # fails

expect_error(maxControl(Adam_momentum1=NA))  # should fail w/'NA in Adam_momentum'


### ---------- printing ----------
### ---------- max.columns, max.rows ----------
loglik <- function(beta) {
   e <- y - X %*% beta
   -crossprod(e)
}
gradlik <- function(beta) {
   e <- y - X %*% beta
   l <- crossprod(e)
   g <- t(-2*t(X) %*% e)
   -g
}
## linear regression with many columns
X <- matrix(rnorm(20*15), 20, 15)
beta <- rep(1, ncol(X))
y <- X %*% beta + rnorm(20, sd=0.3)
m <- maxNR(loglik, gradlik, start=rep(1, ncol(X)), iterlim=1)

## print estimates + gradient, and hessian
## should print only 4 rows for estimates, 4 rows + 2 cols for Hessia
## should give message "reached getOption("max.cols") -- omitted 13 columns" etc
expect_stdout(print(summary(m, hessian=TRUE), max.rows=4, max.cols=2, digits=3),
              pattern=paste0('reached getOption\\("max.rows"\\) -- omitted 11 rows',
                             '.*',
                             'reached getOption\\("max.cols"\\) -- omitted 13 columns',
                             '.*',
                             'reached getOption\\("max.rows"\\) -- omitted 11 rows')
              )
