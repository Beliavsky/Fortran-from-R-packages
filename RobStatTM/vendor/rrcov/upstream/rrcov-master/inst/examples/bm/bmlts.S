##### bm_lts #####
# Benchmark for rrcov::ltsReg() and MASS::ltsreg() on several n and p.
#
#   V.Todorov: 16.08.2004
# 
# For each n and p (specified by the arrays <an> and <ap>) a data set 
#   is generated (see the function gendata() for description of the model)
#   and LTS is computed by calling ltsReg(). All defaults are accepted.
#
# Each test (for given n and p) is performed several times (specified
#   by the parameter <nrep>) and the result is averaged.
#   
# Input argument: 
#   nrep: Number of times the tests are executed (defaults to 3)
#   eps: contamination proportion - see function gendata() (defaults to 0.4)
#   method: "rrcov" or "MASS"
#
bmlts <- function(nrep=1, eps=0.4, method=c("rrcov","MASS", "S")){

    method <- match.arg(method)
    if(method == "rrcov")
        library(rrcov)

    library(MASS)
 
    ap <- c(2, 3, 5, 10)
    an <- c(100, 500, 1000, 10000, 50000)

    set.seed(0)
    btime <- proc.time()

    cat("*** Benchmark for R/S: LTS ***")
    cat("\nThe results are averaged on nrep = ", nrep, "\n")

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    cat("     n   p       Time\n")
    cat("=====================\n")

    for(i in 1:length(an)) {
        for(j in 1:length(ap)) {
            n <- as.integer(an[i])
            p <- as.integer(ap[j])
            if(5*p <= n){
                a <- gendataLTS(n, p, eps)
		X <- a$x
		Y <- a$y
                
                # set mcd=FALSE - we want to time only the LTS algorithm
                ptm <- proc.time()
                for(k in 1:nrep){
                    if(method == "MASS")
                        ltsreg(X,Y)
                    else if(method == "rrcov")
                        ltsReg(X, Y, mcd = FALSE)
		    else
                        ltsreg(X, Y, mcd = FALSE)
                }
                xtime <- proc.time() - ptm
                xtime <- xtime[1]/nrep
#                cat(sprintf("%6d %3d %10.2f\n", n, p, xtime))
                cat(n, p, format(xtime, nsmall=2), "\n")
            }
        } 
    }
    
    tottime <- proc.time() - btime
    cat("=====================\n")
    cat("Total time: ", tottime[1], "\n")
}

#### gendata() ####
# Generates a data set with bad leverage points (outliers in x-space) 
# n observations in p dimensions acording to the model:
#   yi = Xi1+Xi2+...+ei
# where ei - N(0,1) is the error term, Xi,j for j=1...p-1 - N(0,100) are 
# the non-trivial explanatory variables and xip is the intercept term.
# The outliers in the x-space are introduced by replacing eps. percent of
# xi1 by values distributed as N(100,100).
#
# Defaults: eps=0
#
gendataLTS <- function(n, p, eps=0){

    if(eps < 0 || eps >= 0.5)
        stop(message="eps must be in [0,0.5)")

    p <- p-1
    x <- matrix(data=rnorm(n*p,0,100), nrow=n, ncol=p)
    y <-rowSums(x) + 1 + rnorm(n, 0, 1)
    
    nbad <- as.integer(eps * n)
    xind <- sort(sample(n,nbad))
    xbad <- rnorm(nbad,100,100)
    for(i in 1:nbad){
        x[xind[i],1] <- xbad[i] 
    }
    list(x=x, y=y, xind=xind)
}
