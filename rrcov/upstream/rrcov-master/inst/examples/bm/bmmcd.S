##### bmmcd #####
# Benchmark for rrcov::covMcd() and MASS::cov.mcd() on several n and p.
#
#   V.Todorov: 16.08.2004
#
# For each n and p (specified by the arrays <an> and <ap>) a data set
#   is generated (see the function gendata() for description of the model)
#   and MCD is computed by calling covMcd(). All defaults are accepted
#
# Each test (for given n and p) is performed several times (specified
#   by the parameter <nrep>) and the result is averaged.
#
# Input argument:
#   nrep: Number of times the tests are executed (defaults to 3)
#   eps: contamination proportion - see function gendata() (defaults to 0.4)
#   method: "rrcov" or "MASS"
#
bmmcd <- function(nrep=1, eps=0.40, method=c("rrcov", "MASS", "S", "Smcd", "mrcd", "detmcd")){

    method <- match.arg(method)

    if(method == "rrcov"){
       library(rrcov)
    }

    if(method == "S" || method == "Smcd"){
       library(Robust)
    }
	
    library(MASS)

    ap <- c(2, 5, 10, 20, 30)
    an <- c(100, 500, 1000, 10000, 50000)

    set.seed(0)
    btime <- proc.time()

    cat("\n*** Benchmark for R/S MCD estimators ***")
    cat("\nImplementation=",method)
    cat("\nThe results are averaged on nrep = ",nrep," runs.\n")

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    cat("     n   p       Time\n")
    cat("=====================\n")

    for(i in 1:length(an)) {
        for(j in 1:length(ap)) {
            n <- as.integer(an[i])
            p <- as.integer(ap[j])
            if(5*p <= n){
                xx <- gendataMCD(n, p, eps)
		X <- as.data.frame(xx$X)

            ptm <- proc.time()
            for(k in 1:nrep){
                if(method == "MASS")
                    cov.mcd(X)
                else if(method == "rrcov")
                    covMcd(X)
                else if(method == "S")	
                    covRob(X)
                else if(method == "Smcd")	
                    covRob(X,estim="mcd")
                else if(method == "detmcd")	
                    CovMcd(X, nsamp="deterministic")
                else if(method == "mrcd")	
                    CovMrcd(X)
                else
                    stop("No method specified")
            }
            xtime <- proc.time() - ptm
            xtime <- xtime[1]/nrep

#                cat(sprintf("%6d %3d %10.2f\n", n, p, xtime))
                cat(format(n), format(p), format(xtime, nsmall=2),"\n")

            }
        }
    }

    tottime <- proc.time() - btime
    cat("=====================\n")
    cat("Total time: ", tottime[1], "\n")
}

#### gendata() ####
# Generates a location contaminated multivariate
# normal sample of n observations in p dimensions
#    (1-eps)*Np(0,Ip) + eps*Np(m,Ip)
# where
#    m = (b,b,...,b)
# Defaults: eps=0 and b=10
#
gendataMCD <- function(n,p,eps=0,b=10){

    if(missing(n) || missing(p))
        stop("Please specify (n,p)")
    if(eps < 0 || eps >= 0.5)
        stop(message="eps must be in [0,0.5)")
    X <- mvrnorm(n,rep(0,p),diag(1,nrow=p,ncol=p))
    nbad <- as.integer(eps * n)
    if(nbad > 0){
        Xbad <- mvrnorm(nbad,rep(b,p),diag(1,nrow=p,ncol=p))
        xind <- sample(n,nbad)
        X[xind,] <- Xbad
    }
    list(X=X, xind=xind)
}
