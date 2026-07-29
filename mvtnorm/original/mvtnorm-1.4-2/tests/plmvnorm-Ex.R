
library("mvtnorm")

pmvnorm <- function(lower = -Inf, upper = Inf, mean = rep(0, length(lower)), 
                    corr = NULL, sigma = NULL, algorithm = GenzBretz(), ...) {

    if (!inherits(algorithm, "GenzBretz"))
        return(mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, 
                         sigma = sigma, corr = corr, algorithm = algorithm, ...))

    args <- mvtnorm:::checkmvArgs(lower = lower, upper = upper, mean = mean, 
                                  sigma = sigma, corr = corr)
    if (args$uni)
        return(mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, 
                         sigma = sigma, corr = corr, algorithm = algorithm, ...))

    if (!is.null(args$corr)) args$sigma <- args$corr

    Chol <- try(chol(args$sigma), silent = TRUE)
    if (inherits(Chol, "try-error")) 
       return(mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, 
                         sigma = sigma, corr = corr, algorithm = algorithm, ...))
    Chol <- matrix(t(Chol)[lower.tri(Chol, diag = TRUE)], ncol = 1)
    Chol <- ltMatrices(Chol, diag = TRUE, byrow = FALSE)

    args$chol <- Chol
    
    M <- algorithm$maxpts
    if (require("qrng", quietly = TRUE)) {
        w <- t(ghalton(M, d = length(args$lower) - 1))
    } else {
        w <- NULL
    }
    args$w <- w
    args$M <- M
    args$seed <- 290875
    args$logLik <- FALSE
    args$corr <- args$sigma <- args$uni <- NULL
    args$lower <- matrix(args$lower, ncol = 1)
    args$upper <- matrix(args$upper, ncol = 1)

    ret <- exp(do.call("lpmvnorm", args))
    ret[ret < .Machine$double.eps] <- 0
    ret
}

try(source("regtest-TVPACK.R", echo = TRUE))
try(source("test-noisy-root.R", echo = TRUE))
try(source("bugfix-tests.R", echo = TRUE))
