okLongDouble <- function(lambda = 999, verbose = 0L, tol = 1e-15)
{
    ## Author: Martin Maechler, Oct 17, 2019
    stopifnot(length(lambda) == 1L, length(verbose) == 1L, length(tol) == 1L)
    .Call(C_chk_LDouble, lambda, verbose, tol) # -> ../src/ppois-direct.c
}

ppoisD <- function(q, lambda, all.from.0 = TRUE, verbose = 0L)
{
  ## Purpose: ppois() via direct computation ("sum(dpois)"
  ## ----------------------------------------------------------------------
  ## Arguments: as ppois(); not yet "lower.tail" or "log.p"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  2 Mar 2004; .C -> .Call(): 2019-08-15
    if(all.from.0) stopifnot(length(q) == 1L)
    .Call(C_ppoisD, as.double(q), lambda, all.from.0, verbose) # -> ../src/ppois-direct.c
}


ppoisErr <- function(lambda, ppFUN = ppoisD, iP = 1e-15,
                     xM = qpois(iP, lambda=lambda, lower.tail=FALSE),
                     verbose = FALSE) {
    ## Purpose: Given lambda, find "worst case 'x' for ppois(x, lambda)
    ##		and return the relative error of ppois()
    ## ----------------------------------------------------------------------
    ## Arguments: lambda:  parameter lambda
    ##		  iP : 1 - p; where maximal x is the p-th quantile;
    ##			as we find, the worst x we need is for p ~= 0.229
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 1 Mar 2004;  Aug. 2019 ('ppFUN' as argument)

    stopifnot(is.function(ppFUN), length(ff <- formals(ppFUN)) >= 2
	    , names(ff)[2:3] == c("lambda", "all.from.0")
	      )
    ## NB: x will be large for large lambda! (R >= 3.5.0: uses 'compact' ALTREP)
    x <- 0:xM
    i0 <- which.max(abs((p1x <- ppois(x , lambda = lambda)) -
                        (p0x <- ppFUN(xM, lambda = lambda, all.from=TRUE))))
    x0 <-   x[i0] # = i0 - 1
    p1 <- p1x[i0] # = ppois(x0, lambda)
    p2 <- p0x[i0] # ~= sum(dpois(x0:0, lambda)) # supposedly "true" value
    if(verbose)
        cat("x0=",formatC(x0),", ppois()=",formatC(p1),"\n")
    ## return rel.error *and* 'x0' as attribute :
    structure(1 - p1/p2, x = c(x0=x0, xM=xM))
}
