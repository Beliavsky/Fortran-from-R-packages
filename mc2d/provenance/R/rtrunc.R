#' Random Truncated Distributions
#'
#' Provides samples from classical \R distributions and \samp{mc2d} specific
#' distributions truncated between \samp{linf} (excluded) and \samp{lsup} (included).
#'
#' The function 1) evaluates the \samp{p} values corresponding to \samp{linf} and
#' \samp{lsup} using \samp{pdistr}; 2) samples \samp{n} values using
#' \samp{runif(n, min=pinf, max=psup)}, and 3) takes the \samp{n} corresponding
#' quantiles from the specified distribution using \samp{qdistr}.
#'
#' All distributions (but \code{sample}) implemented in the \pkg{stats} library
#' could be used. The arguments in \dots should be named. Do not use \samp{log},
#' \samp{log.p} or \samp{lower.tail}.
#'
#' For discrete distributions, \samp{rtrunc} samples within \samp{(linf, lsup]}.
#'
#' \strong{Warning:} The method is flexible, but can lead to problems linked to rounding
#' errors in some extreme situations. The function checks that the values are in the
#' expected range and returns an error if not.
#'
#' @param distr A function providing random data or its name as character.
#'   The function \samp{rdistr} should have a \samp{qdistr} form (with argument \samp{p})
#'   and a \samp{pdistr} form (with argument \samp{q}). Example: \samp{rnorm}
#'   (has a \samp{qnorm} and a \samp{pnorm} form), \samp{rbeta}, \samp{rbinom},
#'   \samp{rgamma}, ...
#' @param n the size of the sample.
#' @param linf a vector of lower bounds.
#' @param lsup a vector of upper bounds, with \samp{lsup > linf} (strictly).
#' @param ... all arguments to be passed to \samp{pdistr} and \samp{qdistr}.
#' @return A vector of \samp{n} values.
#' @keywords distribution
#' @examples
#' rtrunc("rnorm", n = 10, linf = 0)
#' range(rtrunc(rnorm, n = 1000, linf = 3, lsup = 5, sd = 10))
#' ## Discrete distributions
#' range(rtrunc(rpois, 1000, linf = 2, lsup = 4, lambda = 1))
#' @export
rtrunc <- function(distr=runif, n, linf=-Inf, lsup=Inf,...)
{
    linf <- as.vector(linf)
	lsup <- as.vector(lsup)
	if(!is.character(distr)) distr <- as.character(match.call()$distr)          #retrieve the name of the function
    distr <- substr(distr, 2, 1000)                                             #remove the r

    if(any(linf >= lsup)) stop("linf should be < lsup")  #recycle vectors

    pfun <- get(paste("p",distr,sep=""),mode="function")

    pinf <- as.vector(pfun(q=linf,...))
    psup <- as.vector(pfun(q=lsup,...))


    p <- runif(n,min=pinf,max=psup)

    qfun <- get(paste("q",distr,sep=""),mode="function")

    res <- as.vector(qfun(p,...))
    # Some possible problem (check if you think to others)
    #
    res[pinf <= 0 & res > lsup] <- NaN          #ex: rtrunc("lnorm",10,linf=-2,lsup=-1)
    res[psup>=1 & res < linf] <- NaN          #ex: rtrunc("unif",10,linf=2,lsup=4,max=1)
    res[is.na(linf) | is.na(lsup)] <- NaN   #ex: rtrunc("norm",10,sd=-2)


    #Two tests for extreme situations. None Catch all possibilities. THe error is first to avoid the warning
    if(any(res <= linf | res > lsup, na.rm=TRUE)) stop("Error in rtrunc: some values are not in the expected range (maybe due to rounding errors)")
    if(isTRUE(all.equal(pinf,1)) | isTRUE(all.equal(psup,0)) ) warning("Warning: check the results from rtrunc. It may have reached rounding errors")

    return(res)
}
