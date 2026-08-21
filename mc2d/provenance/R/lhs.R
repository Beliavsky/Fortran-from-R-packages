#' Random Latin Hypercube Sampling
#'
#' Creates a Latin Hypercube Sample (LHS) from the specified distribution.
#'
#' @note The resulting LHS is a latin hypersquare sampling: the LHS structure is provided
#' only in the first 2 dimensions.
#' It is not possible to use truncated distributions directly with \code{\link{rtrunc}}.
#' Use \code{\link{mcstoc}} with \samp{lhs=TRUE} and \samp{rtrunc=TRUE} for that purpose.
#' The \dots arguments will be recycled.
#'
#' @param distr the function for generating random samples, or its name as a character string.
#'   If \samp{distr} is \samp{"rdist"}, the quantile function \samp{"qdist"} must exist with
#'   argument \samp{p} as a vector of probabilities.
#' @param nsv the number of rows of the final matrix.
#' @param nsu the number of columns of the final matrix.
#' @param nvariates the number of variates.
#' @param ... all arguments to be passed to \samp{distr} except the size of the sample.
#' @return A \samp{nsv x nsu} matrix of random variates.
#' @note Adapted from a code by Rob Carnell (package \pkg{lhs}).
#' @seealso \code{\link{mcstoc}}
#' @keywords design
#' @examples
#' ceiling(lhs(runif, nsu=10, nsv=10) * 10)
#' @export
lhs <- function(distr="runif",nsv=ndvar(),nsu=ndunc(),nvariates=1,...)
{
    nsv
    nsu
    if (length(nsv) != 1 | length(nsu) != 1) stop("nsv and nsu may not be vectors")
    if (any(is.na(c(nsv, nsu)))) stop("nsv and nsu may not be NA or NaN")
    if (any(is.infinite(c(nsv, nsu)))) stop("nsv and nsu may not be infinite")
    if (floor(nsv) != nsv | nsv < 1) stop("nsv must be a positive integer\n")
    if (floor(nsu) != nsu | nsu < 1) stop("nsu must be a positive integer\n")

    arg <- list(...)

    if(!is.character(distr)) distr <- as.character(match.call()$distr)                     #retrieve the name of the function
    distr <- substr(distr, 2, 1000)                             #remove the r
    distr <- paste("q",distr,sep="")                            # add the q

    ranperm <- function(X, N) order(runif(N))

    P <- array(dim=c(nsv, nsu, nvariates))
    for(i in 1:nvariates) {
      P[,,i] <- apply(P[,,i,drop=FALSE], 2, ranperm, N = nsv)
      eps <- matrix(runif(nsv * nsu), nrow = nsv, ncol = nsu)
      P[,,i] <- (P[,,i] - 1 + eps) / nsv
      }
    return(as.vector(do.call(distr,c(list(p=P),arg))))
}
