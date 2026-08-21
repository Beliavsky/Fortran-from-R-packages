#' Creates a Stochastic mcnode Object Using a Probability Tree
#'
#' Builds an \samp{mcnode} as a mixture of \samp{mcnode} objects.
#'
#' \samp{mcswitch} may be either:
#' \itemize{
#'   \item a vector of weights (need not sum to one, must be nonneg and not all zero).
#'     Length must equal the number of elements in \samp{mcvalues}.
#'   \item a \samp{"0 mcnode"} to build any type of node.
#'   \item a \samp{"V mcnode"} to build a \samp{"V"} or \samp{"VU"} node.
#'   \item a \samp{"U mcnode"} to build a \samp{"U"} or \samp{"VU"} node.
#'   \item a \samp{"VU mcnode"} to build a \samp{"VU"} node.
#' }
#'
#' Names of elements in \samp{mcvalues} should correspond to values in \samp{mcswitch},
#' specified as character. Elements are evaluated only if needed.
#'
#' @param mcswitch a vector of probabilities/weights or an \samp{mcnode}.
#' @param mcvalues a named list of \samp{mcnode}s, \samp{mcdata} calls, or \samp{mcstoc} calls.
#'   Each element should be or lead to a compatible \samp{mcnode}.
#' @param type the type of \samp{mcnode} to be built. By default, a \samp{"V"} node.
#' @param nsv the number of simulations in the variability dimension.
#' @param nsu the number of simulations in the uncertainty dimension.
#' @param nvariates the number of variates of the final \samp{mcnode}.
#' @param outm the default output for multivariate nodes. See \code{\link{outm}}.
#' @param seed the random seed. If \samp{NULL} the seed is unchanged.
#' @return An \samp{mcnode} object.
#' @seealso \code{\link{mcdata}}, \code{\link{mcstoc}}, \code{\link{switch}}.
#' @keywords methods
#' @examples
#' conc1 <- mcstoc(rnorm, type="VU", mean=10, sd=2)
#' conc2 <- mcstoc(runif, type="VU", min=-6, max=-5)
#' conc3 <- mcdata(0, type="VU")
#' ## Randomly in the cells
#' whichdist <- mcstoc(rempiricalD, type="VU", values=1:3, prob=c(.75, .20, .05))
#' mcprobtree(whichdist, list("1"=conc1, "2"=conc2, "3"=conc3), type="VU")
#' ## Equivalent using weights
#' mcprobtree(c(.75, .20, .05), list("1"=conc1, "2"=conc2, "3"=conc3), type="VU")
#' @export
mcprobtree <- function(mcswitch, mcvalues, type=c("V","U","VU","0"), nsv=ndvar(), nsu=ndunc(), nvariates=1, outm="each", seed=NULL)
{
    if (!is.null(seed)) set.seed(seed)
    if (!is.character(outm) || (outm != "none" && outm != "each" &&
        !all(sapply(outm, exists, mode = "function"))))
        stop("outm should be 'none','each' or the name a valid function")
    type <- match.arg(type)
    nsv
    nsu
    if (type == "V") nsu <- 1 else if (type == "U") nsv <- 1 else if(type=="0") {nsv <- nsu <- 1}
    stoc <- as.list(substitute(mcvalues))
    choixstoc <- as.numeric(names(stoc))
    if (any(is.na(choixstoc[-1])))
        stop("Names of the mcvalues element should be convertible as numeric")
    stoc <- substitute(mcvalues)
    if (inherits(mcswitch, "mcnode")) {
        typem <- attr(mcswitch, "type")
        if ((type == "V" && typem == "U") || (type == "U" &&
            typem == "V") || (type != "VU" && typem == "VU") ||
            (type == "0" && typem != "0"))
            stop("Incompatible type and type of mcswitch")
        mcswitch <- mcdata(mcswitch, type = type, nsv = nsv, nsu = nsu, nvariates = nvariates)
    }
    else if (is.vector(mcswitch)) {
        if (length(mcswitch) != length(mcvalues))
            stop("the vector mcswitch should have the same length as mcvalues")
        if (any(mcswitch < 0) || sum(mcswitch) == 0 || any(!is.finite(mcswitch)))
            stop("mcswitch values should be finite, nonnegative and not all zero")
        mcswitch <- mcswitch/sum(mcswitch)
        dimf <- prod(c(nsv, nsu, nvariates))
        mcswitch <- sample(x = choixstoc[-1], size = dimf, replace = TRUE, prob = mcswitch)
    }
    else stop("mcswitch should be an mcnode or a vector")
    choixswitch <- unique(mcswitch)
    if (!all(choixswitch %in% choixstoc))
        stop("Some values of mcswitch are not a name of elements of mcvalues")
    res <- mcdata(NA, type = type, nsv = nsv, nsu = nsu, nvariates = nvariates)
    for (i in choixswitch) {
        whichswitch <- mcswitch == i
        nbcall <- which(choixstoc == i)
        thecall <- eval(stoc[[nbcall]])
        if (!inherits(thecall, "mcnode"))
            stop("One mcvalues does not lead to a mcnode")
        typen <- attr(thecall, "type")
        if ((type == "V" && typen == "U") || (type == "U" &&
            typen == "V") || (type != "VU" && typen == "VU") ||
            (type == "0" && typen != "0"))
            stop("One element of mcvalues leads to an mcnode of incorrect type")
        thecall <- mcdata(thecall, type = type, nsv = nsv, nsu = nsu, nvariates = nvariates)
        res[whichswitch] <- thecall[whichswitch]
    }
    class(res) <- "mcnode"
    attr(res, "type") <- type
    attr(res, "outm") <- outm
    return(res)
}
