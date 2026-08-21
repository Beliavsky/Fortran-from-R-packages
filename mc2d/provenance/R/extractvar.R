#' Utilities for Multivariate Nodes
#'
#' \samp{extractvar} extracts one or more variates from a multivariate \samp{mcnode}.
#' \samp{addvar} combines consistent \samp{mcnode}s into a single multivariate \samp{mcnode}.
#'
#' The \samp{outm} attribute of the output of \samp{addvar} is taken from the first element.
#'
#' @param x a multivariate \samp{mcnode}.
#' @param which a vector specifying which variate(s) to extract.
#' @param ... \samp{mcnode}s to be combined in a multivariate \samp{mcnode}.
#'   These \samp{mcnode}s should be of the same type and dimension.
#' @return The new \samp{mcnode}.
#' @seealso \code{\link{mcnode}} for \samp{mcnode} objects.
#' @keywords methods
#' @examples
#' x <- mcdata(0:3, "0", nvariates = 4)
#' y <- extractvar(x, c(1, 3))
#' y
#' addvar(x, y)
#' @export
extractvar <- function(x, which = 1)
{
  if(missing(x) || !inherits(x,"mcnode")) stop("extractvar need a mcnode object")
  dimm <- dim(x)
  if(any(which > dimm[3]) || any(which < 1)) stop("Incorrect value of which")
  x <- mcdata(x[,,which,drop=FALSE],type=typemcnode(x),nsv=dimm[1],nsu=dimm[2],nvariates=length(which))
  return(x)
}

#' @rdname extractvar
#' @export
addvar <- function(...)
{
  argsd <- list(...)
  ismc <-  sapply(argsd,inherits,"mcnode")
  if(!all(ismc)) stop("addvar needs mcnodes object")
  dimm <- sapply(argsd,dim)
  if(any(dimm[1,] != dimm[1,1]) || any(dimm[2,]!=dimm[2,1]))
    stop("Arguments should have the same dimension of variability and uncertainty")
  typem <- sapply(argsd,attr,which="type")
  if(any(typem != typem[1])) stop("Arguments should be of same type")
  outm <- attr(argsd[[1]],which="outm")
  argsd <- unlist(argsd)
  dim(argsd) <- c(dimm[1:2,1],sum(dimm[3,]))
  attr(argsd,which="type") <- typem[1]
  attr(argsd,which="outm") <- outm
  class(argsd) <- "mcnode"
  return(argsd)
}
