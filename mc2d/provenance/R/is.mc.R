#' Dimension of mcnode and mc Objects
#'
#' Provides the dimension of a \samp{mcnode} or a \samp{mc} object: the number of simulations
#' in the variability dimension, in the uncertainty dimension, and the maximum number of
#' variates.
#'
#' @note These functions do not test whether the object is correctly built. See
#'   \code{\link{is.mcnode}} and \code{\link{is.mc}}.
#'
#' @param x a \samp{mcnode} or a \samp{mc} object.
#' @param index if \samp{TRUE}, give the index of the type rather than the type (for
#'   \samp{typemcnode} only).
#' @return \samp{dimmcnode}: a vector of three scalars (nsv, nsu, nvariates).
#'   \samp{dimmc}: a vector of three scalars (max nsv, max nsu, max variates).
#'   \samp{typemcnode}: \samp{"0"}, \samp{"V"}, \samp{"U"} or \samp{"VU"}, or the
#'   corresponding index if \samp{index=TRUE}. \samp{NULL} if none found.
#'   \samp{is.mc}: \samp{TRUE} or \samp{FALSE}.
#'   \samp{is.mcnode}: \samp{TRUE} or \samp{FALSE}.
#' @seealso \code{\link{mcnode}}, \code{\link{mc}}.
#' @keywords utilities
#' @name dimmcnode
#' @examples
#' data(total)
#' dimmcnode(xVUM2)
#' dimmc(total)
#' typemcnode(total$xVUM2)
#' is.mcnode(xVU)
#' is.mc(total)
#' @export
dimmcnode <- function(x)
{
  if(!inherits(x,"mcnode")) stop("x is not an mcnode object")
  y <- dim(x)
  names(y) <- c("nsv","nsu","nvariates")
  return(y)}

#' @rdname dimmcnode
#' @export
dimmc <- function(x)
{
  if(!inherits(x,"mc")) stop("x is not an mc object")
  lesdim <- sapply(x,dimmcnode)
  y <- apply(lesdim,1,max)
  names(y) <- c("nsv","nsu","max variates")
  return(y)}

#' @rdname dimmcnode
#' @export
typemcnode <- function(x,index=FALSE)
{
  if(!inherits(x,"mcnode")) stop("x is not an mcnode object")
  type <- attr(x,"type")
  if(index) return(which(c("0", "V","U","VU")==type)) else return(type)
}

#' @rdname dimmcnode
#' @export
is.mc <- function(x)
{
  if (!inherits(x, "mc")) return(FALSE)
  x <- unclass(x)
  if(!is.list(x)) return(FALSE)
  mcn <- sapply(x,is.mcnode)
  if(!all(mcn)) return(FALSE)
  nsim <- sapply(x, dim)
  if(!all(nsim[1,] %in% c(1,max(nsim[1,])))) return(FALSE)
  if(!all(nsim[2,] %in% c(1,max(nsim[2,])))) return(FALSE)
  return(TRUE)}

#' @rdname dimmcnode
#' @export
is.mcnode <- function(x)
{
  if(!inherits(x,"mcnode")) return(FALSE)
  type <- typemcnode(x)
  if(is.null(type)) return(FALSE)
  x <- unclass(x)
  if(!is.numeric(x) && !is.logical(x)) return(FALSE)
  dimx <- dim(x)
  if(type == "0" && (!is.array(x) || dimx[1]!=1 && dimx[2]!=1)) return(FALSE)
  if(type == "V" && (!is.array(x) || dimx[2]!=1)) return(FALSE)
  if(type == "U" && (!is.array(x) || dimx[1]!=1)) return(FALSE)
  if(type == "VU" && !is.array(x)) return(FALSE)
  return (TRUE)}
