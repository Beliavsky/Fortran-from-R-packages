#' Unclass the mc or mcnode Object
#'
#' Unclasses an \samp{mc} object into a list of arrays, or an \samp{mcnode} object into
#' an array.
#'
#' @param x a \samp{mc} or a \samp{mcnode} object.
#' @param drop should dimensions of size 1 be dropped (see \code{\link{drop}}).
#' @return If \samp{x} is an \samp{mc} object: a list of arrays. If \samp{drop=TRUE},
#'   a list of vectors, matrices and arrays.
#'   If \samp{x} is an \samp{mcnode} object: an array. If \samp{drop=TRUE},
#'   a vector, matrix or array.
#' @keywords manip
#' @examples
#' data(total)
#' ## A vector
#' unmc(total$xV, drop=TRUE)
#' ## An array
#' unmc(total$xV, drop=FALSE)
#' @export
unmc <- function(x, drop=TRUE)
{
  unmcnode <- function(y){
    attr(y,"type") <- NULL
    attr(y,"outm") <- NULL
    y <- unclass(y)
    if(drop) y <- drop(y)
    return(y)}

  if(is.mc(x)){
    x <- lapply(x,unmcnode)
    x <- unclass(x)
    return(x)}

  return(unmcnode(x))
}
