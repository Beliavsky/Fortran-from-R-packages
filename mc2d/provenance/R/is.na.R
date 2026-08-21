#' Finite, Infinite, NA and NaN Numbers in mcnode
#'
#' \samp{is.na}, \samp{is.nan}, \samp{is.finite} and \samp{is.infinite} return a logical
#' \samp{mcnode} of the same dimension as \samp{x}.
#'
#' @param x a \samp{mcnode} object.
#' @return a logical \samp{mcnode} object.
#' @seealso \code{\link{is.finite}}, \code{\link{NA}}
#' @keywords NA
#' @name NA.mcnode
#' @examples
#' x <- log(mcstoc(rnorm, nsv=1001))
#' x
#' is.na(x)
#' @export
is.na.mcnode <- function(x)
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}

#' @rdname NA.mcnode
#' @export
is.nan.mcnode <- function(x)
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}

#' @rdname NA.mcnode
#' @export
is.finite.mcnode <- function(x)
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}

#' @rdname NA.mcnode
#' @export
is.infinite.mcnode <- function(x)
{
  y <- NextMethod()
  attr(y,"type") <- attr(x,"type")
  attr(y,"outm") <- attr(x,"outm")
  class(y) <- "mcnode"
  return(y)}
