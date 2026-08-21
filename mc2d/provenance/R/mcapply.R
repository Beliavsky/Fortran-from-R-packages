#' Apply Functions Over mc or mcnode Objects
#'
#' Apply a function on all values or over a given dimension of an \samp{mcnode} object.
#' May be used for all \samp{mcnode} objects of an \samp{mc} object.
#'
#' @param x a \samp{mc} or a \samp{mcnode} object.
#' @param margin the dimension on which to apply the function. May be \samp{"all"} (default)
#'   to apply the function on all values, \samp{"var"} to apply on the variability dimension,
#'   \samp{"unc"} to apply on the uncertainty dimension, or \samp{"variates"} to apply on
#'   the variates. Note: do not use \samp{"var"} for \samp{"variates"}.
#' @param fun the function to be applied. When applied to a vector of length \samp{n},
#'   \samp{fun} should return a vector of length \samp{n} or \samp{1}.
#' @param ... optional arguments to \samp{fun}.
#' @return If \samp{fun} returns a vector of length \samp{n} or if \samp{margin="all"},
#'   the returned \samp{mcnode}s have the same type and dimension as \samp{x}.
#'   Otherwise, the type of \samp{mcnode} is changed accordingly.
#' @seealso \code{\link{apply}}, \code{\link{mc}}, \code{\link{mcnode}}.
#' @keywords misc
#' @examples
#' data(total)
#' xVUM
#' mcapply(xVUM, "unc", sum)
#' mcapply(xVUM, "var", sum)
#' mcapply(xVUM, "all", sum)
#' mcapply(xVUM, "variates", sum)
#' mcapply(total, "all", exp)
#' @export
mcapply <- function(x, margin=c("all","var","unc","variates"), fun, ...)
{
  margin <- match.arg(margin)

  if(is.mc(x)) {
    y <- lapply(x, mcapply, margin=margin, fun=fun, ...)
    class(y) <- "mc"
    return(y)}

  if(!is.mcnode(x)) stop("x should be an mcnode")

  typen <- attr(x,"type")
  oldoutm <- attr(x,which="outm")
  dimn <-  dim(x)

  dime <- switch(margin, "unc" = c(1,3), "var" = c(2,3), "all" = c(1,2,3), "variates" = c(1,2))
  nlo <-  switch(margin, "unc" = prod(dimn[c(1,3)]), "var" = prod(dimn[c(2,3)]), "all" = prod(dimn), "variates" = prod(dimn[c(1,2)]))

  dat <- apply(x, dime, fun, ...)
  nl <- length(dat)

  if(nl == length(x)){		 # fun(n) give n scalar
	if(margin == "unc") 			dat <- aperm(dat,c(2,1,3))
	else if(margin == "variates") 	dat <- aperm(dat,c(2,3,1))
	x[] <- dat
	return(x)
	}

	if(nl != nlo) stop("fun(x) with length(x) = n should return a vector of length n or 1")

  nsu <- ifelse(margin=="unc", 1, dimn[2])
  nsv <- ifelse(margin=="var", 1, dimn[1])
  nva <- ifelse(margin=="variates", 1, dimn[3])

  if(typen=="0") ntype <- "0"

  else if(typen=="V"){
    if(margin=="var") ntype <- "0"
    else ntype <- "V"}

  else if(typen=="U"){
    if(margin=="unc") ntype <- "0"
    else ntype <- "U"}

  else if(margin=="unc") ntype <- "V"
  else if(margin=="var") ntype <- "U"
  else ntype <- "VU"

 return(mcdata(as.vector(dat),type=ntype,nsv=nsv,nsu=nsu,nvariates=nva,outm=oldoutm))
}
