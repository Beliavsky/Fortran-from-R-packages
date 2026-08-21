#' Build mcnode Objects from Data
#'
#' Creates a \samp{mcnode} object from a vector, an array or another \samp{mcnode}.
#'
#' A \samp{mcnode} object is the basic element of a \code{\link{mc}} object.
#' It is an array of dimension \samp{(nsv x nsu x nvariates)}. Four types exist:
#' \describe{
#'   \item{\samp{"V" mcnode}}{for "Variability", arrays of dimension \samp{(nsv x 1 x nvariates)}.
#'     The variability of the data should reflect parameter variability.}
#'   \item{\samp{"U" mcnode}}{for "Uncertainty", arrays of dimension \samp{(1 x nsu x nvariates)}.
#'     The variability reflects parameter uncertainty.}
#'   \item{\samp{"VU" mcnode}}{for "Variability and Uncertainty", arrays of dimension
#'     \samp{(nsv x nsu x nvariates)}.}
#'   \item{\samp{"0" mcnode}}{for "Neither Variability nor Uncertainty", arrays of dimension
#'     \samp{(1 x 1 x nvariates)}.}
#' }
#'
#' Multivariate nodes (\samp{nvariates != 1}) should be used for multivariate distributions
#' (\code{\link{rmultinomial}}, \code{\link{rmultinormal}}, \code{\link{rempiricalD}},
#' \code{\link{rdirichlet}}).
#'
#' Recycling rules are limited: recycling is only permitted to fill a dimension from 1 to
#' the final size of the dimension.
#'
#' @param data the numeric/logical vector/matrix/array of data or a \samp{mcnode} object.
#' @param type the type of node to be built. By default, a \samp{"V"} node.
#' @param nsv the variability dimension (\samp{type="V"} or \samp{type="VU"}) of the node.
#'   By default: the current value in \code{\link{mc.control}}.
#' @param nsu the uncertainty dimension (\samp{type="U"} or \samp{type="VU"}) of the node.
#'   By default: the current value in \code{\link{mc.control}}.
#' @param nvariates the number of variates. By default: 1.
#' @param outm the output of the \samp{mcnode} for multivariate nodes. May be \samp{"each"}
#'   (default), \samp{"none"}, or a vector of function names applied on the variates dimension
#'   before any output.
#' @return An \samp{mcnode} object.
#' @seealso \code{\link{mcstoc}} to build a stochastic \samp{mcnode} object.
#'   \code{\link{mc}} to build a Monte-Carlo object.
#'   \code{\link{is.mcnode}}, \code{\link{dimmcnode}}, \code{\link{typemcnode}}.
#' @keywords methods
#' @name mcnode
#' @examples
#' oldvar <- ndvar()
#' oldunc <- ndunc()
#' ndvar(3)
#' ndunc(5)
#'
#' (x0 <- mcdata(100, type = "0"))
#' (xV <- mcdata(1:ndvar(), type = "V"))
#' (xU <- mcdata(10 * 1:ndunc(), type = "U"))
#' (xVU <- mcdata(1:(ndvar() * ndunc()), type = "VU"))
#'
#' ## Multivariates
#' (x0M <- mcdata(1:2, type = "0", nvariates = 2))
#' (xVM <- mcdata(1:(2 * ndvar()), type = "V", nvariates = 2))
#'
#' ndvar(oldvar)
#' ndunc(oldunc)
#' @export
mcdata <- function(data, type=c("V","U","VU","0"), nsv=ndvar(), nsu=ndunc(),nvariates=1,outm="each")
{
  type <- match.arg(type)

    if(!is.character(outm)  || (outm != "none" && outm != "each" && !all(sapply(outm,exists,mode="function"))))
      stop("outm should be 'none','each' or the name a valid function")

  dimf <- switch(type,"V"=c(nsv,1,nvariates),"U"=c(1,nsu,nvariates),"VU"=c(nsv,nsu,nvariates),"0"=c(1,1,nvariates))

  if(inherits(data,"mcnode") ){             # data = mcnode
    typem <- attr(data,"type")
    dimd <- dim(data)
    err <- paste("The output node dimension is not compatible with the input node dimension. Should be of dim: ",paste(dimf,collapse=" "))
    if(dimd[3] != 1 && dimd[3]!=dimf[3]) stop(err)
    if(typem=="VU") {
      if(dimd[1]!=dimf[1] || dimd[2]!=dimf[2]) stop(err)
      if(type!="VU") stop("Incompatible node type. A 'VU' mcnode makes only a 'VU'mcnode.")
      data <- array(data,dim=dimf)}                                             # recycle on nvariates
    else if(typem=="U") {
      if(dimd[2]!=dimf[2]) stop(err)
      if(type!="VU" && type!="U") stop("Incompatible node type. A 'U' mcnode makes a 'U' or a 'VU' mcnode.")
      data <- array(apply(data,3,rep,each=dimf[1]),dim=dimf)}                   # recycle on nsu and nvariates
    else if(typem=="V") {
      if(dimd[1]!=dimf[1]) stop(err)
      if(type!="VU" && type!="V") stop("Incompatible node type. A 'V' mcnode makes a 'V' or a 'VU' mcnode.")
      data <- array(apply(data,3,rep,times=dimf[2]),dim=dimf)}                  # recycle on nsv and nvariates
    else if(typem=="0") {
      data <- array(apply(data,3,rep,each=dimf[1]*dimf[2]),dim=dimf)}           # recycle on nvariates
  }
  else {                                    # data non mcnode
      if(!is.numeric(data) && !is.logical(data)) stop("data should be numeric or logical")
      if(is.vector(data)){
        l <- length(data)
        if(l != 1 && l != dimf[1]*dimf[2] && l!= prod(dimf)) stop("The vector size is not compatible
          with the node dimension. Length should be 1 or n=",dimf[1]*dimf[2]," or n=",prod(dimf))
        data <- array(data,dim=dimf)}

      else if(is.array(data)){
        dimd <- dim(data)
        l <- length(dimd)
        if(l > 3) stop("Maximum accepted dim of arrays is 3.")
        if(l == 2) dimd <- c(dimd,0)  # just to simplify the following tests
        if( (type=="VU" && dimd[1]== nsv && dimd[2]==nsu && dimd[3] %in% c(0,1,nvariates)) ||
            (type=="V" &&  dimd[1]== nsv && dimd[2]==1 && dimd[3] %in% c(0,1,nvariates)) ||
            (type=="U" &&  dimd[1]== 1 && dimd[2]==nsu && dimd[3] %in% c(0,1,nvariates)) ||
            (type=="0" &&  dimd[1]== 1 && dimd[1]==1 && dimd[3] %in% c(0,1,nvariates)) )
        data <- array(data,dim=dimf)
        else stop("The array size is not compatible with the node dimension. Should be of dim: ",paste(dimf,collapse=" "))
        }

      else stop("data should be a vector, a matrix, an array or a mcnode")}

  class(data) <- "mcnode"
  attr(data,"type") <- type
  attr(data,"outm") <- outm
  return(data)
}
