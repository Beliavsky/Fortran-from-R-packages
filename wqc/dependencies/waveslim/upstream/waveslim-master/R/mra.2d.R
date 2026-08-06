#' Multiresolution Analysis of an Image
#' 
#' This function performs a level \eqn{J} additive decomposition of the input
#' matrix or image using the pyramid algorithm (Mallat 1989).
#' 
#' This code implements a two-dimensional multiresolution analysis by
#' performing the one-dimensional pyramid algorithm (Mallat 1989) on the rows
#' and columns of the input matrix.  Either the DWT or MODWT may be used to
#' compute the multiresolution analysis, which is an additive decomposition of
#' the original matrix (image).
#' 
#' @param x A matrix or image containing the data be to decomposed.  This must
#' be have dyadic length in both dimensions (but not necessarily the same) for
#' \code{method="dwt"}.
#' @param wf Name of the wavelet filter to use in the decomposition.  By
#' default this is set to \code{"la8"}, the Daubechies orthonormal compactly
#' supported wavelet of length L=8 least asymmetric family.
#' @param J Specifies the depth of the decomposition.  This must be a number
#' less than or equal to log(length(x),2).
#' @param method Either \code{"dwt"} or \code{"modwt"}.
#' @param boundary Character string specifying the boundary condition.  If
#' \code{boundary=="periodic"} the default, then the matrix you decompose is
#' assumed to be periodic on its defined interval,\cr if
#' \code{boundary=="reflection"}, the matrix beyond its boundaries is assumed
#' to be a symmetric reflection of itself.
#' @return Basically, a list with the following components \item{LH?}{Wavelet
#' detail image in the horizontal direction.} \item{HL?}{Wavelet detail image
#' in the vertical direction.} \item{HH?}{Wavelet detail image in the diagonal
#' direction.} \item{LLJ}{Wavelet smooth image at the coarsest resolution.}
#' \item{J}{Depth of the wavelet transform.} \item{wavelet}{Name of the wavelet
#' filter used.} \item{boundary}{How the boundaries were handled.}
#' @author B. Whitcher
#' @seealso \code{\link{dwt.2d}}, \code{\link{modwt.2d}}
#' @references Mallat, S. G. (1989) A theory for multiresolution signal
#' decomposition: the wavelet representation, \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, \bold{11}, No. 7, 674-693.
#' 
#' Mallat, S. G. (1998) \emph{A Wavelet Tour of Signal Processing}, Academic
#' Press.
#' @keywords ts
#' @examples
#' 
#' ## Easy check to see if it works...
#' ## --------------------------------
#' 
#' x <- matrix(rnorm(32*32), 32, 32)
#' # MODWT
#' x.mra <- mra.2d(x, method="modwt")
#' x.mra.sum <- x.mra[[1]]
#' for(j in 2:length(x.mra))
#'   x.mra.sum <- x.mra.sum + x.mra[[j]]
#' sum((x - x.mra.sum)^2)
#' 
#' # DWT
#' x.mra <- mra.2d(x, method="dwt")
#' x.mra.sum <- x.mra[[1]]
#' for(j in 2:length(x.mra))
#'   x.mra.sum <- x.mra.sum + x.mra[[j]]
#' sum((x - x.mra.sum)^2)
#' 
#' @export mra.2d
mra.2d <-
  function(x, wf="la8", J=4, method="modwt", boundary="periodic")
{
  m <- dim(x)[1]
  n <- dim(x)[2]

  switch(boundary,
    "periodic" = invisible(),
    stop("Invalid boundary rule in mra"))

  if(method == "modwt") {
    x.wt <- modwt.2d(x, wf, J, "periodic")
  } else {
    x.wt <- dwt.2d(x, wf, J, "periodic")
  }
  
  x.mra <- vector("list", 3*J+1)

  ## Smooth
  zero <- vector("list", 3*J+1)
  names(zero) <-
    c(matrix(rbind(paste("LH", 1:J, sep=""), paste("HL", 1:J, sep=""),
                   paste("HH", 1:J, sep="")), nrow=1), paste("LL", J, sep=""))
  attr(zero, "J") <- J
  attr(zero, "wavelet") <- wf
  attr(zero, "boundary") <- boundary
  zero[[3*J+1]] <- x.wt[[3*J+1]]
  if(method == "modwt") {
    for(k in 1:(3*J))
      zero[[k]] <- matrix(0, m, n)
    x.mra[[3*J+1]] <- imodwt.2d(zero)
  } else {
    for(k in 1:J)
      zero[[3*(k-1)+1]] <- zero[[3*(k-1)+2]] <- zero[[3*k]] <-
        matrix(0, m/2^k, n/2^k)
    x.mra[[3*J+1]] <- idwt.2d(zero)
  }

## Details
  for(j in (3*J):1) {
    Jj <- ceiling(j/3)
    zero <- vector("list", 3*Jj+1)
    names(zero) <-
      c(matrix(rbind(paste("LH", 1:Jj, sep=""), paste("HL", 1:Jj, sep=""),
                     paste("HH", 1:Jj, sep="")), nrow=1),
        paste("LL", Jj, sep=""))
    attr(zero, "J") <- Jj
    attr(zero, "wavelet") <- wf
    attr(zero, "boundary") <- boundary
    zero[[j]] <- x.wt[[j]]
    if(method == "modwt") {
      for(k in names(zero)[-charmatch(names(zero)[j], names(zero))])
        zero[[k]] <- matrix(0, m, n)
      x.mra[[j]] <- imodwt.2d(zero)
    } else {
      for(k in 1:Jj)
        zero[[3*(k-1)+1]] <- zero[[3*(k-1)+2]] <- zero[[3*k]] <-
          matrix(0, m/2^k, n/2^k)
      zero[[3*Jj+1]] <- matrix(0, m/2^Jj, n/2^Jj)
      zero[[j]] <- x.wt[[j]]
      x.mra[[j]] <- idwt.2d(zero)
    }
  }

  names(x.mra) <-
    c(matrix(rbind(paste("LH", 1:J, sep=""), paste("HL", 1:J, sep=""),
                   paste("HH", 1:J, sep="")), nrow=1), paste("LL", Jj, sep=""))
   return(x.mra)
}

