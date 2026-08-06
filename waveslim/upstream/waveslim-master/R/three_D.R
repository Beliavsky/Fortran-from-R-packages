###########################################################################
###########################################################################
###########################################################################



#' Three Dimensional Separable Discrete Wavelet Transform
#' 
#' Three-dimensional separable discrete wavelet transform (DWT).
#' 
#' 
#' @usage dwt.3d(x, wf, J = 4, boundary = "periodic")
#' @usage idwt.3d(y)
#' @aliases dwt.3d idwt.3d
#' @param x input array
#' @param wf name of the wavelet filter to use in the decomposition
#' @param J depth of the decomposition, must be a number less than or equal to
#' log(minZ,Y,Z,2)
#' @param boundary only \code{"periodic"} is currently implemented
#' @param y an object of class \code{dwt.3d}
#' @author B. Whitcher
#' @keywords ts
#' @export dwt.3d
dwt.3d <- function(x, wf, J=4, boundary="periodic")
{
  nx <- dim(x)[1]
  storage.mode(nx) <- "integer"
  ny <- dim(x)[2]
  storage.mode(ny) <- "integer"
  nz <- dim(x)[3]
  storage.mode(nz) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  z <- array(0, dim=c(nx,ny,nz)/2)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 7*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C(C_three_D_dwt, "cube"=as.double(x), "NX"=nx, "NY"=ny,
              "NZ"=nz, "filter.length"=L, "hpf"=h, "lpf"=g, "LLL"=z,
              "HLL"=z, "LHL"=z, "LLH"=z, "HHL"=z, "HLH"=z, "LHH"=z,
              "HHH"=z)[8:15]
    if(j < J) {
      index <- (7*(j-1)+1):(7*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      nx <- dim(x)[1]
      storage.mode(nx) <- "integer"
      ny <- dim(x)[2]
      storage.mode(ny) <- "integer"
      nz <- dim(x)[3]
      storage.mode(nz) <- "integer"
      z <- array(0, dim=c(nx,ny,nz)/2)
      storage.mode(z) <- "double"
    }
    else {
      index <- (7*(j-1)+1):(7*j+1)
      x.wt[index] <- out[c(2:8,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:8,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  class(x.wt) <- "dwt.3d"
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  return(x.wt)
}

###########################################################################
###########################################################################
###########################################################################

idwt.3d <- function(y)
{
  J <- attributes(y)$J
  LLL <- paste("LLL", J, sep="")

  wf <- attributes(y)$wavelet
  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  y.in <- y$LLL

  for(j in J:1) {
    HLL <- paste("HLL", j, sep="")
    LHL <- paste("LHL", j, sep="")
    LLH <- paste("LLH", j, sep="")
    HHL <- paste("HHL", j, sep="")
    HLH <- paste("HLH", j, sep="")
    LHH <- paste("LHH", j, sep="")
    HHH <- paste("HHH", j, sep="")
    
    nx <- dim(y.in)[1]
    storage.mode(nx) <- "integer"
    ny <- dim(y.in)[2]
    storage.mode(ny) <- "integer"
    nz <- dim(y.in)[3]
    storage.mode(nz) <- "integer"

    z <- array(0, dim=2*c(nx, ny, nz))
    storage.mode(z) <- "double"

    out <- .C(C_three_D_idwt, as.double(y.in), as.double(y[[HLL]]),
              as.double(y[[LHL]]), as.double(y[[LLH]]),
              as.double(y[[HHL]]), as.double(y[[HLH]]),
              as.double(y[[LHH]]), as.double(y[[HHH]]), 
              nx, ny, nz, L, h, g, "Y"=z)

    y.in <- out$Y
  }
  zapsmall(y.in)
}



#' Three Dimensional Separable Maximal Ovelrap Discrete Wavelet Transform
#' 
#' Three-dimensional separable maximal overlap discrete wavelet transform
#' (MODWT).
#' 
#' 
#' @usage modwt.3d(x, wf, J = 4, boundary = "periodic")
#' @usage imodwt.3d(y)
#' @aliases modwt.3d imodwt.3d
#' @param x input array
#' @param wf name of the wavelet filter to use in the decomposition
#' @param J depth of the decomposition
#' @param boundary only \code{"periodic"} is currently implemented
#' @param y an object of class \code{modwt.3d}
#' @author B. Whitcher
#' @keywords ts
#' @export modwt.3d 
#' @export imodwt.3d
modwt.3d <- function(x, wf, J=4, boundary="periodic")
{
  nx <- dim(x)[1]
  storage.mode(nx) <- "integer"
  ny <- dim(x)[2]
  storage.mode(ny) <- "integer"
  nz <- dim(x)[3]
  storage.mode(nz) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  z <- array(0, dim=c(nx,ny,nz))
  storage.mode(z) <- "double"

  x.wt <- vector("list", 7*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C(C_three_D_modwt, "cube"=as.double(x), "NX"=nx, "NY"=ny,
              "NZ"=nz, "J"=j, "filter.length"=L, "hpf"=h, "lpf"=g,
              "LLL"=z, "HLL"=z, "LHL"=z, "LLH"=z, "HHL"=z, "HLH"=z,
              "LHH"=z, "HHH"=z)[9:16]
    if(j < J) {
      index <- (7*(j-1)+1):(7*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      nx <- dim(x)[1]
      storage.mode(nx) <- "integer"
      ny <- dim(x)[2]
      storage.mode(ny) <- "integer"
      nz <- dim(x)[3]
      storage.mode(nz) <- "integer"
      z <- array(0, dim=c(nx,ny,nz))
      storage.mode(z) <- "double"
    }
    else {
      index <- (7*(j-1)+1):(7*j+1)
      x.wt[index] <- out[c(2:8,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:8,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  class(x.wt) <- "modwt.3d"
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  return(x.wt)
}

###########################################################################
###########################################################################
###########################################################################

imodwt.3d <- function(y)
{
  J <- attributes(y)$J
  LLL <- paste("LLL", J, sep="")

  wf <- attributes(y)$wavelet
  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  y.in <- y$LLL

  for(j in J:1) {
    HLL <- paste("HLL", j, sep="")
    LHL <- paste("LHL", j, sep="")
    LLH <- paste("LLH", j, sep="")
    HHL <- paste("HHL", j, sep="")
    HLH <- paste("HLH", j, sep="")
    LHH <- paste("LHH", j, sep="")
    HHH <- paste("HHH", j, sep="")
    
    nx <- dim(y.in)[1]
    storage.mode(nx) <- "integer"
    ny <- dim(y.in)[2]
    storage.mode(ny) <- "integer"
    nz <- dim(y.in)[3]
    storage.mode(nz) <- "integer"

    z <- array(0, dim=c(nx, ny, nz))
    storage.mode(z) <- "double"

    out <- .C(C_three_D_imodwt, as.double(y.in), as.double(y[[HLL]]),
              as.double(y[[LHL]]), as.double(y[[LLH]]),
              as.double(y[[HHL]]), as.double(y[[HLH]]),
              as.double(y[[LHH]]), as.double(y[[HHH]]), 
              nx, ny, nz, j, L, h, g, "Y"=z)

    y.in <- out$Y
  }
  zapsmall(y.in)
}

###########################################################################
###########################################################################
###########################################################################



#' Three Dimensional Multiresolution Analysis
#' 
#' This function performs a level \eqn{J} additive decomposition of the input
#' array using the pyramid algorithm (Mallat 1989).
#' 
#' This code implements a three-dimensional multiresolution analysis by
#' performing the one-dimensional pyramid algorithm (Mallat 1989) on each
#' dimension of the input array.  Either the DWT or MODWT may be used to
#' compute the multiresolution analysis, which is an additive decomposition of
#' the original array.
#' 
#' @param x A three-dimensional array containing the data be to decomposed.
#' This must be have dyadic length in all three dimensions (but not necessarily
#' the same) for \code{method="dwt"}.
#' @param wf Name of the wavelet filter to use in the decomposition.  By
#' default this is set to \code{"la8"}, the Daubechies orthonormal compactly
#' supported wavelet of length \eqn{L=8} least asymmetric family.
#' @param J Specifies the depth of the decomposition.  This must be a number
#' less than or equal to \eqn{\log(\mbox{length}(x),2)}{log(length(x),2)}.
#' @param method Either \code{"dwt"} or \code{"modwt"}.
#' @param boundary Character string specifying the boundary condition.  If
#' \code{boundary=="periodic"} the default and only method implemented, then
#' the matrix you decompose is assumed to be periodic on its defined interval.
#' @return List structure containing the filter triplets associated with the
#' multiresolution analysis.
#' @author B. Whitcher
#' @seealso \code{\link{dwt.3d}}, \code{\link{modwt.3d}}
#' @references Mallat, S. G. (1989) A theory for multiresolution signal
#' decomposition: the wavelet representation, \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, \bold{11}, No. 7, 674-693.
#' 
#' Mallat, S. G. (1998) \emph{A Wavelet Tour of Signal Processing}, Academic
#' Press.
#' @keywords ts
#' @export mra.3d
mra.3d <- function(x, wf="la8", J=4, method="modwt", boundary="periodic")
{
  nx <- dim(x)[1]
  ny <- dim(x)[2]
  nz <- dim(x)[3]

  if(method == "modwt") {
    x.wt <- modwt.3d(x, wf, J, "periodic")
  } else {
    x.wt <- dwt.3d(x, wf, J, "periodic")
  }
  
  x.mra <- vector("list", 7*J+1)
  names(x.mra) <-
    c(matrix(rbind(paste("HLL", 1:J, sep=""), paste("LHL", 1:J, sep=""),
                   paste("LLH", 1:J, sep=""), paste("HHL", 1:J, sep=""),
                   paste("HLH", 1:J, sep=""), paste("LHH", 1:J, sep=""),
                   paste("HHH", 1:J, sep="")), nrow=1),
      paste("LLL", J, sep=""))

  ## Smooth
  zero <- vector("list", 7*J+1)
  names(zero) <- names(x.mra)
  attr(zero, "J") <- J
  attr(zero, "wavelet") <- wf
  attr(zero, "boundary") <- boundary
  zero[[7*J+1]] <- x.wt[[7*J+1]]
  if(method == "modwt") {
    class(x.wt) <- "modwt.3d"
    for(k in 1:(7*J))
      zero[[k]] <- array(0, dim=c(nx,ny,nz))
    x.mra[[7*J+1]] <- imodwt.3d(zero)
  } else {
    class(x.wt) <- "dwt.3d"
    for(k in 1:J)
      zero[[7*(k-1)+1]] <- zero[[7*(k-1)+2]] <- zero[[7*(k-1)+3]] <-
        zero[[7*(k-1)+4]] <- zero[[7*(k-1)+5]] <- zero[[7*(k-1)+6]] <-
          zero[[7*k]] <- array(0, dim=c(nx,ny,nz)/2^k)
    x.mra[[7*J+1]] <- idwt.3d(zero)
  }

  ## Details
  for(j in (7*J):1) {
    Jj <- ceiling(j/7)
    zero <- vector("list", 7*Jj+1)
    names(zero) <-
      c(matrix(rbind(paste("HLL", 1:Jj, sep=""), paste("LHL", 1:Jj, sep=""),
                     paste("LLH", 1:Jj, sep=""), paste("HHL", 1:Jj, sep=""),
                     paste("HLH", 1:Jj, sep=""), paste("LHH", 1:Jj, sep=""),
                     paste("HHH", 1:Jj, sep="")), nrow=1),
        paste("LLL", Jj, sep=""))
    attr(zero, "J") <- Jj
    attr(zero, "wavelet") <- wf
    attr(zero, "boundary") <- boundary
    zero[[j]] <- x.wt[[j]]
    if(method == "modwt") {
      for(k in names(zero)[-charmatch(names(zero)[j], names(zero))])
        zero[[k]] <- array(0, dim=c(nx,ny,nz))
      x.mra[[j]] <- imodwt.3d(zero)
    } else {
      for(k in 1:Jj)
        zero[[7*(k-1)+1]] <- zero[[7*(k-1)+2]] <- zero[[7*(k-1)+3]] <-
          zero[[7*(k-1)+4]] <- zero[[7*(k-1)+5]] <- zero[[7*(k-1)+6]] <-
            zero[[7*k]] <- array(0, dim=c(nx,ny,nz)/2^k)
      zero[[7*Jj+1]] <- array(0, dim=c(nx,ny,nz)/2^Jj)
      zero[[j]] <- x.wt[[j]]
      x.mra[[j]] <- idwt.3d(zero)
    }
  }
  return(x.mra)
}
