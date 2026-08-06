#' Discrete Wavelet Transform (DWT)
#' 
#' This function performs a level \eqn{J} decomposition of the input vector or
#' time series using the pyramid algorithm (Mallat 1989).
#' 
#' The code implements the one-dimensional DWT using the pyramid algorithm
#' (Mallat, 1989).  The actual transform is performed in C using pseudocode
#' from Percival and Walden (2001).  That means convolutions, not inner
#' products, are used to apply the wavelet filters.
#' 
#' For a non-dyadic length vector or time series, \code{dwt.nondyadic} pads
#' with zeros, performs the orthonormal DWT on this dyadic length series and
#' then truncates the wavelet coefficient vectors appropriately.
#' 
#' @usage dwt(x, wf = "la8", n.levels = 4, boundary = "periodic")
#' @usage dwt.nondyadic(x)
#' @usage idwt(y)
#' @aliases dwt dwt.nondyadic idwt
#' @param x a vector or time series containing the data be to decomposed.  This
#' must be a dyadic length vector (power of 2).
#' @param wf Name of the wavelet filter to use in the decomposition.  By
#' default this is set to \code{"la8"}, the Daubechies orthonormal compactly
#' supported wavelet of length L=8 (Daubechies, 1992), least asymmetric family.
#' @param n.levels Specifies the depth of the decomposition.  This must be a
#' number less than or equal to log(length(x),2).
#' @param boundary Character string specifying the boundary condition.  If
#' \code{boundary=="periodic"} the default, then the vector you decompose is
#' assumed to be periodic on its defined interval,\cr if
#' \code{boundary=="reflection"}, the vector beyond its boundaries is assumed
#' to be a symmetric reflection of itself.
#' @param y An object of S3 class \code{dwt}.
#' @return Basically, a list with the following components 
#' \item{d?}{Wavelet coefficient vectors.} 
#' \item{s?}{Scaling coefficient vector.}
#' \item{wavelet}{Name of the wavelet filter used.} 
#' \item{boundary}{How the boundaries were handled.}
#' @author B. Whitcher
#' @seealso \code{\link{modwt}}, \code{\link{mra}}.
#' @references Daubechies, I. (1992) \emph{Ten Lectures on Wavelets}, CBMS-NSF
#' Regional Conference Series in Applied Mathematics, SIAM: Philadelphia.
#' 
#' Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An Introduction to
#' Wavelets and Other Filtering Methods in Finance and Economics}, Academic
#' Press.
#' 
#' Mallat, S. G. (1989) A theory for multiresolution signal decomposition: the
#' wavelet representation, \emph{IEEE Transactions on Pattern Analysis and
#' Machine Intelligence}, \bold{11}(7), 674--693.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @examples
#' 
#' ## Figures 4.17 and 4.18 in Gencay, Selcuk and Whitcher (2001).
#' data(ibm)     
#' ibm.returns <- diff(log(ibm))
#' ## Haar
#' ibmr.haar <- dwt(ibm.returns, "haar")
#' names(ibmr.haar) <- c("w1", "w2", "w3", "w4", "v4")
#' ## plot partial Haar DWT for IBM data
#' par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
#' plot.ts(ibm.returns, axes=FALSE, ylab="", main="(a)")
#' for(i in 1:4)
#'   plot.ts(up.sample(ibmr.haar[[i]], 2^i), type="h", axes=FALSE,
#'           ylab=names(ibmr.haar)[i])
#' plot.ts(up.sample(ibmr.haar$v4, 2^4), type="h", axes=FALSE,
#'         ylab=names(ibmr.haar)[5])
#' axis(side=1, at=seq(0,368,by=23), 
#'      labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
#' ## LA(8)
#' ibmr.la8 <- dwt(ibm.returns, "la8")
#' names(ibmr.la8) <- c("w1", "w2", "w3", "w4", "v4")
#' ## must shift LA(8) coefficients
#' ibmr.la8$w1 <- c(ibmr.la8$w1[-c(1:2)], ibmr.la8$w1[1:2])
#' ibmr.la8$w2 <- c(ibmr.la8$w2[-c(1:2)], ibmr.la8$w2[1:2])
#' for(i in names(ibmr.la8)[3:4])
#'   ibmr.la8[[i]] <- c(ibmr.la8[[i]][-c(1:3)], ibmr.la8[[i]][1:3])
#' ibmr.la8$v4 <- c(ibmr.la8$v4[-c(1:2)], ibmr.la8$v4[1:2])
#' ## plot partial LA(8) DWT for IBM data
#' par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
#' plot.ts(ibm.returns, axes=FALSE, ylab="", main="(b)")
#' for(i in 1:4)
#'   plot.ts(up.sample(ibmr.la8[[i]], 2^i), type="h", axes=FALSE,
#'           ylab=names(ibmr.la8)[i])
#' plot.ts(up.sample(ibmr.la8$v4, 2^4), type="h", axes=FALSE,
#'         ylab=names(ibmr.la8)[5])
#' axis(side=1, at=seq(0,368,by=23), 
#'   labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
#' 
#' @export dwt
dwt <- function(x, wf="la8", n.levels=4, boundary="periodic")
{
  switch(boundary,
    "reflection" =  x <- c(x, rev(x)),
    "periodic" = invisible(),
    stop("Invalid boundary rule in dwt"))
  N <- length(x)
  J <- n.levels
  if(N/2^J != trunc(N/2^J))
    stop("Sample size is not divisible by 2^J")
  if(2^J > N)
    stop("wavelet transform exceeds sample size in dwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  y <- vector("list", J+1)
  names(y) <- c(paste("d", 1:J, sep=""), paste("s", J, sep=""))
  for(j in 1:J) {
    W <- V <- numeric(N/2^j)
    out <- .C(C_dwt, as.double(x), as.integer(N/2^(j-1)), L, h, g, 
              W=as.double(W), V=as.double(V))[6:7]
    y[[j]] <- out$W
    x <- out$V
  }
  y[[J+1]] <- x
  class(y) <- "dwt"
  attr(y, "wavelet") <- wf
  attr(y, "boundary") <- boundary
  return(y)
}

dwt.nondyadic <- function(x)
{
  M <- length(x)
  N <- 2^(ceiling(log(M, 2)))
  xx <- c(x, rep(0, N - M))
  y <- dwt(xx)
  
  J <- length(y) - 1
  for(j in 1:J)
    y[[j]] <- y[[j]][1:trunc(M/2^j)]
  return(y)
}

idwt <- function(y)
{
  ctmp <- class(y)
  if(is.null(ctmp) || all(ctmp != "dwt"))
    stop("argument `y' is not of class \"dwt\"")

  J <- length(y) - 1

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  jj <- paste("s", J, sep="")
  X <- y[[jj]]
  for(j in J:1) {
    jj <- paste("d", j, sep="")
    N <- length(X)
    XX <- numeric(2 * length(y[[jj]]))
    X <- .C(C_idwt, as.double(y[[jj]]), as.double(X), as.integer(N), L, 
            h, g, out=as.double(XX))$out
  }
  if(attr(y, "boundary") == "reflection") return(X[1:N])
  else return(X)
}



#' (Inverse) Maximal Overlap Discrete Wavelet Transform
#' 
#' This function performs a level \eqn{J} decomposition of the input vector
#' using the non-decimated discrete wavelet transform. The inverse transform
#' performs the reconstruction of a vector or time series from its maximal
#' overlap discrete wavelet transform.
#' 
#' The code implements the one-dimensional non-decimated DWT using the pyramid
#' algorithm.  The actual transform is performed in C using pseudocode from
#' Percival and Walden (2001).  That means convolutions, not inner products,
#' are used to apply the wavelet filters.
#' 
#' The MODWT goes by several names in the statistical and engineering
#' literature, such as, the ``stationary DWT'', ``translation-invariant DWT'',
#' and ``time-invariant DWT''.
#' 
#' The inverse MODWT implements the one-dimensional inverse transform using the
#' pyramid algorithm (Mallat, 1989).
#' 
#' @usage modwt(x, wf = "la8", n.levels = 4, boundary = "periodic")
#' @usage imodwt(y)
#' @aliases modwt imodwt
#' @param x a vector or time series containing the data be to decomposed.
#' There is \bold{no} restriction on its length.
#' @param y Object of class \code{"modwt"}.
#' @param wf Name of the wavelet filter to use in the decomposition.  By
#' default this is set to \code{"la8"}, the Daubechies orthonormal compactly
#' supported wavelet of length L=8 (Daubechies, 1992), least asymmetric family.
#' @param n.levels Specifies the depth of the decomposition.  This must be a
#' number less than or equal to log(length(x),2).
#' @param boundary Character string specifying the boundary condition.  If
#' \code{boundary=="periodic"} the defaulTRUE, then the vector you decompose is
#' assumed to be periodic on its defined interval,\cr if
#' \code{boundary=="reflection"}, the vector beyond its boundaries is assumed
#' to be a symmetric reflection of itself.
#' @param y an object of class \code{"modwt"}
#' @return Basically, a list with the following components 
#' \item{d?}{Wavelet coefficient vectors.} 
#' \item{s?}{Scaling coefficient vector.} 
#' \item{wavelet}{Name of the wavelet filter used.}
#' \item{boundary}{How the boundaries were handled.}
#' @author B. Whitcher
#' @seealso \code{\link{dwt}}, \code{\link{idwt}}, \code{\link{mra}}.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' 
#' Percival, D. B. and P. Guttorp (1994) Long-memory processes, the Allan
#' variance and wavelets, In \emph{Wavelets and Geophysics}, pages 325-344,
#' Academic Press.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @examples
#' 
#' ## Figure 4.23 in Gencay, Selcuk and Whitcher (2001)
#' data(ibm)     
#' ibm.returns <- diff(log(ibm))
#' # Haar
#' ibmr.haar <- modwt(ibm.returns, "haar")
#' names(ibmr.haar) <- c("w1", "w2", "w3", "w4", "v4")
#' # LA(8)
#' ibmr.la8 <- modwt(ibm.returns, "la8")
#' names(ibmr.la8) <- c("w1", "w2", "w3", "w4", "v4")
#' # shift the MODWT vectors
#' ibmr.la8 <- phase.shift(ibmr.la8, "la8")
#' ## plot partial MODWT for IBM data
#' par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
#' plot.ts(ibm.returns, axes=FALSE, ylab="", main="(a)")
#' for(i in 1:5)
#'   plot.ts(ibmr.haar[[i]], axes=FALSE, ylab=names(ibmr.haar)[i])
#' axis(side=1, at=seq(0,368,by=23), 
#'   labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
#' par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
#' plot.ts(ibm.returns, axes=FALSE, ylab="", main="(b)")
#' for(i in 1:5)
#'   plot.ts(ibmr.la8[[i]], axes=FALSE, ylab=names(ibmr.la8)[i])
#' axis(side=1, at=seq(0,368,by=23), 
#'   labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
#' 
#' @export modwt
modwt <- function(x, wf="la8", n.levels=4, boundary="periodic")
{
  switch(boundary,
    "reflection" =  x <- c(x, rev(x)),
    "periodic" = invisible(),
    stop("Invalid boundary rule in modwt"))
  N <- length(x)
  storage.mode(N) <- "integer"
  J <- n.levels
  if(2^J > N)
    stop("wavelet transform exceeds sample size in modwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  ht <- dict$hpf / sqrt(2)
  storage.mode(ht) <- "double"
  gt <- dict$lpf / sqrt(2)
  storage.mode(gt) <- "double"

  y <- vector("list", J+1)
  names(y) <- c(paste("d", 1:J, sep=""), paste("s", J, sep=""))
  W <- V <- numeric(N)
  storage.mode(W) <- "double"
  storage.mode(V) <- "double"
  
  for(j in 1:J) {
    out <- .C(C_modwt, as.double(x), N, as.integer(j), L, ht, gt, 
              W=W, V=V)[7:8]
    y[[j]] <- out$W
    x <- out$V
  }
  y[[J+1]] <- x
  class(y) <- "modwt"
  attr(y, "wavelet") <- wf
  attr(y, "boundary") <- boundary
  return(y)
}

imodwt <- function(y)
{
  ctmp <- class(y)
  if(is.null(ctmp) || all(ctmp != "modwt"))
    stop("argument `y' is not of class \"modwt\"")

  J <- length(y) - 1

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  ht <- dict$hpf / sqrt(2)
  storage.mode(ht) <- "double"
  gt <- dict$lpf / sqrt(2)
  storage.mode(gt) <- "double"

  jj <- paste("s", J, sep="")
  X <- y[[jj]]
  N <- length(X)
  storage.mode(N) <- "integer"
  XX <- numeric(N)
  storage.mode(XX) <- "double"
  for(j in J:1) {
    jj <- paste("d", j, sep="")
    X <- .C(C_imodwt, as.double(y[[jj]]), as.double(X), N, as.integer(j), 
            L, ht, gt, out=XX)$out
  }
  if(attr(y, "boundary") == "reflection") return(X[1:(N/2)])
  else return(X)
}



#' Replace Boundary Wavelet Coefficients with Missing Values
#' 
#' Sets the first \eqn{n} wavelet coefficients to \code{NA}.
#' 
#' The fact that observed time series are finite causes boundary issues.  One
#' way to get around this is to simply remove any wavelet coefficient computed
#' involving the boundary.  This is done here by replacing boundary wavelet
#' coefficients with \code{NA}.
#' 
#' @usage brick.wall(x, wf, method = "modwt")
#' @usage dwpt.brick.wall(x, wf, n.levels, method = "modwpt")
#' @usage brick.wall.2d(x, method = "modwt")
#' @aliases brick.wall dwpt.brick.wall brick.wall.2d
#' @param x DWT/MODWT/DWPT/MODWPT object
#' @param wf Character string; name of wavelet filter
#' @param n.levels Specifies the depth of the decomposition. This must be a
#' number less than or equal to log(length(x),2).
#' @param method Either \code{\link{dwt}} or \code{\link{modwt}} for
#' \code{brick.wall}, or either \code{\link{dwpt}} or \code{\link{modwpt}} for
#' \code{dwpt.brick.wall}
#' @return Same object as \code{x} only with some missing values.
#' @author B. Whitcher
#' @references Lindsay, R. W., D. B. Percival and D. A. Rothrock (1996). The
#' discrete wavelet transform and the scale anlaysis of the surface properties
#' of sea ice, \emph{IEEE Transactions on Geoscience and Remote Sensing},
#' \bold{34}, No. 3, 771-787.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @export brick.wall
brick.wall <- function(x, wf, method = "modwt")
{
  m <- wave.filter(wf)$length
  for (j in 1:(length(x) - 1)) {
    if (method == "dwt") {
      n <- ceiling((m - 2) * (1 - 1 / 2 ^ j))
    } else {
      n <- (2^j - 1) * (m - 1)
    }
    n <- min(n, length(x[[j]]))
    x[[j]][1:n] <- NA
  }
  x[[j + 1]][1:n] <- NA
  return(x)
}



#' Phase Shift Wavelet Coefficients
#' 
#' Wavelet coefficients are circularly shifted by the amount of phase shift
#' induced by the wavelet transform.
#' 
#' The center-of-energy argument of Hess-Nielsen and Wickerhauser (1996) is
#' used to provide a flexible way to circularly shift wavelet coefficients
#' regardless of the wavelet filter used.  The results are not identical to
#' those used by Percival and Walden (2000), but are more flexible.
#' 
#' \code{phase.shift.packet} is not yet implemented fully.
#' 
#' @usage phase.shift(z, wf, inv = FALSE)
#' @usage phase.shift.packet(z, wf, inv = FALSE)
#' @aliases phase.shift phase.shift.packet
#' @param z DWT object
#' @param wf character string; wavelet filter used in DWT
#' @param inv Boolean variable; if \code{inv=TRUE} then the inverse phase shift
#' is applied
#' @return DWT (DWPT) object with coefficients circularly shifted.
#' @author B. Whitcher
#' @references Hess-Nielsen, N. and M. V. Wickerhauser (1996) Wavelets and
#' time-frequency analysis, \emph{Proceedings of the IEEE}, \bold{84}, No. 4,
#' 523-540.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @export phase.shift
phase.shift <- function(z, wf, inv = FALSE)
{
  coe <- function(g)
    sum(0:(length(g)-1) * g^2) / sum(g^2)

  J <- length(z) - 1
  g <- wave.filter(wf)$lpf
  h <- wave.filter(wf)$hpf

  if(!inv) {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(z[[j]])
      z[[j]] <- c(z[[j]][(ph + 1):Nj], z[[j]][1:ph])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    z[[J]] <- c(z[[J]][(ph + 1):Nj], z[[J]][1:ph])
  } else {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(z[[j]])
      z[[j]] <- c(z[[j]][(Nj - ph + 1):Nj], z[[j]][1:(Nj - ph)])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    z[[J]] <- c(z[[J]][(Nj - ph + 1):Nj], z[[J]][1:(Nj - ph)])
  }
  return(z)
}



#' Multiresolution Analysis of Time Series
#' 
#' This function performs a level \eqn{J} additive decomposition of the input
#' vector or time series using the pyramid algorithm (Mallat 1989).
#' 
#' This code implements a one-dimensional multiresolution analysis introduced
#' by Mallat (1989).  Either the DWT or MODWT may be used to compute the
#' multiresolution analysis, which is an additive decomposition of the original
#' time series.
#' 
#' @param x A vector or time series containing the data be to decomposed.  This
#' must be a dyadic length vector (power of 2) for \code{method="dwt"}.
#' @param wf Name of the wavelet filter to use in the decomposition.  By
#' default this is set to \code{"la8"}, the Daubechies orthonormal compactly
#' supported wavelet of length L=8 least asymmetric family.
#' @param J Specifies the depth of the decomposition.  This must be a number
#' less than or equal to log(length(x), 2).
#' @param method Either \code{"dwt"} or \code{"modwt"}.
#' @param boundary Character string specifying the boundary condition.  If
#' \code{boundary=="periodic"} the default, then the vector you decompose is
#' assumed to be periodic on its defined interval,\cr if
#' \code{boundary=="reflection"}, the vector beyond its boundaries is assumed
#' to be a symmetric reflection of itself.
#' @return Basically, a list with the following components \item{D?}{Wavelet
#' detail vectors.} \item{S?}{Wavelet smooth vector.} \item{wavelet}{Name of
#' the wavelet filter used.} \item{boundary}{How the boundaries were handled.}
#' @author B. Whitcher
#' @seealso \code{\link{dwt}}, \code{\link{modwt}}.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' 
#' Mallat, S. G. (1989) A theory for multiresolution signal decomposition: the
#' wavelet representation, \emph{IEEE Transactions on Pattern Analysis and
#' Machine Intelligence}, \bold{11}, No. 7, 674-693.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' @keywords ts
#' @examples
#' 
#' ## Easy check to see if it works...
#' x <- rnorm(32)
#' x.mra <- mra(x)
#' sum(x - apply(matrix(unlist(x.mra), nrow=32), 1, sum))^2
#' 
#' ## Figure 4.19 in Gencay, Selcuk and Whitcher (2001)
#' data(ibm)     
#' ibm.returns <- diff(log(ibm))
#' ibm.volatility <- abs(ibm.returns)
#' ## Haar
#' ibmv.haar <- mra(ibm.volatility, "haar", 4, "dwt")
#' names(ibmv.haar) <- c("d1", "d2", "d3", "d4", "s4")
#' ## LA(8)
#' ibmv.la8 <- mra(ibm.volatility, "la8", 4, "dwt")
#' names(ibmv.la8) <- c("d1", "d2", "d3", "d4", "s4")
#' ## plot multiresolution analysis of IBM data
#' par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
#' plot.ts(ibm.volatility, axes=FALSE, ylab="", main="(a)")
#' for(i in 1:5)
#'   plot.ts(ibmv.haar[[i]], axes=FALSE, ylab=names(ibmv.haar)[i])
#' axis(side=1, at=seq(0,368,by=23), 
#'   labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
#' par(mfcol=c(6,1), pty="m", mar=c(5-2,4,4-2,2))
#' plot.ts(ibm.volatility, axes=FALSE, ylab="", main="(b)")
#' for(i in 1:5)
#'   plot.ts(ibmv.la8[[i]], axes=FALSE, ylab=names(ibmv.la8)[i])
#' axis(side=1, at=seq(0,368,by=23), 
#'   labels=c(0,"",46,"",92,"",138,"",184,"",230,"",276,"",322,"",368))
#' 
#' @export mra
mra <- function(x, wf = "la8", J = 4, method = "modwt", boundary = "periodic")
{
  switch(boundary,
         "reflection" =  x <- c(x, rev(x)),
         "periodic" = invisible(),
         stop("Invalid boundary rule in mra"))
  n <- length(x)
  
  if(method == "modwt")
    x.wt <- modwt(x, wf, J, "periodic")
  else
    x.wt <- dwt(x, wf, J, "periodic")
  x.mra <- vector("list", J+1)
  
  ## Smooth
  zero <- vector("list", J+1)
  names(zero) <- c(paste("d", 1:J, sep = ""), paste("s", J, sep = ""))
  class(zero) <- method
  attr(zero, "wavelet") <- wf
  attr(zero, "boundary") <- boundary
  zero[[J+1]] <- x.wt[[J+1]]
  if(method == "modwt") {
    for(k in 1:J)
      zero[[k]] <- numeric(n)
    x.mra[[J+1]] <- imodwt(zero)
  } else {
    for(k in 1:J)
      zero[[k]] <- numeric(n/2^k)
    x.mra[[J+1]] <- idwt(zero)
  }

  ## Details
  for(j in J:1) {
    zero <- vector("list", j+1)
    names(zero) <- c(paste("d", 1:j, sep = ""), paste("s", j, sep = ""))
    class(zero) <- method
    attr(zero, "wavelet") <- wf
    attr(zero, "boundary") <- boundary
    zero[[j]] <- x.wt[[j]]
    if(method == "modwt") {
      if(j != 1) {
        for(k in c(j+1,(j-1):1))
          zero[[k]] <- numeric(n)
      } else {
        zero[[j+1]] <- numeric(n)
      }
      x.mra[[j]] <- imodwt(zero)
    } else {
      zero[[j+1]] <- numeric(n/2^j)
      if(j != 1) {
        for(k in (j-1):1)
          zero[[k]] <- numeric(n/2^k)
      }
      x.mra[[j]] <- idwt(zero)
    }
  }

  names(x.mra) <- c(paste("D", 1:J, sep = ""), paste("S", J, sep = ""))
  if(boundary == "reflection") { 
    for(j in (J+1):1)
      x.mra[[j]] <- x.mra[[j]][1:(n/2)]
    return(x.mra)
  } else {
    return(x.mra)
  }
}

