#' Two-Dimensional Discrete Wavelet Transform
#' 
#' Performs a separable two-dimensional discrete wavelet transform (DWT) on a
#' matrix of dyadic dimensions.
#' 
#' See references.
#' 
#' @usage dwt.2d(x, wf, J = 4, boundary = "periodic")
#' @usage idwt.2d(y)
#' @aliases dwt.2d idwt.2d
#' @param x input matrix (image)
#' @param wf name of the wavelet filter to use in the decomposition
#' @param J depth of the decomposition, must be a number less than or equal to
#' log(minM,N,2)
#' @param boundary only \code{"periodic"} is currently implemented
#' @param y an object of class \code{dwt.2d}
#' @return List structure containing the \eqn{3J+1} sub-matrices from the
#' decomposition.
#' @author B. Whitcher
#' @seealso \code{\link{modwt.2d}}.
#' @references Mallat, S. (1998) \emph{A Wavelet Tour of Signal Processing},
#' Academic Press.
#' 
#' Vetterli, M. and J. Kovacevic (1995) \emph{Wavelets and Subband Coding},
#' Prentice Hall.
#' @keywords ts
#' @examples
#' 
#' ## Xbox image
#' data(xbox)
#' xbox.dwt <- dwt.2d(xbox, "haar", 3)
#' par(mfrow=c(1,1), pty="s")
#' plot.dwt.2d(xbox.dwt)
#' par(mfrow=c(2,2), pty="s")
#' image(1:dim(xbox)[1], 1:dim(xbox)[2], xbox, xlab="", ylab="",
#'       main="Original Image")
#' image(1:dim(xbox)[1], 1:dim(xbox)[2], idwt.2d(xbox.dwt), xlab="", ylab="",
#'       main="Wavelet Reconstruction")
#' image(1:dim(xbox)[1], 1:dim(xbox)[2], xbox - idwt.2d(xbox.dwt),
#'       xlab="", ylab="", main="Difference")
#' 
#' ## Daubechies image
#' data(dau)
#' par(mfrow=c(1,1), pty="s")
#' image(dau, col=rainbow(128))
#' sum(dau^2)
#' dau.dwt <- dwt.2d(dau, "d4", 3)
#' plot.dwt.2d(dau.dwt)
#' sum(plot.dwt.2d(dau.dwt, plot=FALSE)^2)
#' 
#' @export dwt.2d
dwt.2d <- function(x, wf, J=4, boundary="periodic")
{
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  z <- matrix(0, m/2, n/2)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 3*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C(C_two_D_dwt, "Image"=as.double(x), "Rows"=m, "Cols"=n, 
                "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z, "LH"=z,
                "HL"=z, "HH"=z)[7:10]
    if(j < J) {
      index <- (3*j-2):(3*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      m <- dim(x)[1]
      storage.mode(m) <- "integer"
      n <- dim(x)[2]
      storage.mode(n) <- "integer"
      z <- matrix(0, m/2, n/2)
      storage.mode(z) <- "double"
    }
    else {
      index <- (3*j):(3*(j+1)) - 2
      x.wt[index] <- out[c(2:4,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  attr(x.wt, "class") <- "dwt.2d"
  x.wt
}

idwt.2d <- function(y)
{
  J <- attributes(y)$J

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  LL <- paste("LL", J, sep="")
  y.in <- y[[LL]]

  for(j in J:1) {
    LH <- paste("LH", j, sep="")
    HL <- paste("HL", j, sep="")
    HH <- paste("HH", j, sep="")
    
    m <- dim(y.in)[1]
    storage.mode(m) <- "integer"
    n <- dim(y.in)[2]
    storage.mode(n) <- "integer"
    x <- matrix(0, 2*m, 2*n)
    storage.mode(x) <- "double"

    out <- .C(C_two_D_idwt, as.double(y.in), as.double(y[[LH]]),
              as.double(y[[HL]]), as.double(y[[HH]]), m, n, L, h, g,
              "Y"=x)
    y.in <- out$Y
  }
  zapsmall(y.in)
}



#' Two-Dimensional Maximal Overlap Discrete Wavelet Transform
#' 
#' Performs a separable two-dimensional maximal overlap discrete wavelet
#' transform (MODWT) on a matrix of arbitrary dimensions.
#' 
#' See references.
#' 
#' @usage modwt.2d(x, wf, J = 4, boundary = "periodic")
#' @usage imodwt.2d(y)
#' @aliases modwt.2d imodwt.2d
#' @param x input matrix
#' @param wf name of the wavelet filter to use in the decomposition
#' @param J depth of the decomposition
#' @param boundary only \code{"periodic"} is currently implemented
#' @param y an object of class \code{dwt.2d}
#' @return List structure containing the \eqn{3J+1} sub-matrices from the
#' decomposition.
#' @author B. Whitcher
#' @seealso \code{\link{dwt.2d}}, \code{\link{shift.2d}}.
#' @references Liang, J. and T. W. Parks (1994) A two-dimensional translation
#' invariant wavelet representation and its applications, \emph{Proceedings
#' ICIP-94}, Vol. 1, 66-70.
#' 
#' Liang, J. and T. W. Parks (1994) Image coding using translation invariant
#' wavelet transforms with symmetric extensions, \emph{IEEE Transactions on
#' Image Processing}, \bold{7}, No. 5, 762-769.
#' @keywords ts
#' @examples
#' 
#' ## Xbox image
#' data(xbox)
#' xbox.modwt <- modwt.2d(xbox, "haar", 2)
#' ## Level 1 decomposition
#' par(mfrow=c(2,2), pty="s")
#' image(xbox.modwt$LH1, col=rainbow(128), axes=FALSE, main="LH1")
#' image(xbox.modwt$HH1, col=rainbow(128), axes=FALSE, main="HH1")
#' frame()
#' image(xbox.modwt$HL1, col=rainbow(128), axes=FALSE, main="HL1")
#' ## Level 2 decomposition
#' par(mfrow=c(2,2), pty="s")
#' image(xbox.modwt$LH2, col=rainbow(128), axes=FALSE, main="LH2")
#' image(xbox.modwt$HH2, col=rainbow(128), axes=FALSE, main="HH2")
#' image(xbox.modwt$LL2, col=rainbow(128), axes=FALSE, main="LL2")
#' image(xbox.modwt$HL2, col=rainbow(128), axes=FALSE, main="HL2")
#' sum((xbox - imodwt.2d(xbox.modwt))^2)
#' 
#' data(dau)
#' par(mfrow=c(1,1), pty="s")
#' image(dau, col=rainbow(128), axes=FALSE, main="Ingrid Daubechies")
#' sum(dau^2)
#' dau.modwt <- modwt.2d(dau, "d4", 2)
#' ## Level 1 decomposition
#' par(mfrow=c(2,2), pty="s")
#' image(dau.modwt$LH1, col=rainbow(128), axes=FALSE, main="LH1")
#' image(dau.modwt$HH1, col=rainbow(128), axes=FALSE, main="HH1")
#' frame()
#' image(dau.modwt$HL1, col=rainbow(128), axes=FALSE, main="HL1")
#' ## Level 2 decomposition
#' par(mfrow=c(2,2), pty="s")
#' image(dau.modwt$LH2, col=rainbow(128), axes=FALSE, main="LH2")
#' image(dau.modwt$HH2, col=rainbow(128), axes=FALSE, main="HH2")
#' image(dau.modwt$LL2, col=rainbow(128), axes=FALSE, main="LL2")
#' image(dau.modwt$HL2, col=rainbow(128), axes=FALSE, main="HL2")
#' sum((dau - imodwt.2d(dau.modwt))^2)
#' 
#' @export modwt.2d
modwt.2d <- function(x, wf, J=4, boundary="periodic")
{
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  z <- matrix(0, m, n)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 3*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("two_D_modwt", "Image"=as.double(x), "Rows"=m, "Cols"=n,
              "Level"=j, "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z,
              "LH"=z, "HL"=z, "HH"=z, PACKAGE="waveslim")[8:11]
    if(j < J) {
      index <- (3*j-2):(3*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out$LL
    }
    else {
      index <- (3*j):(3*(j+1)) - 2
      x.wt[index] <- out[c(2:4,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  x.wt
}


imodwt.2d <- function(y)
{
  J <- attributes(y)$J

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf / sqrt(2)
  storage.mode(h) <- "double"
  g <- dict$lpf / sqrt(2)
  storage.mode(g) <- "double"

  LL <- paste("LL", J, sep="")
  y.in <- y[[LL]]

  for(j in J:1) {
    LH <- paste("LH", j, sep="")
    HL <- paste("HL", j, sep="")
    HH <- paste("HH", j, sep="")
    
    m <- dim(y.in)[1]
    storage.mode(m) <- "integer"
    n <- dim(y.in)[2]
    storage.mode(n) <- "integer"
    x <- matrix(0, m, n)
    storage.mode(x) <- "double"

    out <- .C(C_two_D_imodwt, as.double(y.in), as.double(y[[LH]]),
              as.double(y[[HL]]), as.double(y[[HH]]), m, n, j, L,
              h, g, "Y"=x)
    y.in <- out$Y
  }
  zapsmall(y.in)
}



#' Plot Two-dimensional Discrete Wavelet Transform
#' 
#' Organizes the wavelet coefficients from a 2D DWT into a single matrix and
#' plots it.  The coarser resolutions are nested within the lower-lefthand
#' corner of the image.
#' 
#' The wavelet coefficients from the DWT object (a list) are reorganized into a
#' single matrix of the same dimension as the original image and the result is
#' plotted.
#' 
#' @param x input matrix (image)
#' @param cex.axis \code{par} plotting parameter that controls the size of the
#' axis text
#' @param plot if \code{plot = FALSE} then the matrix of wavelet coefficients
#' is returned, the default is \code{plot = TRUE}
#' @param ... additional graphical parameters if necessary
#' @return Image plot.
#' @author B. Whitcher
#' @seealso \code{\link{dwt.2d}}.
#' @keywords ts
#' @export plot.dwt.2d
plot.dwt.2d <- function(x, cex.axis=1, plot=TRUE, ...)
{
  J <- attributes(x)$J
  X <- x[[paste("LL", J, sep="")]]
  for(j in J:1) {
    x.names <- sapply(c("LH","HL","HH"), paste, j, sep="")
    X <- rbind(cbind(X, x[[x.names[2]]]),
               cbind(x[[x.names[1]]], x[[x.names[3]]]))
  }
  M <- dim(X)[1]; N <- dim(X)[2]
  if(plot) {
    image(1:M, 1:N, X, col=rainbow(128), axes=FALSE, xlab="", ylab="", ...)
    x.label <- NULL
    lines(c(0,N,N,0,0) + 0.5, c(0,0,M,M,0) + 0.5)
    for(j in J:1) {
      lines(c(M/2^j,M/2^j) + 0.5, 2*c(0,N/2^j) + 0.5)
      lines(2*c(0,M/2^j) + 0.5, c(N/2^j,N/2^j) + 0.5)
    }
    at <- c((3*N+2)/2^(1:J+1),(N+2)/2^(J+1))
    labs <- c(paste("H",1:J,sep=""), paste("L",J,sep=""))
    axis(side=1, at=at, labels=labs, tick=FALSE, cex.axis=cex.axis)
    axis(side=2, at=at, labels=labs, tick=FALSE, cex.axis=cex.axis)
  }
  else
    return(X)
  invisible()
}



#' Denoise an Image via the 2D Discrete Wavelet Transform
#' 
#' Perform simple de-noising of an image using the two-dimensional discrete
#' wavelet transform.
#' 
#' See \code{\link{Thresholding}}.
#' 
#' @aliases denoise.dwt.2d denoise.modwt.2d
#' @param x input matrix (image)
#' @param wf name of the wavelet filter to use in the decomposition
#' @param J depth of the decomposition, must be a number less than or equal to
#' log(minM,N,2)
#' @param method character string describing the threshold applied, only
#' \code{"universal"} and \code{"long-memory"} are currently implemented
#' @param H self-similarity or Hurst parameter to indicate spectral scaling,
#' white noise is 0.5
#' @param noise.dir number of directions to estimate background noise standard
#' deviation, the default is 3 which produces a unique estimate of the
#' background noise for each spatial direction
#' @param rule either a \code{"hard"} or \code{"soft"} thresholding rule may be
#' used
#' @return Image of the same dimension as the original but with high-freqency
#' fluctuations removed.
#' @author B. Whitcher
#' @seealso \code{\link{Thresholding}}
#' @references See \code{\link{Thresholding}} for references concerning
#' de-noising in one dimension.
#' @keywords ts
#' @examples
#' 
#' ## Xbox image
#' data(xbox)
#' n <- nrow(xbox)
#' xbox.noise <- xbox + matrix(rnorm(n*n, sd=.15), n, n)
#' par(mfrow=c(2,2), cex=.8, pty="s")
#' image(xbox.noise, col=rainbow(128), main="Original Image")
#' image(denoise.dwt.2d(xbox.noise, wf="haar"), col=rainbow(128),
#'       zlim=range(xbox.noise), main="Denoised image")
#' image(xbox.noise - denoise.dwt.2d(xbox.noise, wf="haar"), col=rainbow(128),
#'       zlim=range(xbox.noise), main="Residual image")
#' 
#' ## Daubechies image
#' data(dau)
#' n <- nrow(dau)
#' dau.noise <- dau + matrix(rnorm(n*n, sd=10), n, n)
#' par(mfrow=c(2,2), cex=.8, pty="s")
#' image(dau.noise, col=rainbow(128), main="Original Image")
#' dau.denoise <- denoise.modwt.2d(dau.noise, wf="d4", rule="soft")
#' image(dau.denoise, col=rainbow(128), zlim=range(dau.noise),
#'       main="Denoised image")
#' image(dau.noise - dau.denoise, col=rainbow(128), main="Residual image")
#' 
denoise.dwt.2d <- function(x, wf = "la8", J = 4, method = "universal", 
                           H = 0.5, noise.dir = 3, rule = "hard")
{
  soft <- function(x, delta) sign(x) * pmax(abs(x) - delta, 0)
  hard <- function(x, delta) ifelse(abs(x) > delta, x, 0)

  n <- length(x)
  x.dwt <- dwt.2d(x, wf, J)
  if(noise.dir == 3)
    sigma.mad <- list(HH = mad(x.dwt$HH1), HL = mad(x.dwt$HL1), 
  		      LH = mad(x.dwt$LH1))
  else {
    noise <- x.dwt$jj
    sigma.mad <- list(HH = mad(noise), HL = mad(noise), LH = mad(noise))
  }
    thresh <- list(HH = rep(sqrt(2 * sigma.mad$HH^2 * log(n)), J), 
		   HL = rep(sqrt(2 * sigma.mad$HL^2 * log(n)), J),
		   LH = rep(sqrt(2 * sigma.mad$LH^2 * log(n)), J))

  if(method == "long-memory")
    thresh <- lapply(thresh, function(x,J,H) 2^(0:(J-1)*(H-1/2))*x, J=J, H=H)
  for(j in 1:J) {
    jj <- paste("HL", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$HL[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$HL[j])
    jj <- paste("LH", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$LH[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$LH[j])
    jj <- paste("HH", j, sep = "")
    if(rule == "hard")
      x.dwt[[jj]] <- hard(x.dwt[[jj]], thresh$HH[j])
    else 
      x.dwt[[jj]] <- soft(x.dwt[[jj]], thresh$HH[j])
  }
  idwt.2d(x.dwt)
}


denoise.modwt.2d <- function(x, wf = "la8", J = 4, method = "universal", 
  H = 0.5, rule = "hard")
{
  soft <- function(x, delta) sign(x) * pmax(abs(x) - delta, 0)
  hard <- function(x, delta) ifelse(abs(x) > delta, x, 0)
  n <- length(x)
  x.modwt <- modwt.2d(x, wf, J)
  sigma.mad <- list(HH = sqrt(2) * mad(x.modwt$HH1),
		    HL = sqrt(2) * mad(x.modwt$HL1),
		    LH = sqrt(2) * mad(x.modwt$LH1))
    thresh <- list(HH = rep(sqrt(2 * sigma.mad$HH^2 * log(n))/2^(1:J), J), 
		   HL = rep(sqrt(2 * sigma.mad$HL^2 * log(n))/2^(1:J), J), 
		   LH = rep(sqrt(2 * sigma.mad$LH^2 * log(n))/2^(1:J), J))
  if(method == "long-memory")
    thresh <- lapply(thresh, function(x,J,H) 2^(0:(J-1)*(H-1/2))*x, J=J, H=H)
  for(j in 1:J) {
    jj <- paste("HL", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$HL[j])
    else 
      x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$HL[j])
    jj <- paste("LH", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$LH[j])
    else 
      x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$LH[j])
    jj <- paste("HH", j, sep = "")
    if(rule == "hard")
      x.modwt[[jj]] <- hard(x.modwt[[jj]], thresh$HH[j])
    else 
     x.modwt[[jj]] <- soft(x.modwt[[jj]], thresh$HH[j])
  }
  imodwt.2d(x.modwt)
}



#' (Inverse) Discrete Wavelet Packet Transforms in Two Dimensions
#' 
#' All possible filtering combinations (low- and high-pass) are performed to
#' decompose a matrix or image.  The resulting coefficients are associated with
#' a quad-tree structure corresponding to a partitioning of the two-dimensional
#' frequency plane.
#' 
#' The code implements the two-dimensional DWPT using the pyramid algorithm of
#' Mallat (1989).
#' 
#' @usage dwpt.2d(x, wf = "la8", J = 4, boundary = "periodic")
#' @usage idwpt.2d(y, y.basis)
#' @aliases dwpt.2d idwpt.2d
#' @param x a matrix or image containing the data be to decomposed.  This
#' ojbect must be dyadic (power of 2) in length in each dimension.
#' @param wf Name of the wavelet filter to use in the decomposition.  By
#' default this is set to \code{"la8"}, the Daubechies orthonormal compactly
#' supported wavelet of length \eqn{L=8} (Daubechies, 1992), least asymmetric
#' family.
#' @param J Specifies the depth of the decomposition.  This must be a number
#' less than or equal to \eqn{\log(\mbox{length}(x),2)}.
#' @param boundary Character string specifying the boundary condition.  If
#' \code{boundary=="periodic"} the default, then the vector you decompose is
#' assumed to be periodic on its defined interval,\cr if
#' \code{boundary=="reflection"}, the vector beyond its boundaries is assumed
#' to be a symmetric reflection of itself.
#' @param y \code{dwpt.2d} object (list-based structure of matrices)
#' @param y.basis Boolean vector, the same length as \eqn{y}, where \code{TRUE}
#' means the basis tensor should be used in the reconstruction.
#' @return Basically, a list with the following components
#' \item{w?.?-w?.?}{Wavelet coefficient matrices (images).  The first index is
#' associated with the scale of the decomposition while the second is
#' associated with the frequency partition within that level.  The left and
#' right strings, separated by the dash `-', correspond to the first \eqn{(x)}
#' and second \eqn{(y)} dimensions.} \item{wavelet}{Name of the wavelet filter
#' used.} \item{boundary}{How the boundaries were handled.}
#' @author B. Whitcher
#' @seealso \code{\link{dwt.2d}}, \code{\link{modwt.2d}},
#' \code{\link{wave.filter}}.
#' @references Mallat, S. G. (1989) A theory for multiresolution signal
#' decomposition: the wavelet representation, \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, \bold{11}, No. 7, 674-693.
#' 
#' Wickerhauser, M. V. (1994) \emph{Adapted Wavelet Analysis from Theory to
#' Software}, A K Peters.
#' @keywords ts
#' @export dwpt.2d
dwpt.2d <- function(x, wf="la8", J=4, boundary="periodic")
{
  ## x <- xbox
  ## Define image dimensions (assign mode for C) and perform simple
  ## diagnostics.
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"
  if(log(m, 2) != trunc(log(m, 2)) | log(n, 2) != trunc(log(n, 2)))
    stop("One dimension is not a power of 2")
  if(2^J > m | 2^J > n)
    stop("Wavelet transform exceeds sample size in one dimension of DWPT")

  ## Extract wavelet and scaling filter coefficients, along with filter
  ## length, from the filter name provided.
  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  ## Create names for wavelet packet nodes (quad-tree structure).
  N <- sum(4^(1:J))
  level <- rep(1:J, 4^(1:J))
  x.wpt <- vector("list", N)
  c1 <- rep(1:J, 2^(1:J))
  c2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  cry <- paste("w", c1, ".", c2, sep="")
  x.wpt.names <- NULL
  for(j in 1:J) {
    xx <- matrix(cry[c1 == j], 2^j, 2^j)
    yy <- matrix(cry[c1 == j], 2^j, 2^j, byrow=TRUE)
    x.wpt.names <- c(x.wpt.names, as.matrix(paste(xx, "-", yy, sep="")))
  }
  names(x.wpt) <- x.wpt.names
  rm(j,xx,yy,c1,c2,cry)
  ## Define initial zero matrix to store wavelet sub-images.
  z <- matrix(0, m/2, n/2)
  storage.mode(z) <- "double"

  ## Implement the 2D DWPT in a nested loop structure.
  for(j in 1:J) {
    ## cat("j =", j, fill=TRUE)
    for(k in 0:(4^j/4-1)) {
      if(j > 1) {
        ## if j > 1, grab wavelet coefficient image and also its name.
        index <- min((1:N)[level == j-1]) + k
        parent <- x.wpt.names[index]
        ## cat("parent =", parent, fill=TRUE)
        x <- x.wpt[[parent]]
        tmp <- unlist(strsplit(parent, "\\-"))
      }
      else
        tmp <- c("w0.0", "w0.0")
      ## Deconstruct name into nodes for the x and y dimensions.
      node <- unlist(strsplit(tmp, "\\."))
      node <- as.integer(node[-c(1,3)])
      ## Preliminary assignments in order to keep wavelet coefficient
      ## sub-images in sequency order.
      if(node[1] %% 2 == 0) {
        Xlow <- paste("w", j, ".", 2 * node[1], sep="")
        Xhigh <- paste("w", j, ".", 2 * node[1] + 1, sep="")
      }
      else {
        Xlow <- paste("w", j, ".", 2 * node[1] + 1, sep="")
        Xhigh <- paste("w", j, ".", 2 * node[1], sep="")
      }
      if(node[2] %% 2 == 0) {
        Ylow <- paste("w", j, ".", 2 * node[2], sep="")
        Yhigh <- paste("w", j, ".", 2 * node[2] + 1, sep="")
      }
      else {
        Ylow <- paste("w", j, ".", 2 * node[2] + 1, sep="")
        Yhigh <- paste("w", j, ".", 2 * node[2], sep="")
      }
      ## Create names for the new wavelet coefficient images.
      LL <- paste(Xlow, "-", Ylow, sep="")
      LH <- paste(Xlow, "-", Yhigh, sep="")
      HL <- paste(Xhigh, "-", Ylow, sep="")
      HH <- paste(Xhigh, "-", Yhigh, sep="")
      ## cat(matrix(c(LH,LL,HH,HL), 2, 2), fill=TRUE)
      ## Perform the DWPT
      out <- .C(C_two_D_dwt, "Image"=as.double(x), "Rows"=m, "Cols"=n, 
                "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z, "LH"=z,
                "HL"=z, "HH"=z)[7:10]
      ## Pass wavelet coefficient images into the DWPT object.
      x.wpt[[LL]] <- out[["LL"]]
      x.wpt[[LH]] <- out[["LH"]]
      x.wpt[[HL]] <- out[["HL"]]
      x.wpt[[HH]] <- out[["HH"]]
    }
    ## Redefine zero matrix to its new (decimated) size.
    m <- dim(out[["LL"]])[1]
    storage.mode(m) <- "integer"
    n <- dim(out[["LL"]])[2]
    storage.mode(n) <- "integer"
    z <- matrix(0, m/2, n/2)
    storage.mode(z) <- "double"
  }
  attr(x.wpt, "J") <- J
  attr(x.wpt, "wavelet") <- wf
  attr(x.wpt, "boundary") <- boundary
  return(x.wpt)
}


idwpt.2d <- function(y, y.basis)
{
  ## Error checking
  if(length(y) != length(y.basis))
    stop("DWPT object and basis selection must be the same length")
  ## Number of wavelet scales
  J <- attributes(y)$J
  ## Define wavelet/scaling filter coefficients and length
  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"
  ## Nested for loops
  names(y.basis) <- names(y)
  for(j in J:1) {
    for(nx in seq(0, 2^j - 1, by = 2)) {
      for(ny in seq(0, 2^j - 1, by = 2)) {
        ## Name the four wavelet coefficients sub-images
        LL <- paste("w", j, ".", nx, "-", "w", j, ".", ny, sep="")
        LH <- paste("w", j, ".", nx, "-", "w", j, ".", ny+1, sep="")
        HL <- paste("w", j, ".", nx+1, "-", "w", j, ".", ny, sep="")
        HH <- paste("w", j, ".", nx+1, "-", "w", j, ".", ny+1, sep="")
        if(any(y.basis[LL], y.basis[LH], y.basis[HL], y.basis[HH])) {
          m <- nrow(y[[LL]])
          storage.mode(m) <- "integer"
          n <- ncol(y[[LL]])
          storage.mode(n) <- "integer"
          XX <- matrix(0, 2*m, 2*n)
          storage.mode(XX) <- "double"
          ## parent indices to construct string
          pnx <- floor(nx / 2)
          pny <- floor(ny / 2)
          if((pnx %% 2 != 0) & (pny %% 2 != 0))
            ## Upper right-hand corner
            out <- .C(C_two_D_idwt, as.double(y[[HH]]),
                      as.double(y[[HL]]), as.double(y[[LH]]),
                      as.double(y[[LL]]), m, n, L, h, g, "Y"=XX)$Y
          else {
            ## Upper left-hand corner
            if((pnx %% 2 == 0) & (pny %% 2 != 0))
              out <- .C(C_two_D_idwt, as.double(y[[LH]]),
                        as.double(y[[LL]]), as.double(y[[HH]]),
                        as.double(y[[HL]]), m, n, L, h, g, "Y"=XX)$Y
            else {
              ## Lower right-hand corner
              if((pnx %% 2 != 0) & (pny %% 2 == 0))
                out <- .C(C_two_D_idwt, as.double(y[[HL]]),
                          as.double(y[[HH]]), as.double(y[[LL]]),
                          as.double(y[[LH]]), m, n, L, h, g, "Y"=XX)$Y
              else {
                ## Lower left-hand corner
                if((pnx %% 2 == 0) & (pny %% 2 == 0))
                  out <- .C(C_two_D_idwt, as.double(y[[LL]]),
                            as.double(y[[LH]]), as.double(y[[HL]]),
                            as.double(y[[HH]]), m, n, L, h, g, "Y"=XX)$Y
                else
                  stop("Ouch!")
              }
            }
          }
          if(j > 1) {
            pname <- paste("w", j-1, ".", pnx, "-", "w", j-1, ".", pny, sep="")
            y[[pname]] <- out
            y.basis[pname] <- 1
          }
        }
      }
    }
  }
  return(out)
}

