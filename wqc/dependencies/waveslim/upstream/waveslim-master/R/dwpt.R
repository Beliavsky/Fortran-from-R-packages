#' (Inverse) Discrete Wavelet Packet Transforms
#' 
#' All possible filtering combinations (low- and high-pass) are performed to
#' decompose a vector or time series.  The resulting coefficients are
#' associated with a binary tree structure corresponding to a partitioning of
#' the frequency axis.
#' 
#' The code implements the one-dimensional DWPT using the pyramid algorithm
#' (Mallat, 1989).
#' 
#' @usage dwpt(x, wf = "la8", n.levels = 4, boundary = "periodic")
#' @usage idwpt(y, y.basis)
#' @aliases dwpt idwpt modwpt
#' @param x a vector or time series containing the data be to decomposed. This
#' must be a dyadic length vector (power of 2).
#' @param wf Name of the wavelet filter to use in the decomposition. By
#' default this is set to \code{"la8"}, the Daubechies orthonormal compactly
#' supported wavelet of length L=8 (Daubechies, 1992), least asymmetric family.
#' @param n.levels Specifies the depth of the decomposition.This must be a
#' number less than or equal to
#' \eqn{\log(\mbox{length}(x),2)}{log2[length(x)]}.
#' @param boundary Character string specifying the boundary condition. If
#' \code{boundary=="periodic"} the default, then the vector you decompose is
#' assumed to be periodic on its defined interval,\cr if
#' \code{boundary=="reflection"}, the vector beyond its boundaries is assumed
#' to be a symmetric reflection of itself.
#' @param y Object of S3 class \code{dwpt}.
#' @param y.basis Vector of character strings that describe leaves on the DWPT 
#' basis tree.
#' @return Basically, a list with the following components 
#' \item{w?.?}{Wavelet coefficient vectors.  The first index is associated with 
#' the scale of the decomposition while the second is associated with the 
#' frequency partition within that level.} 
#' \item{wavelet}{Name of the wavelet filter used.}
#' \item{boundary}{How the boundaries were handled.}
#' @author B. Whitcher
#' @seealso \code{\link{dwt}}, \code{\link{modwpt}}, \code{\link{wave.filter}}.
#' @references Mallat, S. G. (1989) A theory for multiresolution signal
#' decomposition: the wavelet representation, \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, \bold{11}(7), 674--693.
#' 
#' Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
#' Series Analysis}, Cambridge University Press.
#' 
#' Wickerhauser, M. V. (1994) \emph{Adapted Wavelet Analysis from Theory to
#' Software}, A K Peters.
#' @keywords ts
#' @examples
#' 
#' data(mexm)
#' J <- 4
#' mexm.mra <- mra(log(mexm), "mb8", J, "modwt", "reflection")
#' mexm.nomean <- ts(
#'   apply(matrix(unlist(mexm.mra), ncol=J+1, byrow=FALSE)[,-(J+1)], 1, sum), 
#'   start=1957, freq=12)
#' mexm.dwpt <- dwpt(mexm.nomean[-c(1:4)], "mb8", 7, "reflection")
#' 
#' @export dwpt
dwpt <- function(x, wf="la8", n.levels=4, boundary="periodic") {
  N <- length(x)
  J <- n.levels
  if(N/2^J != trunc(N/2^J))
    stop("Sample size is not a power of 2")
  if(2^J > N)
    stop("wavelet transform exceeds sample size in dwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  y <- vector("list", sum(2^(1:J)))
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  names(y) <- paste("w", crystals1, ".", crystals2, sep="")

  for(j in 1:J) {
    jj <- min((1:length(crystals1))[crystals1 == j])
    for(n in 0:(2^j/2-1)) {
      if(j > 1)
        x <- y[[(1:length(crystals1))[crystals1 == j-1][n+1]]]
      W <- V <- numeric(N/2^j)
      if(n %% 2 == 0) {
        z <- .C(C_dwt, as.double(x), as.integer(N/2^(j-1)), L, h, g, 
	  W=as.double(W), V=as.double(V))
        y[[jj + 2*n + 1]] <- z$W
        y[[jj + 2*n]] <- z$V
      }
      else {
        z <- .C(C_dwt, as.double(x), as.integer(N/2^(j-1)), L, h, g,
                W=as.double(W), V=as.double(V))
        y[[jj + 2*n]] <- z$W
        y[[jj + 2*n + 1 ]] <- z$V
      }
    }
  }
  attr(y, "wavelet") <- wf
  return(y)
}

idwpt <- function(y, y.basis)
{
  J <- trunc(log(length(y), 2))

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  for(j in J:1) {
    a <- min((1:length(rep(1:J, 2^(1:J))))[rep(1:J, 2^(1:J)) == j])
    b <- max((1:length(rep(1:J, 2^(1:J))))[rep(1:J, 2^(1:J)) == j])
    n <- a
    while(n <= b) {
      if(y.basis[n]) {
        m <- length(y[[n]])
        XX <- numeric(2 * m)
        if(floor((n-a)/2) %% 2 == 0)
          X <- .C(C_idwt, as.double(y[[n+1]]), as.double(y[[n]]),
                  as.integer(m), L, h, g, out=as.double(XX))$out
        else
          X <- .C(C_idwt, as.double(y[[n]]), as.double(y[[n+1]]), 
                  as.integer(m), L, h, g, out=as.double(XX))$out
        if(j != 1) {
          y[[a-(b-a+1)/2 + (n-a)/2]] <- X
          y.basis[[a-(b-a+1)/2 + (n-a)/2]] <- 1
        }
        n <- n + 2
      }
      else { n <- n + 1 }
    }
  }
  return(X)
}

##plot.dwpt <- function(x, n.levels, pgrid=TRUE)
##{
##  J <- n.levels
##  scales <- rep(1:J, 2^(1:J))
##  y <- matrix(NA, 2*length(x[[1]]), J)
##  for(j in 1:J) {
##    a <- min((1:length(scales))[scales == j])
##    b <- max((1:length(scales))[scales == j])
##    y[, j] <- unlist(x[a:b])
##    x.length <- length(y[, j])
##  }
##  plot(ts(y), ylim=c(-.45,.45))
##  if(pgrid) {
##    lines(x.length * c(0,1), c(0,0), lty=2)
##    for(j in 1:J) {
##      lines(x.length * c(0,1), c(-j,-j), lty=2)
##      for(n in 0:2^j) lines(x.length * c(n/2^j, n/2^j), c(-j,-(j-1)), lty=2)
##    }
##  }
##  title(ylab="Level")
##}



#' Produce Boolean Vector from Wavelet Basis Names
#' 
#' Produce a vector of zeros and ones from a vector of basis names.
#' 
#' None.
#' 
#' @param x Output from the discrete wavelet package transfrom (DWPT).
#' @param basis.names Vector of character strings that describe leaves on the
#' DWPT basis tree.  See the examples below for appropriate syntax.
#' @return Vector of zeros and ones.
#' @seealso \code{\link{dwpt}}.
#' @keywords ts
#' @examples
#' 
#' data(acvs.andel8)
#' \dontrun{
#' x <- hosking.sim(1024, acvs.andel8[,2])
#' x.dwpt <- dwpt(x, "la8", 7)
#' ## Select orthonormal basis from wavelet packet tree
#' x.basis <- basis(x.dwpt, c("w1.1","w2.1","w3.0","w4.3","w5.4","w6.10",
#'                            "w7.22","w7.23"))
#' for(i in 1:length(x.dwpt))
#'   x.dwpt[[i]] <- x.basis[i] * x.dwpt[[i]]
#' ## Resonstruct original series using selected orthonormal basis
#' y <- idwpt(x.dwpt, x.basis)
#' par(mfrow=c(2,1), mar=c(5-1,4,4-1,2))
#' plot.ts(x, xlab="", ylab="", main="Original Series")
#' plot.ts(y, xlab="", ylab="", main="Reconstructed Series")
#' }
#' 
#' @export basis
basis <- function(x, basis.names)
{
  m <- length(x)
  n <- length(basis.names)
  y <- numeric(m)
  for(i in 1:n) { y <- y + as.integer(names(x) == basis.names[i]) }
  return(y)
}



#' Derive Orthonormal Basis from Wavelet Packet Tree
#' 
#' An orthonormal basis for the discrete wavelet transform may be characterized
#' via a disjoint partitioning of the frequency axis that covers
#' \eqn{[0,\frac{1}{2})}{[0,1/2)}.  This subroutine produces an orthonormal
#' basis from a full wavelet packet tree.
#' 
#' A wavelet packet tree is a binary tree of Boolean variables.  Parent nodes
#' are removed if any of their children exist.
#' 
#' @param xtree is a vector whose entries are associated with a wavelet packet
#' tree.
#' @return Boolean vector describing the orthonormal basis for the DWPT.
#' @author B. Whitcher
#' @keywords ts
#' @examples
#' 
#' data(japan)
#' J <- 4
#' wf <- "mb8"
#' japan.mra <- mra(log(japan), wf, J, boundary="reflection")
#' japan.nomean <-
#'   ts(apply(matrix(unlist(japan.mra[-(J+1)]), ncol=J, byrow=FALSE), 1, sum),
#'      start=1955, freq=4)
#' japan.nomean2 <- ts(japan.nomean[42:169], start=1965.25, freq=4)
#' plot(japan.nomean2, type="l")
#' japan.dwpt <- dwpt(japan.nomean2, wf, 6)
#' japan.basis <-
#'   ortho.basis(portmanteau.test(japan.dwpt, p=0.01, type="other"))
#' # Not implemented yet
#' # par(mfrow=c(1,1))
#' # plot.basis(japan.basis)
#' 
#' @export ortho.basis
ortho.basis <- function(xtree) {
  J <- trunc(log(length(xtree), 2))
  X <- vector("list", J)
  X[[1]] <- xtree[rep(1:J, 2^(1:J)) == 1]
  for(i in 2:J) {
    for(j in i:J) {
      if(i == 2) X[[j]] <- xtree[rep(1:J, 2^(1:J)) == j]
        X[[j]] <- X[[j]] + 2 * c(apply(matrix(xtree[rep(1:J, 2^(1:J)) == i-1]),
                                       1, rep, 2^(j-i+1)))
    }
  }
  X[[J]][X[[J]] == 0] <- 1
  ifelse(unlist(X) == 1, 1, 0)
}

##plot.basis <- function(xtree)
##{
##  J <- trunc(log(length(xtree), base=2))
##  j <- rep(1:J, 2^(1:J))
##  n <- unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0))
##  basis <- ifelse(xtree, paste("w", j, ".", n, sep=""), NA)
##  pgrid.plot(basis[basis != "NA"])
##  invisible()
##}

phase.shift.packet <- function(z, wf, inv=FALSE)
{
  ## Center of energy
  coe <- function(g)
    sum(0:(length(g)-1) * g^2) / sum(g^2)

  J <- length(z) - 1
  g <- wave.filter(wf)$lpf
  h <- wave.filter(wf)$hpf

  if(!inv) {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(z[[j]])
      z[[j]] <- c(z[[j]][(ph+1):Nj], z[[j]][1:ph])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    z[[J]] <- c(z[[J]][(ph+1):Nj], z[[J]][1:ph])
  } else {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(z[[j]])
      z[[j]] <- c(z[[j]][(Nj-ph+1):Nj], z[[j]][1:(Nj-ph)])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    z[[J]] <- c(z[[j]][(Nj-ph+1):Nj], z[[j]][1:(Nj-ph)])
  }
  return(z)
}

modwpt <- function(x, wf="la8", n.levels=4, boundary="periodic")
{
  N <- length(x); storage.mode(N) <- "integer"
  J <- n.levels
  if(2^J > N) stop("wavelet transform exceeds sample size in modwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  ht <- dict$hpf/sqrt(2)
  storage.mode(ht) <- "double"
  gt <- dict$lpf/sqrt(2)
  storage.mode(gt) <- "double"

  y <- vector("list", sum(2^(1:J)))
  yn <- length(y)
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  names(y) <- paste("w", crystals1, ".", crystals2, sep="")

  W <- V <- numeric(N)
  storage.mode(W) <- storage.mode(V) <- "double"
  for(j in 1:J) {
    index <- 0
    jj <- min((1:yn)[crystals1 == j])
    for(n in 0:(2^j / 2 - 1)) {
      index <- index + 1
      if(j > 1)
        x <- y[[(1:yn)[crystals1 == j-1][index]]]
      if(n %% 2 == 0) {
        z <- .C(C_modwt, as.double(x), N, as.integer(j), L, ht, gt, 
                W = W, V = V)[7:8]
        y[[jj + 2*n + 1]] <- z$W
        y[[jj + 2*n]] <- z$V
      }
      else {
        z <- .C(C_modwt, as.double(x), N, as.integer(j), L, ht, gt, 
                W = W, V = V)[7:8]
        y[[jj + 2*n]] <- z$W
        y[[jj + 2*n + 1 ]] <- z$V
      }
    }
  }
  attr(y, "wavelet") <- wf
  return(y)
}

dwpt.brick.wall <- function(x, wf, n.levels, method="modwpt")
{
  N <- length(x[[1]])
  m <- wave.filter(wf)$length
  J <- n.levels
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))

  if(method=="dwpt") {
    ## for DWPT
    for(j in 1:J) {
      jj <- min((1:length(crystals1))[crystals1 == j])
      L <- switch(j,
                  (m-2)/2,
                  ((m-2)/2 + floor(m/4)),
                  ((m-2)/2 + floor((m/2 + floor(m/4))/2)))
      if(is.null(L)) L <- (m-2)
      for(n in 0:(2^j-1))
        x[[jj+n]][1:L] <- NA
    }
  }
  else {
    ## for MODWPT
    for(j in 1:J) {
      jj <- min((1:length(crystals1))[crystals1 == j])
      L <- min((2^j - 1) * (m - 1), N)
      for(n in 0:(2^j-1))
        x[[jj+n]][1:L] <- NA
    }
  }
  return(x)
}



#' Testing the Wavelet Packet Tree for White Noise
#' 
#' A wavelet packet tree, from the discrete wavelet packet transform (DWPT), is
#' tested node-by-node for white noise.  This is the first step in selecting an
#' orthonormal basis for the DWPT.
#' 
#' Top-down recursive testing of the wavelet packet tree is
#' 
#' @usage cpgram.test(y, p = 0.05, taper = 0.1)
#' @usage css.test(y)
#' @usage entropy.test(y)
#' @usage portmanteau.test(y, p = 0.05, type = "Box-Pierce")
#' @aliases cpgram.test css.test entropy.test portmanteau.test
#' @param y wavelet packet tree (from the DWPT)
#' @param p significance level
#' @param taper weight of cosine bell taper (\code{cpgram.test} only)
#' @param type \code{"Box-Pierce"} and \code{other} recognized
#' (\code{portmanteau.test} only)
#' @return Boolean vector of the same length as the number of nodes in the
#' wavelet packet tree.
#' @author B. Whitcher
#' @seealso \code{\link{ortho.basis}}.
#' @references Brockwell and Davis (1991) \emph{Time Series: Theory and
#' Methods}, (2nd. edition), Springer-Verlag.
#' 
#' Brown, Durbin and Evans (1975) Techniques for testing the constancy of
#' regression relationships over time, \emph{Journal of the Royal Statistical
#' Society B}, \bold{37}, 149-163.
#' 
#' Percival, D. B., and A. T. Walden (1993) \emph{Spectral Analysis for
#' Physical Applications: Multitaper and Conventional Univariate Techniques},
#' Cambridge University Press.
#' @keywords ts
#' @examples
#' 
#' data(mexm)
#' J <- 6
#' wf <- "la8"
#' mexm.dwpt <- dwpt(mexm[-c(1:4)], wf, J)
#' ## Not implemented yet
#' ## plot.dwpt(x.dwpt, J)
#' mexm.dwpt.bw <- dwpt.brick.wall(mexm.dwpt, wf, 6, method="dwpt")
#' mexm.tree <- ortho.basis(portmanteau.test(mexm.dwpt.bw, p=0.025))
#' ## Not implemented yet
#' ## plot.basis(mexm.tree)
#'
css.test <- function(y) 
{
  K <- length(y)
  test <- numeric(K)

  for(k in 1:K) {
    x <- y[[k]]
    x <- x[!is.na(x)]
    n <- length(x)
    plus <- 1:n/(n - 1) - cumsum(x^2)/sum(x^2)
    minus <- cumsum(x^2)/sum(x^2) - 0:(n - 1)/(n - 1)
    D <- max(abs(plus), abs(minus))
    if(D < 1.224/(sqrt(n) + 0.12 + 0.11/sqrt(n))) test[k] <- 1
  }
  return(test)
}

entropy.test <- function(y)
{
  K <- length(y)
  test <- numeric(K)

  for(k in 1:K) {
    x <- y[[k]]
    test[k] <- sum(x^2 * log(x^2), na.rm=TRUE)
  }
  return(test)
}

cpgram.test <- function(y, p=0.05, taper=0.1)
{
  K <- length(y)
  test <- numeric(K)
  
  for(k in 1:K) {
    x <- y[[k]]
    x <- x[!is.na(x)]
    x <- spec.taper(scale(x, center=TRUE, scale=FALSE), p=taper)
    y <- Mod(fft(x))^2/length(x)
    y[1] <- 0
    n <- length(x)
    x <- (0:(n/2))/n
    if(length(x) %% 2 == 0) {
      n <- length(x) - 1
      y <- y[1:n]
      x <- x[1:n]
    }
    else y <- y[1:length(x)]
    mp <- length(x) - 1
    if(p == 0.05)
      crit <- 1.358/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
    else {
      if(p == 0.01) crit <- 1.628/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
      else stop("critical value is not known")
    }
    D <- abs(cumsum(y)/sum(y) -  0:mp/mp)
    if(max(D) < crit) test[k] <- 1
  }
  return(test)
}

portmanteau.test <- function(y, p = 0.05, type = "Box-Pierce")
{
  K <- length(y)
  test <- numeric(K)

  for(k in 1:K) {
  x <- y[[k]]
  x <- x[!is.na(x)]
  n <- length(x)
  h <- trunc(n/2)
  x.acf <- my.acf(x)[1:(h+1)]
  x.acf <- x.acf / x.acf[1];
  if(type == "Box-Pierce")
    test[k] <- ifelse(n * sum((x.acf[-1])^2) > qchisq(1-p, h), 0, 1)
  else 
    test[k] <- ifelse(n*(n+2) * sum((x.acf[-1])^2 / (n - h:1)) > 
                      qchisq(1-p, h), 0, 1)
  }
  return(test)
}

