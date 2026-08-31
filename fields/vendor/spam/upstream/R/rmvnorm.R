# HEADER ####################################################
# This is file spam/R/rmvnorm.R.                            #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



# Draw from a multivariate normal:
# (Algorithm 2.3 from Rue and Held, 2005)
rmvnorm <- function(n, mu = rep.int(0, dim(Sigma)[1]), Sigma, ..., mean, sigma) {

  if (!missing(mean))   mu <- mean
  if (!missing(sigma))  Sigma <- sigma

  if (is( Sigma, "spam")) return( rmvnorm.spam( n, mu, Sigma, ...))

  R <- chol( Sigma, ...)
  N <- dim( Sigma)[1]

  return(sweep( as.matrix( ( array(rnorm(n*N),c(n,N)) %*% R)), 2, mu, "+"))
}

rmvnorm.spam <- function(n, mu = rep.int(0, dim(Sigma)[1]), Sigma, Rstruct = NULL, ..., mean, sigma) {
  # taken from ?chol.spam
  if (!missing(mean))   mu <- mean
  if (!missing(sigma))  Sigma <- sigma

  if (!is(Sigma, "spam"))  return(  rmvnorm(n, mu, Sigma, ...))  # just in case!

  if (is(Rstruct,"spam.chol.NgPeyton"))
    cholS <- update.spam.chol.NgPeyton( Rstruct, Sigma, ...)
  else
    cholS <- chol.spam( Sigma,...)
  # cholS is the upper triangular part of the permutated matrix Sigma
  iord <- ordering(cholS, inv=TRUE)

  N <- dim(Sigma)[1]

  R <- as.spam(cholS)
  retval <- as.matrix( ( array(rnorm(n*N),c(n,N)) %*% R)[,iord,drop=F])
     # It is often better to order the sample than the matrix
     # R itself.
  return(sweep(retval, 2, mu, "+"))
}


rmvt <- function(n, Sigma, df = 1, delta = rep(0, nrow(Sigma)),
                 type = c("shifted", "Kshirsagar"), ..., sigma) {

  if (!missing(sigma))  Sigma <- sigma

  return( rmvt.spam(n = n, Sigma, df, delta, type, ...))
}


rmvt.spam <- function (n, Sigma, df = 1, delta = rep(0, nrow(Sigma)),
    type = c("shifted", "Kshirsagar"), ..., sigma) {

    if (!missing(sigma))  Sigma <- sigma
    if (length(delta) != nrow(Sigma))
        stop("'delta' and 'Sigma' have non-conforming size")
    if (hasArg(mean))
        stop("Providing 'mean' does *not* sample from a multivariate t-distribution!")
    if (df == 0 || ( (df > 0) & is.infinite(df)))
        return(rmvnorm.spam(n, mean = delta, Sigma = Sigma, ...))
    type <- match.arg(type)
    switch(type, Kshirsagar = {
        return(rmvnorm.spam(n, mean = delta, Sigma = Sigma, ...)/sqrt(rchisq(n,
            df)/df))
    }, shifted = {
        return( sweep( rmvnorm.spam(n, Sigma = Sigma, ...)/sqrt(rchisq(n,
            df)/df), 2, delta, "+"))
    }, stop("wrong 'type'"))
}






# Draw from a multivariate normal given a precision matrix:
# (Algorithm 2.4 from Rue and Held, 2005)
rmvnorm.prec <- function(n, mu = rep.int(0, dim(Q)[1]), Q, Rstruct = NULL, ...) {

  if (is(Rstruct,"spam.chol.NgPeyton"))
    R <- update.spam.chol.NgPeyton( Rstruct, Q, ...)
  else
    R <- chol(Q, ...)
  # R is the upper triangular part of the permutated matrix Sigma

  N <- dim(Q)[1]
  nu <- backsolve(R, array(rnorm(n*N), c(N, n)), k = N)
  return(t(nu + mu))
}


# Draw from the canonical representation of a multivariate normal:
# (Algorithm 2.5 from Rue and Held, 2005)
rmvnorm.canonical <- function(n, b, Q, Rstruct=NULL, ...) {
  N=dim(Q)[1]
  if (is(Rstruct,"spam.chol.NgPeyton"))
    R <- update.spam.chol.NgPeyton( Rstruct, Q, ...)
  else
    R <- chol(Q,...)

  if(is(R,"spam.chol.NgPeyton")){
     mu <- drop(solve.spam( R, b))
  } else {
     mu <- backsolve( R, forwardsolve( t(R), b), k=N)
  }
  nu <- backsolve(R, array( rnorm(n*N), c(N, n)), k=N)
  return(t(nu+mu))
}



rmvnorm.const <- function (n, mu = rep.int(0, dim(Sigma)[1]), Sigma,
                                Rstruct = NULL,
                                A = array(1, c(1,dim(Sigma)[1])), a=0, U=NULL,  ...)
{
    N <- dim(Sigma)[1]
    if (!identical(dim(A)[2], N)) stop("Incorrect constraint specification")

    if (is(Rstruct, "spam.chol.NgPeyton"))
        cholS <- update.spam.chol.NgPeyton(Rstruct, Sigma, ...)
    else cholS <- chol.spam(Sigma, ...)
    iord <- ordering(cholS, inv = TRUE)
    N <- dim(Sigma)[1]
    R <- as.spam(cholS)
    x <- sweep( (array(rnorm(n * N), c(n, N)) %*% R)[, iord, drop = F], 2, mu, "+")

    if (is.null(U)){
      V <- backsolve( R, forwardsolve( R, t(A)), k=N)
      W <- A %*% V
      U <- solve(W, t(V))
    }
    correct <- A %*% t(x) - a
    return(x - t( t(U)%*% correct))
}


rmvnorm.prec.const <- function (n, mu = rep.int(0, dim(Q)[1]), Q,
                                Rstruct = NULL,
                                A = array(1, c(1,dim(Q)[1])), a=0, U=NULL,  ...)
{
    N = dim(Q)[1]
    if (!identical(dim(A)[2], N)) stop("Incorrect constraint specification")

    if (is(Rstruct, "spam.chol.NgPeyton"))
        R <- update.spam.chol.NgPeyton(Rstruct, Q, ...)
    else R <- chol(Q, ...)
    x <- backsolve(R, array(rnorm(n * N), c(N, n)), k=N) + mu


    if (is.null(U)){
        tV <- t( backsolve( R, forwardsolve( R, t(A)), k=N))
        W <- tcrossprod(A, tV)
        U <- solve(W, tV)
    }
    correct <- A %*% x - a
    return(t(x- t(U) %*% correct))
}



rmvnorm.canonical.const <- function (n, b, Q, Rstruct = NULL,
                                     A = array(1, c(1,dim(Q)[1])), a=0, U=NULL, ...)
{
    N = dim(Q)[1]
    if (!identical(dim(A)[2], N)) stop("Incorrect constraint specification")

    if (is(Rstruct, "spam.chol.NgPeyton"))
        R <- update.spam.chol.NgPeyton(Rstruct, Q, ...)
    else R <- chol(Q, ...)
    if (is(R, "spam.chol.NgPeyton")) {
        mu <- drop(solve.spam(R, b))
    }
    else {
        mu <- backsolve(R, forwardsolve(t(R), b))
    }
    x <- backsolve(R, array(rnorm(n * N), c(N, n)), k=N) + mu


    if (is.null(U)){
        tV <- t( backsolve( R, forwardsolve( R, t(A)), k=N))
        W <- tcrossprod(A, tV)
        U <- solve(W, tV)
    }
    correct <- A %*% x - a
    return(t(x- t(U) %*% correct))


}



### conditional simulation
rmvnorm.conditional <- function(n, y, mu = rep.int(0, dim(SigmaXX)[1]+dim(SigmaYY)[1]),
                                   SigmaXX, SigmaYY, SigmaXY,
                                   noise,
                                   RstructYY = NULL, ...) {

  if (!is(SigmaXX, "spam"))  stop("Implemented only for 'spam' covariance matrices")

  nx <- dim(SigmaXX)[1]
  if (missing(SigmaYY)) {
    if (missing(noise)) stop("Either 'noise' or 'SigmaYY' are required")
    SigmaYY <- SigmaXX+diag.spam(noise, nx, nx)
  } else if (!missing(noise)) stop("Either 'noise' or 'SigmaYY' are required, not both")

  ny <- dim(SigmaYY)[1]
  if (length(y)!=ny)  stop("Length of 'y' does not match dimension of 'SigmaYY'")

  if (missing(SigmaXY)) {  # save one transpose
    if (nx!=ny) stop("'SigmaXY' needs to be specified")
    Sigma <- rbind.spam(cbind.spam(SigmaXX, SigmaXX), cbind(SigmaXX, SigmaYY))
    SigmaXY <- SigmaXX   # alternative would be to have another if below...
  } else {
    if (!all(dim(SigmaXY) == c(nx,ny))) stop("Wrong dimension for 'SigmaXY'")
    Sigma <- rbind.spam(cbind.spam(SigmaXX, SigmaXY), cbind(t(SigmaXY), SigmaYY))
  }
  #  Rstruct <- chol.spam(Sigma)
  rsample <- t(rmvnorm.spam(n = n, mu = mu, Sigma = Sigma)) #, Rstruct = Rstruct)

  t(rsample[1:nx,] + SigmaXY %*%
      solve.spam(a = SigmaYY, b = y-rsample[(nx+1):(nx+ny),],
                 Rstruct = RstructYY))
}
rmvnorm.cond <- rmvnorm.conditional
