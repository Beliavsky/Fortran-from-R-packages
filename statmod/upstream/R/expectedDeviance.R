expectedDeviance <- function(mu, family="binomial", binom.size, nbinom.size, gamma.shape)
# Expectation and variance of the unit deviance for linear exponential families
# Lizhong Chen and Gordon Smyth
# Created 02 October 2022, last revised 28 December 2022.
{
# For simplicity, NA inputs or invalid arguments are not allowed and will generate an error.

# Output will preserve dimensions and attributes of `mu`.
  out <- list(mean=mu,variance=mu)
  m <- as.numeric(mu)
  length.m <- length(m)
  if(!length.m) return(out)

# Check family
  if(identical(family,"Poisson")) family <- "poisson"
  if(identical(family,"gamma")) family <- "Gamma"
  family <- match.arg(family,c("binomial","gaussian","Gamma","inverse.gaussian","poisson","negative.binomial"))

  if(identical(family,"binomial")) {
#   Check binom.size
    binom.size <- as.integer(binom.size)
    length.n <- length(binom.size)
    if(!identical(length.n,1L) && !identical(length.n,length.m))
      stop("binom.size must have length 1 or length must agree with mu")
    min.n <- min(binom.size)
    if(is.na(min.n)) stop("NAs not allowed in binom.size")
    if(min.n < 1) stop("binom.size must be >= 1")

#   Check for permissable mu
    min.m <- min(m)
    max.m <- max(m)
    if(is.na(min.m)) stop("NAs not allowed in mu")
    if(min.m < 0 || max.m > 1) stop("binomial mu must be between 0 and 1")

    big.n <- 200L
    C.out <- .C("mbinomdev", m, binom.size,
      mean = double(length.m), variance = double(length.m),
      length.m, length.n, big.n)[c(3,4)]
    out$mean[] <- C.out$mean
    out$variance[] <- C.out$variance
  }

  if(identical(family,"gaussian")) {
    out$mean[] <- 1
    out$variance[] <- 2
  }

  if(identical(family,"Gamma")) {
    gamma.shape <- as.numeric(gamma.shape)
    length.s <- length(gamma.shape)
    if(!identical(length.s,1L) && !identical(length.s,length.m))
      stop("gamma.shape must have length 1 or length must agree with mu")
    min.s <- min(gamma.shape)
    if(is.na(min.s)) stop("NAs not allowed in gamma.shape")
    if(min.s <=0) stop("gamma.shape should be positive")
    out$mean[] <- meanval.digamma(-gamma.shape) * gamma.shape
    out$variance[] <- 2*d2cumulant.digamma(-gamma.shape) * gamma.shape * gamma.shape
  }

  if(identical(family,"inverse.gaussian")) {
    min.m <- min(m)
    if(is.na(min.m)) stop("NAs not allowed in mu")
    if(min.m < 0) stop("inverse.gaussian mu must be non-negative")
    out$mean[] <- 1
    out$variance[] <- 2
  }

  if(identical(family,"poisson")) {
    min.m <- min(m)
    if(is.na(min.m)) stop("NAs not allowed in mu")
    if(min.m < 0) stop("Poisson mu must be non-negative")
    m <- pmin(m,1e7)
    C.out <- .C("mpoisdev", m,
      mean = double(length.m), variance = double(length.m),
      length.m)[c(2,3)]
    out$mean[] <- C.out$mean
    out$variance[] <- C.out$variance
  }

  if(identical(family,"negative.binomial")) {
#   Check nbinom.size
    nbinom.size <- as.numeric(nbinom.size)
    length.s <- length(nbinom.size)
    if(!identical(length.s,1L) && !identical(length.s,length.m))
      stop("nbinom.size must have length 1 or length must agree with mu")
    min.s <- min(nbinom.size)
    if(is.na(min.s)) stop("NAs not allowed in nbinom.size")
    if(min.s <= 0) stop("nbinom.size must be positive")

#   Large size corresponds to Poisson
    if(min.s > 1e7) return(Recall(mu=mu,family="poisson"))
    nbinom.size <- pmin(nbinom.size,1e7)

#   Need to avoid very small size parameters
#   Chebychev approximation works for size > 1/4. For a limited range of
#   smaller values, direct summation is used.
    if(min.s <= 0.25) {
      limit.size <- pmin(mu*mu/(1e5 - 10),0.25) + 1e-10
      nbinom.size <- pmax(nbinom.size,limit.size)
    }

#   Check for permissable mu
    min.m <- min(m)
    if(is.na(min.m)) stop("NAs not allowed in mu")
    if(min.m < 0) stop("Negative binomial mu must be non-negative")
    m <- pmin(m,1e7)
    
    C.out <- .C("mnbinomdev", m, nbinom.size,
      mean = double(length.m), variance = double(length.m),
      length.m, length.s)[c(3,4)]
    out$mean[] <- C.out$mean
    out$variance[] <- C.out$variance
  }

  out
}
