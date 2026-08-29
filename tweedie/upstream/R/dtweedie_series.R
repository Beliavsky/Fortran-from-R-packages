#' @title Series Evaluation for the Tweedie Probability Function
#' @name dtweedie_series
#' @description
#' Evaluates the probability density function (\acronym{pdf}) for Tweedie distributions using an infinite series, 
#' for given values of the dependent variable \code{y}, the mean \code{mu}, dispersion \code{phi}, and power parameter \code{power}.
#' \emph{Not usually called by general users}, but can be used in the case of evaluation problems.
#'
#' @usage dtweedie_series(y, power, mu, phi, details = FALSE)
#'
#' @param y vector of quantiles.
#' @param power scalar; the value of \eqn{p}{power} such that the variance is \eqn{\mbox{var}[Y]=\phi\mu^{p}}{var[Y] = phi * mu^power}.
#' @param mu vector of mean \eqn{\mu}{mu}.
#' @param phi vector of dispersion parameters \eqn{\phi}{phi}.
#' @param details logical; if \code{TRUE}, returns the value of the distribution function and some details.
#' 
#' @return A numeric vector of densities.
#' 
#' @references
#' Dunn, Peter K and Smyth, Gordon K (2005).
#' Series evaluation of Tweedie exponential dispersion model densities
#' \emph{Statistics and Computing},
#' \bold{15}(4). 267--280.
#' \doi{10.1007/s11222-005-4070-y}
#' 

#' @examples 
#' # Plot a Tweedie density
#' y <- seq(0.01, 4, length = 50)
#' fy <- dtweedie_series(y, power = 1.1, mu = 1, phi = 1)
#' plot(y, fy, type = "l", lwd = 2, ylab = "Density")

#' @importFrom stats dgamma dpois 
#' 
#' @keywords distribution
#'
#' @export
dtweedie_series <- function(y, power, mu, phi, details = FALSE){ 
  # Evaluates the Tweedie density using a series expansion
  
  if ( power < 1) stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y < 0) ) stop("y must be a non-negative vector.")
  if ( any(mu <= 0) ) stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  y0 <- (y == 0)
  yp <- (y != 0)
  density <- lo <- hi <-array( dim = length(y))
  
  if ( (power == 2) | (power == 1) ) { # Special cases
    if ( power == 2 ){
      density <- dgamma( y, 
                         shape = 1 / phi, 
                         rate = 1 / (phi * mu ) )
    }
    
    if ( (power == 1) ){
      density <- dpois(x = y / phi, 
                       lambda = mu / phi ) / phi 
        # Using identity: f(y; mu, phi) = c f(cy; c mu, c phi) for p=1.  Now set c = 1/phi
      if ( !all(phi == 1)){
        warnings("The density computations when phi=1 may be subject to errors using floating-point arithmetic\n")
      }
    }
    # CHANGED: 31 October, thanks to Dina Farkas' email 18 October 2012.  WAS:
    #   if ( (power == 1) & (all(phi==1)) ){
    #    density <- dpois(x=y/phi, lambda=mu/phi )
    #   }
  } else{
    
    if ( any(y == 0) ) {
      if (power > 2) {
        density[y0] <- 0 * y[y0]
      }
      if ( (power > 1) && (power < 2) ) {
        lambda <- mu[y0] ^ (2 - power) / ( phi[y0] * (2 - power) )
        density[y0] <- exp( -lambda )
      }
    }
    
    if ( any( y!= 0 ) ) { 
      if (power > 2) {
        out <- dtweedie_series_bigp(power = power,
                                    mu = mu[yp], 
                                    y = y[yp],
                                    phi = phi[yp])
        density[yp] <- out$density 
        lo[yp] <- out$lo
        hi[yp] <- out$hi
      }
      
      if ( ( power > 1 ) && ( power < 2 ) ) {
        out <- dtweedie_series_smallp(power = power,
                                      mu = mu[yp], 
                                      y = y[yp],
                                      phi = phi[yp])
        density[yp] <- out$density 
        lo[yp] <- out$lo
        hi[yp] <- out$hi
      }
    }
  }
  
  if (details) {
    return(list( density = density,
                 lo = lo,
                 hi = hi) )
  } else {
    return(density)
  }
  
}




#' @rdname dtweedie_series
#' @export
dtweedie.series <- function(y, power, mu, phi){ 
  lifecycle::deprecate_warn(when = "3.0.5", 
                            what = "dtweedie.series()", 
                            with = "dtweedie_series()")
  dtweedie_series(y, power, mu, phi)
}




#############################################################################


#' @noRd
dtweedie_series_smallp <- function(power, y, mu, phi){ 
  
  #
  # Peter K Dunn
  # 02 Feb 2000
  #
  
  #
  # Error traps
  #
  
  if ( power < 1) stop("power must be between 1 and 2.")
  if ( power > 2) stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) & (power >= 2) ) stop("y must be a strictly positive vector.")
  if ( any(y < 0) & (power > 1) & (power < 2) ) stop("y must be a non-negative vector.")
  if ( any(mu <= 0) ) stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
  }
  
  result <- dtweedie_logw_smallp( y = y, 
                                  power = power, 
                                  phi = phi)
  logw <- result$logw
  
  
  tau <- phi * (power - 1) * mu ^ ( power - 1 )
  lambda <- mu ^ (2 - power) / ( phi * (2 - power) )
  logf <- -y / tau - lambda - log(y) + logw
  f <- exp( logf )
  
  list(density = f, 
       logw = logw, 
       hi = result$hi, 
       lo = result$lo)
  
}


#############################################################################


#' @noRd
dtweedie_jw_smallp <- function(y, phi, power){ 
  #
  # Peter K Dunn
  # 18 Jun 2002
  #
  
  #
  # Error traps
  #
  
  if ( power < 1) stop("power must be between 1 and 2.")
  if ( power > 2) stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be strictly positive.")
  if ( any(y <= 0) ) stop("y must be a strictly positive vector.")
  
  #
  # Set up
  #
  p <- power
  a <- (2 - p) / (1 - p)         # Note that a<0 for 1<p<2
  
  a1 <- 1 - a
  r <- -a * log(y) + a * log(p - 1) - a1 * log(phi) - log(2 - p)        # All terms to power j
  
  drop <- 37             # Accuracy of terms: exp(-37)
  #
  # Find limits of summation using Stirling's approximation
  # to approximate gamma terms, and find max as j.max
  #
  logz <- max(r)             # To find largest  j  needed
  j.max <- max( y^(2 - p) / ( phi * (2 - p) ) )
  j <- max( 1, j.max )
  
  c <- logz + a1 + a * log(-a)
  wmax <- a1 * j.max
  estlogw <- wmax
  
  # First, the upper limit of j
  
  while(estlogw > (wmax - drop) ){
    j <- j + 2
    estlogw <- j * (c - a1 * log(j))
  }
  
  hi.j <- ceiling(j)
  
  # Now the lower limit of j
  logz <- min(r) 
  j.max <- min( y ^ ( 2 - power ) / ( phi * (2 - power) ) )
  
  j <- max( 1, j.max)
  wmax <- a1 * j.max 
  estlogw <- wmax 
  
  while ( ( estlogw > (wmax - drop) ) && ( j >= 2) ) {
    j <- max(1, j - 2)
    estlogw <- j*(c - a1 * log(j))
  }
  
  lo.j <- max(1, floor(j))
  # 
  # Now sum the series between established limits.
  # We ensure it works for vector y.
  
  j <- seq( lo.j, hi.j)                 # sequence of j terms needed
  o <- matrix( 1, nrow = length(y))     # matrix of ones
  
  g <- matrix(lgamma( j + 1 ) + lgamma( -a * j ),
              nrow = 1, 
              ncol = hi.j - lo.j + 1) 
  logj <- matrix(log(j),
                 nrow = 1, 
                 ncol = hi.j - lo.j + 1) 
  og <- o %*% g                             # matrix of gamma terms
  ologj <- o %*% logj 
  A <- outer(r, j) - og + ologj        # the series, almost ready to sum
  m <- apply(A, 1, max)                    # avoid overflow; find maximum values
  we <- exp( A - m )                    # evaluate terms, less max.
  sum.we <- apply( we, 1, sum)            # sum terms
  jw <- sum.we * exp( m )                # now restore max.
  # Since derivs may be negative, can't use log-scale
  
  list(lo = lo.j, 
       hi = hi.j, 
       jw = jw, 
       j.max = j.max )
  
}



#############################################################################


#' @noRd
dtweedie_kv_bigp <- function(y, phi, power){ 
  # 
  # Peter K Dunn 
  # 18 Jun 2002
  # 
  
  #
  # Error traps
  #
  
  if ( power < 2) stop("power must be greater than 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) ) stop("y must be a strictly positive vector.")
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  
  #
  # Set up
  #
  
  p <- power
  a <- (2 - p) / (1 - p)
  a1 <- 1 - a
  
  r <- -a1 * log(phi) - log(p - 2) - a * log(y) + a * log(p - 1)
  drop <- 37
  
  #
  # Now we find limits of summation, using Stirling's approximation
  # to approximate the gamma terms, and find max as k.max.
  #
  
  logz <- max(r)
  k.max <- max( y ^ (2 - p) / ( phi * (p - 2) ) )
  k <- max( 1, k.max )
  
  c <- logz + a1 + a * log(a)
  vmax <- k.max * a1
  estlogv <- vmax
  #
  # Now we search either side for when we can
  # ignore terms as they are negligible.
  # Remember we are using the envelope as the tests,
  # not individual terms in the series.
  
  # First:  the upper limit of k
  
  while ( estlogv > (vmax - drop) ) {
    k <- k + 2
    estlogv <- k * ( c - a1 * log(k) )
  }
  
  hi.k <- ceiling(k)
  #
  # Now the lower limit of k
  #
  logz <- min(r)
  k.max <- min( y ^ (2 - p) / ( phi * (p - 2) ) ) 
  k <- max( 1, k.max )
  c <- logz + a1 + a * log(a)
  vmax <- k.max * a1
  estlogv <- vmax
  
  while ( (estlogv > (vmax - drop) ) && ( k >= 2) ) {
    k <- max(1, k - 2)
    estlogv <- k * ( c - a1 * log(k) )
  }
  
  lo.k <- max(1, floor(k) )
  
  #
  # Now sum the series  between established limits.
  # We ensure it works for vector y.
  #
  k <- seq(lo.k, hi.k) 
  o <- matrix( 1, nrow = length(y)) 
  
  g <- matrix( lgamma( 1 + a * k) - lgamma(1 + k),
               nrow = 1, 
               ncol  =length(k) )
  logk <- matrix( log(k),
                  nrow = 1, 
                  ncol = length(k) )
  
  og <- o %*% g
  ologk <- o %*% logk
  
  A <- outer(r, k) + og + ologk
  
  C <- matrix( sin( -a * pi * k ) * (-1) ^ k,
               nrow = 1, 
               ncol = length(k) )
  C <- o %*% C
  
  m <- apply(A, 1, max)
  ve <- exp(A - m)
  sum.ve <- apply( ve*C, 1, sum )
  kv <- sum.ve * exp( m )
  # Since derivs may be negative, can't use log-scale
  
  list(lo = lo.k, 
       hi = hi.k, 
       kv = kv, 
       k.max = k.max )
  
}
#############################################################################

#' @noRd
dtweedie_logv_bigp <- function( y, phi, power){ 
  # Peter K Dunn 
  # 02 Feb 2000 
  # 
  
  #
  # Error traps
  #
  
  if ( power < 2) stop("power must be greater than 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) ) stop("y must be a strictly positive vector.")
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  #
  # Set up
  #
  
  p <- power
  a <- (2 - p) / (1 - p)
  a1 <- 1 - a
  
  r <- -a1 * log(phi) - log(p - 2) - a * log(y) + a * log(p - 1)
  drop <- 37
  
  #
  # Now we find limits of summation, using Stirling's approximation
  # to approximate the gamma terms, and find max as k.max.
  #
  
  logz <- max(r)
  k.max <- max( y ^ (2 - p) / ( phi * (p - 2) ) )
  k <- max( 1, k.max )
  
  c <- logz + a1 + a * log(a)
  vmax <- k.max * a1
  estlogv <- vmax
  #
  # Now we search either side for when we can
  # ignore terms as they are negligible.
  # Remember we are using the envelope as the tests,
  # not individual terms in the series.
  
  # First:  the upper limit of k
  
  while ( estlogv > (vmax - drop) ) {
    k <- k + 2
    estlogv <- k*( c - a1 * log(k) )
  }
  
  hi.k <- ceiling(k)
  #
  # Now the lower limit of k
  #
  logz <- min(r)
  k.max <- min( y ^ (2 - p) / ( phi * (p - 2) ) )
  k <- max( 1, k.max )
  c <- logz + a1 + a * log(a)
  vmax <- k.max * a1
  estlogv <- vmax
  
  while ( (estlogv > (vmax - drop) ) && ( k >= 2) ) {
    k <- max(1, k - 2)
    estlogv <- k*( c - a1 * log(k) )
  }
  
  lo.k <- max(1, floor(k) )
  
  #
  # Now sum the series between established limits.
  # We ensure it works for vector y.
  #
  k <- seq(lo.k, hi.k) 
  o <- matrix( 1, nrow = length(y)) 
  
  g <- matrix( lgamma( 1 + a * k) - lgamma(1 + k),
               nrow = 1, 
               ncol = length(k) )
  
  og <- o %*% g
  A <- outer(r, k) + og
  C <- matrix( sin( -a * pi * k ) * (-1)^k,
               nrow = 1, 
               ncol = length(k) )
  
  C <- o %*% C
  m <- apply(A, 1, max)
  ve <- exp(A - m)
  sum.ve <- apply( ve * C, 1, sum )
  
  # Now be careful!  Because of the +/- nature of the sin term,
  # sum.ve can be very small but negative  due to subtractive
  # cancellation.  Treat those carefully and separately.
  
  neg.sum.ve <- (sum.ve <= 0)
  pos.sum.ve <- (sum.ve > 0)
  
  logv <- sum.ve
  sum.ve[neg.sum.ve] <- 0
  logv[neg.sum.ve] <- -Inf
  logv[pos.sum.ve] <- log( sum.ve[pos.sum.ve] ) + m[pos.sum.ve]
  
  list(lo = lo.k, 
       hi = hi.k, 
       logv = logv, 
       k.max = k.max )
  
}


#############################################################################

#' @noRd
dtweedie_logw_smallp <- function(y, phi, power){ 
  #
  # Peter K Dunn
  # 02 Feb 2000
  #
  
  #
  # Error traps
  #
  
  if ( power < 1) stop("power must be between 1 and 2.")
  if ( power > 2) stop("power must be between 1 and 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) ) stop("y must be a strictly positive vector.")
  
  #
  # Set up
  #
  p <- power
  a <- (2 - p) / (1 - p)         # Note that a<0 for 1<p<2
  
  a1 <- 1 - a
  r <- -a * log(y) + a * log(p - 1) - a1 * log(phi) - log(2 - p)        # All terms to power j
  
  drop <- 37             # Accuracy of terms: exp(-37)
  #
  # Find limits of summation using Stirling's approximation
  # to approximate gamma terms, and find max as j.max
  #
  logz <- max(r)             # To find largest  j  needed
  j.max <- max( y ^ (2 - p) / ( phi * (2 - p) ) )
  j <- max( 1, j.max )
  
  cc <- logz + a1 + a * log(-a)    #
  # This is all the terms to j-power, not needing any other j terms.
  # The other terms are introduced when we know the values of j
  
  wmax <- a1 * j.max
  estlogw <- wmax
  
  # First, the upper limit of j
  while(estlogw > (wmax - drop) ){
    j <- j + 2
    estlogw <- j * (cc - a1 * log(j))
  }
  
  hi.j <- ceiling(j)
  
  #
  #
  # Now the lower limit of j
  #
  #
  logz <- min(r) 
  j.max <- min( y ^ (2 - power) / ( phi * (2 - power) ) )
  
  j <- max( 1, j.max)
  wmax <- a1 * j.max 
  estlogw <- wmax 
  
  # First, optimize to find the location of the maximum
  while ( ( estlogw > (wmax - drop) ) && ( j >= 2) ) {
    j <- max(1, j - 2)
    oldestlogw <- estlogw
    estlogw <- j * (cc - a1 * log(j))
  }
  
  lo.j <- max(1, floor(j))
  # 
  # Now sum the series between established limits.
  # We ensure it works for vector y.
  
  j <- seq( lo.j, hi.j)      # sequence of j terms needed
  
  o <- matrix( 1, nrow = length(y))     # matrix of ones
  
  g <- matrix(lgamma( j + 1 ) + lgamma( -a * j ),
              nrow = 1, 
              ncol = hi.j - lo.j + 1)
  
  og <- o %*% g                # matrix of gamma terms
  A <- outer(r, j) - og        # the series, almost ready to sum
  m <- apply(A, 1, max)        # avoid overflow; find maximum values
  we <- exp( A - m )           # evaluate terms, less max.
  sum.we <- apply( we, 1, sum) # sum terms
  logw <- log( sum.we ) + m    # now restore max.
  
  list(lo = lo.j, 
       hi = hi.j, 
       logw = logw, 
       j.max = j.max )
  
}



#############################################################################

#' @noRd
dtweedie_series_bigp <- function(power, y, mu, phi){ 
  
  # 
  # Peter K Dunn 
  # 02 Feb 2000 
  # 
  
  #
  # Error traps
  #
  
  if ( power < 2) stop("power must be greater than 2.")
  if ( any(phi <= 0) ) stop("phi must be positive.")
  if ( any(y <= 0) ) stop("y must be a strictly positive vector.")
  if ( any(mu <= 0) ) stop("mu must be positive.")
  if ( length(mu) > 1) {
    if ( length(mu) != length(y) ) stop("mu must be scalar, or the same length as y.")
  } else {
    mu <- array( dim = length(y), mu )
    # A vector of all mu's
  }
  if ( length(phi) > 1) {
    if ( length(phi) != length(y) ) stop("phi must be scalar, or the same length as y.")
  } else {
    phi <- array( dim = length(y), phi )
    # A vector of all phi's
  }
  
  
  result <- dtweedie_logv_bigp(power = power, 
                               y = y, 
                               phi = phi)
  logv <- result$logv
  
  theta <- mu ^ (1 - power) / ( 1 - power )
  kappa <- mu ^ (2 - power) / ( 2 - power )
  
  logfnew <- (y * theta - kappa) / phi - log( pi * y) + logv
  f <- exp( logfnew )
  
  return( list(density = f,
               logv = logv, 
               lo = result$lo, 
               hi = result$hi ) )
  
}



