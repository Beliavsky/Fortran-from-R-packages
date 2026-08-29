
#' @noRd
dtweedie_dldphi_saddle <- function(phi, mu, power, y){
  # Calculates the derivative of log f wrt phi
  # where the density is the saddlepoint density
  
  # Peter Dunn
  # 13 August 2002
  
  dev <- tweedie_dev( power = power, 
                      y = y, 
                      mu = mu)
  l <-  (-1) / (2 * phi) + dev / (2 * phi ^ 2)
  
  -2 * sum(l)
}




#############################################################################

#' @noRd
dtweedie_logl <- function(phi, y, mu, power) {
  # Computes the log-likelihood for
  # a Tweedie density.  
  
  # Peter Dunn
  # 26 April 2001
  sum( log( dtweedie( y = y, 
                      mu = mu, 
                      phi = phi, 
                      power = power) ) )
  
}



#############################################################################

#' @noRd
dtweedie_logl_saddle <- function( phi, power, y, mu, eps=0){
  # Calculates the log likelihood of Tweedie densities
  # where the density is the saddlepoint density
  
  # Peter Dunn
  # 01 May 2001
  sum( log( dtweedie_saddle(power = power, 
                            phi = phi, 
                            y = y, 
                            mu = mu,
                            eps = eps) ) )
  
}



#############################################################################

#' @noRd
dtweedie_dldphi <- function(phi, mu, power, y ){
  # Calculates the log-likelihood
  # function, wrt phi, for p>2.  In particular, it returns
  #    sum{ d(log f)/d(phi) } = d( log-likelihood )/d(phi).
  # The mle of phi can be found, therefore, by setting this to zero.
  # y  is generally a vector of observed values.
  
  #
  # Peter Dunn
  # 31 Jan 2001
  #

  if ( (power != 2 ) & ( power != 1 ) ) {
    
    k <- phi ^ (1 / (power - 2))
    # cat("k=",k,"\n")
    # cat("phi=",phi,"\n")
    # cat("power=",power,"\n")
    if ( k < 1  & k > 0 ) {
      # Use the transform f(y; mu, phi) = c f(c*y; c*mu, c^(2-p)*phi)
      # and differentiate with c=phi^(1/(p-2)):
      #    d log f / d phi = c^(2-p) * {df(cy; c*mu, 1)/dphi} / f(cy; c*mu, 1)
      f <- dtweedie( y = k * y, 
                     power = power, 
                     mu = k * mu, 
                     phi = 1 )
      d <- dtweedie_dlogfdphi( y = k * y, 
                               power = power, 
                               mu = k * mu, 
                               phi = 1 )
      # Note:  We need dlogf/dphi = dlogf.dphi * f
      top <- d * f
      d <- -2* sum( top / f * k ^ (2 - power) )
      
    } else{
      # Compute directly
      d <- -2 * sum( dtweedie_dlogfdphi(y = y, 
                                        power = power, 
                                        mu = mu, 
                                        phi = phi) )
    }
  } else{
    # Cases p == 1 and  p == 2 
    d <- -2 * sum( dtweedie_dlogfdphi(y = y, 
                                      power = power, 
                                      mu = mu, 
                                      phi = phi) )
  }
  d
}



#############################################################################

#' @noRd
dtweedie_dlogfdphi <- function(y, mu, phi, power)
{
  #
  # Calculates d(log f)/d(phi) for the Tweedie
  # densities.
  # It is used, for example, in mle fitting of phi.  We would then
  # sum over  y  and set this function to 0.
  #
  #
  # Peter Dunn
  # 31 Jan 2001
  #
  p <- power
  a <- (2 - p) / (1 - p)
  if(length(phi) == 1) {
    phi <- array(dim = length(y), phi)
  }
  if(length(mu) == 1) {
    mu <- array(dim = length(y), mu)
  }
  A <- (y * mu ^ (1 - p)) / (phi ^ 2 * (p - 1))
  B <- mu^(2 - p)/(phi ^ 2 * (2 - p))
  
  if(power > 2) {
    f <- array(dim = c(length(y)))
    # Here, we evaluate logv and everything as normal.
    # If logv has infinite values then we resort to other tactics.
    
    kv <- dtweedie_kv_bigp(power = power, 
                           phi = phi, 
                           y = y)$kv
    dv.dphi <- (kv * (a - 1)) / phi
    out.logv <- dtweedie_logv_bigp(power = power, 
                                   phi = phi, 
                                   y = y)
    
    # Now see if this causes problems.
    logv <- out.logv$logv
    
    # Now detect problem computing  logv  and remedy them
    probs <- (is.infinite(logv)) | (is.nan(logv)) | (y < 1)
    
    if(any(probs)) {
      
      # OK then:  Troubles computing log(V).
      # best we can do is use definition I think.
      delta <- 1.0e-5
      a1 <- dtweedie(power = power, 
                     phi = phi[probs], 
                     mu = mu[probs], 
                     y = y[probs])
      a2 <- dtweedie(power = power, 
                     phi = phi[probs] + delta, 
                     mu = mu[probs], 
                     y = y[probs])
      f[probs] <- (log(a2) - log(a1) ) / delta
      
    }
    
    f[!probs] <- A[!probs] + B[!probs] + dv.dphi[ !probs] / exp(logv[!probs])
  }
  #  END p>2
  
  if(power == 2) {
    f <-   -log(y)  + ( y / mu ) + digamma(1 / phi) - 1 + log( mu * phi )
    f <- f / (phi ^ 2)
  }
  
  if(power == 1) {
    
    f <- mu - y - y * log(mu / phi) + y * digamma(1 + (y / phi))
    f <- f / (phi ^ 2)
  }
  
  if((power > 1) && (power < 2)) {
    # We need to treat y==0 separately.
    # We fill  f  with the Y=0 case
    f <- array(  dim = length(y), 
                 mu ^ (2 - power) / ( phi ^ 2 * (2 - power) )  )
    jw <- dtweedie_jw_smallp(power = power, 
                             phi = phi[y > 0], 
                             y = y[y > 0])$jw
    dw.dphi <- (jw * (a - 1)) / phi[y > 0]
    logw <- dtweedie_logw_smallp(power = power, 
                                 phi = phi[y > 0], 
                                 y = y[y > 0])$logw
    f[y>0] <- A[y > 0] + B[y > 0] + dw.dphi / exp(logw)
  }
  f
}


