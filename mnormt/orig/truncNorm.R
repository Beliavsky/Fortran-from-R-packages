# Code providing support for the truncated normal and "t" distributions.
# A.Azzalini, 2020-2022.
#--------------------------------------------------------------------------

dmtruncnorm <- function(x, mean, varcov, lower, upper, log= FALSE, ...) {
  d <- if(is.matrix(varcov)) ncol(varcov) else 1
  if(missing(lower)) lower <- rep(-Inf,d)
  if(missing(upper)) upper <- rep(Inf,d)
  if(length(lower) != d | length(upper) != d) stop("dimension mismatch")
  if(!all(lower < upper)) stop("lower<upper is required (componentwise)")         
  if(d==1) {
    p.Int <- diff(pnorm(c(lower, upper), mean, sqrt(varcov)))
    pdf <- rep(0, length(c(x)))
    inside <- (x>lower & x<upper)
    pdf[!inside] <- if(log) -Inf else 0
    pdf.Int <- dnorm(x[inside], mean, sqrt(varcov), log=log)
    pdf[inside] <- if(log) {pdf.Int -log(p.Int)} else {pdf.Int/p.Int}
    return(pdf)
    }
  x <- if (is.vector(x))  t(matrix(x))   else data.matrix(x)
  if(ncol(x) != d)  stop("mismatch of the dimensions of 'x' and 'varcov'")
  if(d > 20) return(rep(NA, NROW(x)))
  if(is.matrix(mean)) {
    if((nrow(x) != nrow(mean)) || (ncol(mean) != d)) 
      stop("mismatch of dimensions of 'x' and 'mean'") }
  ok <- apply((t(x)-lower)>0 & (upper-t(x))>0, 2, all)
  pdf <- rep(0, NROW(x))
  if(sum(ok) > 0) {
    prob <- sadmvn(lower, upper, mean, varcov, ...)
    tmp <- dmnorm(x[ok,], mean, varcov, log=log)
    pdf[ok] <- if(log) {tmp - log(prob)} else {tmp/prob}
    }
  return(pdf)
}

pmtruncnorm <- function(x, mean, varcov, lower, upper, ...) {
  d <- if(is.matrix(varcov)) ncol(varcov) else 1
  if(missing(lower)) lower <- rep(-Inf,d)
  if(missing(upper)) upper <- rep(Inf,d)
  if(length(lower) != d | length(upper) != d) stop("dimension mismatch")
  if(!all(lower < upper)) stop("lower<upper is required (componentwise)")         
  if(d==1) {
    pL <- pnorm(c(lower, upper), mean, sqrt(varcov))
    p <- rep(0, length(c(x)))
    p <- replace(p, x >= upper, 1)
    inside <- (x>lower & x<upper)
    x0 <- x[inside]
    p[inside] <- (pnorm(x0, mean, sqrt(varcov)) - pL[1])/diff(pL)
    return(p)
    }
  x <- if (is.vector(x))  t(matrix(x))     else data.matrix(x)
  if(d > 20) return(rep(NA, NROW(x)))
  if (ncol(x) != d) stop("mismatch of dimensions of 'x' and 'varcov'")
  if (is.matrix(mean)) {
        if ((nrow(x) != nrow(mean)) || (ncol(mean) != d)) 
            stop("mismatch of dimensions of 'x' and 'mean'") }
  n <- NROW(x)  
  p <- numeric(n)
  for(i in 1:n)  p[i] <- if(any(x[i,] < lower)) 0 else 
     sadmvn(lower, pmin(x[i,], upper), mean, varcov)
  return(p/sadmvn(lower, upper, mean, varcov, ...))
}  


rmtruncnorm <- function(n, mean, varcov, lower, upper, start, burnin=5, thinning=1) {
  d <- as.integer(if(is.matrix(varcov)) ncol(varcov) else 1)
  if(missing(lower)) lower <- rep(-Inf, d)
  if(missing(upper)) upper <- rep(Inf, d)
  if(!all(c(length(mean), length(lower), length(upper)) == rep(d,3))) 
    stop("dimension mismatch")
  if(!all(upper>lower)) stop("upper>lower is required")
  if(d == 1) {
      sigma <- c(sqrt(varcov))
      p1 <- pnorm((lower - mean)/sigma)
      p2 <- pnorm((upper - mean)/sigma)
      u <- runif(n)
      return(mean + sigma * qnorm(p1 + u*(p2-p1)))
    }
  if(missing(start)) start <- 
    mom.mtruncnorm(powers=1, mean, varcov, lower, upper, cum=TRUE)$cum1
  regr <- matrix(0, d, d-1)
  sd_c <- numeric(d)
  for(j in 1:d) { 
    r <- c(varcov[j, -j, drop=FALSE] %*% solve(varcov[-j, -j, drop=FALSE]))
    regr[j,] <- r
    sd_c[j] <- sqrt(varcov[j,j] - sum(r * varcov[-j, j]))
  }
  nplus<- as.integer(burnin + n * thinning)
  x <- matrix(0, nplus, d)
  a <- .Fortran("rtmng", nplus, d, mean, regr, sd_c, lower, upper, x, start, 
                NAOK=TRUE, PACKAGE="mnormt")
  x <- a[[8]]
  x[burnin+(1:n)*thinning, , drop=TRUE]
}

#------------------

mom.mtruncnorm <- function(powers=4, mean, varcov, lower, upper, cum=TRUE, ...)
{
  d <- if(is.matrix(varcov)) ncol(varcov) else 1
  if(d > 20) return(NA) # stop("maximal dimension is 20")
  if(any(powers < 0) | any(powers != round(powers)))
    stop("'powers' must be non-negative integers")
  if(length(powers) == 1) powers <- rep(powers, d)  
  if(missing(lower)) lower <- rep(-Inf,d)
  if(missing(upper)) upper <- rep(Inf,d)
  if(!all(lower < upper)) stop("lower<upper is required")   
  if(!all(c(length(powers) == d, length(lower) == d, length(upper) == d, 
     length(mean) == d, dim(varcov) == c(d,d))))  stop("dimension mismatch")
  if(any(lower >= upper)) stop("non-admissible bounds")   
  M <- recintab(kappa=powers, a=lower, b=upper, mean, varcov, ...)
  mom <- M/M[1]
  out <- list(mom=mom)
  # if(mom[1] != 1) { warning("mom[1] != 1"); return(c(out, NA)) }
  cum <- if(cum) mom2cum(mom) else NULL
  return(c(out, cum))
}

mom2cum <- function(mom){
# convert array of multivariate moments to cumulants, up to 4th order maximum
  get.entry <- function(array, subs, val) {
     # get entries with subscripts 'subs' equal to 'val' of 'array' (char)
     x <- get(array)
     ind <- rep(1, length(dim(x)))
     ind[subs] <- val + 1
     subs.char <- paste(as.character(ind), collapse=",")
     eval(str2expression(paste(array, "[", subs.char, "]", sep="")))
     } 
  if(is.na(mom[1])) return(list(cum=NA, message="mom[1] must be 1"))
  if(mom[1] != 1) return(list(cum=NA, message="mom[1] must be 1"))  
  if(is.vector(mom)) { # case d=1 treated separately
    m <- mom[-1]
    powers <- length(m)
    cum <- cmom <- g1 <- g2 <- NULL
    if(powers >= 1) {
       cum <- m[1]
       cmom <- 0
       }
    if(powers >= 2) {
       cum <- c(cum, m[2] - m[1]^2)
       if(cum[2] <= 0 ) warning("cum[2] <= 0")
       cmom <- c(cmom, cum[2])
       }
    if(powers >= 3) {
       cum <- c(cum, m[3] -3*m[1]*m[2] + 2*m[1]^3)
       cmom <- c(cmom, cum[3])
       g1 <- cum[3]/cum[2]^1.5
       }
    if(powers >= 4) {
       cum <- c(cum, m[4] -3*m[2]^2 - 4*m[1]*m[3] + 12*m[1]^2*m[2] -6*m[1]^4)
       cmom <- c(cmom, cum[4] + 3*cum[2]^2)
       g2 <- cum[4]/cum[2]^2
       }
    out <- list(cum=cum, centr.mom=cmom, std.cum=c(gamma1=g1, gamma2=g2))
    return(out)
    } # end of case d=1
  # now case d>1:
  powers <- dim(mom) - 1  
  d <- length(powers)
  out <- list()
  if(all(powers >= 1)) {
    m1 <- numeric(d)
    for(i in 1:d)  m1[i] <- get.entry("mom", i, 1)
    out$cum1  <- m1
    }   
  if(all(powers >= 2)) {
    m2 <- matrix(0, d, d)            # moments of 2nd order
    for(i in 1:d) for(j in 1:d) 
      m2[i,j] <- if(i == j) 
        get.entry("mom", i, 2) else get.entry("mom", c(i, j), c(1,1))
    vcov <- cum2 <- (m2 - m1 %*% t(m1))
    if(any(eigen(cum2, symmetric=TRUE, only.values=TRUE)$values <= 0))
      warning("matrix 'cum2' not positive-definite")
    conc <- pd.solve(vcov, silent=TRUE, log.det=TRUE)
    log.det <- attr(conc, "log.det")
    attr(conc, 'log.det') <- NULL
    out$order2 <- list(m2=m2, cum2=vcov, conc.matrix=conc, log.det.cum2=log.det)
    if(is.null(conc)) { 
      out$message <- "Warning: input array 'mom' appears problematic"  
      return(out)
      }  
    }
  if(all(powers >= 3)) {
    mom2 <- m2[cbind(1:d,1:d)]        # 2nd order marginal moments  
    cmom2 <- vcov[cbind(1:d,1:d)]     # 2nd order marginal central moments  
    # sd <- sqrt(cmom2)
    cmom3 <- mom3 <- numeric(d)       # 3rd order marginal (central) moments 
    m3 <- array(NA, rep(d,3))         # array of 3rd order moments 
    for(i in 1:d) for (j in 1:d) for(k in 1:d) { 
      if(i==j & j==k) { 
         subs <- i
         val <- 3
         mom3[i] <- get.entry("mom", subs, val)
         # next line uses (15.4.4) of Cramér (1946, p.175)
         cmom3[i] <- mom3[i] - 3*m1[i]*mom2[i] + 2*m1[i]^3
         }
      else {
        if (i==j | i==k | j==k) {
          val <- c(2,1)
          if(i==j) subs <- c(i,k) 
          if(i==k) subs <- c(i,j) 
          if(j==k) subs <- c(j,i) 
        } else {
          subs <- c(i,j,k) 
          val <- c(1,1,1)
        } 
      }
      m3[i,j,k] <- get.entry("mom", subs, val)
      }
    #----
    # compute 3rd order cumulants using (2.7) of McCullagh (1987)
    cum3 <- array(NA, rep(d, 3))         # 3rd order cumulants 
    for(i in 1:d) for (j in 1:d) for(k in 1:d) 
      cum3[i,j,k] <- (m3[i,j,k] 
            - (m1[i]*m2[j,k] + m1[j]*m2[i,k] + m1[k]*m2[i,j])
            + 2 * m1[i]*m1[j]*m1[k])
    # compute Mardia`s gamma_{1,d} as \rho^2_{23} in (2.15) of McCullagh (1987)
    g1 <- 0  
    for(i in 1:d) for (j in 1:d) for(k in 1:d) 
      for(l in 1:d) for (m in 1:d) for(n in 1:d)  
        g1 <- g1 + cum3[i,j,k]*cum3[l,m,n]*conc[i,l]*conc[j,m]*conc[k,n]
    out$order3 <- list(m3=m3, cum3=cum3, m3.marginal=mom3,
      centr.mom3.marginal=cmom3, gamma1.marginal=cmom3/cmom2^(3/2),
      gamma1.Mardia=g1, beta1.Mardia=g1)             
    }  
  if(all(powers >= 4)) {
    cmom4 <- mom4 <- numeric(d)       # marginal 4th order (central) moments
    m4 <- array(NA, rep(d,4))         # array of 4th order moments 
    for(i in 1:d) for (j in 1:d) for(k in 1:d) for(l in 1:d) {    
      if(i==j & j==k & k == l) {
        val <- 4
        subs <- i
        mom4[i] <- get.entry("mom", subs, val)
        # next line uses (15.4.4) of Cramér (1946, p.175)
        cmom4[i] <- mom4[i] - 4*m1[i]*mom3[i] + 6*m1[i]^2*mom2[i] - 3*m1[i]^4
       } 
      else { if(i==j & j==k | i==k & k==l | i==j & j==l | j==k & k==l) {
        val <- c(3, 1)
        if(i==j & j==k) subs <- c(i,l)
        if(i==k & k==l) subs <- c(i,j)
        if(i==j & j==l) subs <- c(i,k)
        if(j==k & k==l) subs <- c(j,i)
        }
      else { if(i==j & k==l | i==k & j==l | i==l & j==k) {
        val <- c(2, 2)
        if(i==j & k==l) subs <- c(i,k)
        if(i==k & j==l) subs <- c(i,j)
        if(i==l & j==k) subs <- c(i,j)
        }  
      else { if(i==j | i==k | i==l | j==k | j==l | k==l) {
        val <- c(2, 1, 1)
        if(i==j) subs <- c(i,k,l)
        if(i==k) subs <- c(i,j,l)
        if(i==l) subs <- c(i,j,k)
        if(j==k) subs <- c(j,i,l)
        if(j==l) subs <- c(j,i,k) 
        if(k==l) subs <- c(k,i,j)
        } 
      else {
        val <- c(1,1,1,1)
        subs <- c(i,j,k,l)
        }}}}
      m4[i,j,k,l] <- get.entry("mom", subs, val)    
      }
    # compute 4th order cumulants using (2.7) of McCullagh (1987)
    cum4 <- array(NA, rep(d, 4))      
    for(i in 1:d) for (j in 1:d) for(k in 1:d) for(l in 1:d)  
      cum4[i,j,k,l] <- ( m4[i,j,k,l] 
        -(m1[i]*m3[j,k,l] + m1[j]*m3[i,k,l] + m1[k]*m3[i,j,l] + m1[l]*m3[i,j,k])
        -(m2[i,j]*m2[k,l] + m2[i,k]*m2[j,l] + m2[i,l]*m2[j,k])
        +2 * (m1[i]*m1[j]*m2[k,l] + m1[i]*m1[k]*m2[j,l] + m1[i]*m1[l]*m2[j,k] +
              m1[j]*m1[k]*m2[i,l] + m1[j]*m1[l]*m2[i,k] + m1[k]*m1[l]*m2[i,j])
        -6 * m1[i]*m1[j]*m1[k]*m1[l]  ) # end cum4[i,j,k,l] 
    # compute Mardia`s gamma_{2,d} as \rho_4 in (2.16) of McCullagh (1987),
    g2 <- 0
    for(i in 1:d) for (j in 1:d)  g2 <- g2  + cum4[i,j,,] * conc * conc[i,j]
    g2 <- sum(g2)
    # for(i in 1:d) for (j in 1:d) for(k in 1:d) for(l in 1:d) 
    #  g2 <- g2 + cum4[i,j,k,l]*conc[i,j]*conc[k,l]     
    b2 <- g2 + d*(d+2)
    out$order4 <- list(m4=m4, cum4=cum4, m4.marginal=mom4,
        centr.mom4.marginal=cmom4, gamma2.marginal=(cmom4/cmom2^2 - 3),
        gamma2.Mardia=g2, beta2.Mardia=b2) 
    }    
  return(out)
}

#------------------
recintab <- function(kappa, a, b, mu, S, ...)
{# R translation of Matlab code in 'recintab.m'
if(!all(a < b)) stop("a<b is required")
d <- n <-  length(kappa);
if (n == 1) {
   M <- rep(0, kappa+1)
   s1 <- sqrt(S);
   aa <- (a - mu)/s1;
   bb <- (b - mu)/s1;
   M[1] <- pnorm(bb) - pnorm(aa);
   if (kappa > 0) {
      pdfa <- s1*dnorm(aa);
      pdfb <- s1*dnorm(bb);
      M[2] <- mu*M[1] + pdfa - pdfb;
      if(is.infinite(a)) a <- 0;
      if(is.infinite(b)) b <- 0;
      if(kappa > 1) for(i in 2:kappa) {
          pdfa <- pdfa*a;
          pdfb <- pdfb*b;
          M[i+1] <- mu*M[i] + (i-1)*S*M[i-1] + pdfa - pdfb;
      }}}
else {
#
#  Create a matrix M, with its nu-th element correpsonds to F_{nu-2}^n(mu,S).
#
   M <- array(0, dim=kappa+1);
   pk <- prod(kappa+1);
#
#  We create two long vectors G and H to store the two different sets
#  of integrals with dimension n-1.
#
   nn <- round(pk/(kappa+1));   
   begind <- cumsum(c(0, nn));  
   pk1 <- begind[n+1];   # number of (n-1)-dimensional integrals
#  Each row of cp corresponds to a vector that allows us to map the subscripts
#  to the index in the long vectors G and H
   cp <- matrix(0, n, n);           
   for(i in 1:n) {
       kk <- kappa;
       kk[i] <- 0;
       cp[i,] <-  c(1, cumprod(kk[1:(n-1)] + 1));   
       }
   G <- rep(0, pk1);
   H <- rep(0, pk1);
   s <- sqrt(diag(S));
   pdfa <- dnorm(a, mu, s)
   pdfb <- dnorm(b, mu, s)
   for(i in 1:n) {   
       ind2 <- (1:n)[-i];
       # ind2(i) <- [];
       kappai <- kappa[ind2];
       ai <- a[ind2];
       bi <- b[ind2];
       mui <- mu[ind2];
       Si <- S[ind2,i];
       SSi <- S[ind2,ind2] - Si %*% t(Si)/S[i,i];      
       ind <- (begind[i]+1):begind[i+1];
       if(a[i] > -Inf) {
          mai <- mui + Si/S[i,i] * (a[i]-mu[i]);
          G[ind] <- pdfa[i] * recintab(kappai, ai, bi, mai, SSi);
          }
       if(b[i] < Inf ) {
          mbi <- mui + Si/S[i,i] * (b[i]-mu[i]);
          H[ind] <- pdfb[i] * recintab(kappai, ai, bi, mbi, SSi);
          }
   } 
#
#  Use recursion to obtain M(nu).
   M[1] <- if(n < 3) biv.nt.prob(Inf, a, b, mu, S) else sadmvn(a, b, mu, S, ...) 
   a[is.infinite(a)] <- 0;
   b[is.infinite(b)] <- 0;    
   cp1 <- t(cp[n,,drop=FALSE]);
   for(i in 2:pk) {
       kk <- arrayInd(i, kappa+1)   
       ii <- (kk-1) %*% cp1 + 1;
       i1 <- min(which(kk>1));  # find a nonzero element to start the recursion
       kk1 <- kk;
       kk1[i1] <- kk1[i1] - 1;
       ind3 <- ii - cp1[i1];
       M[ii] <- mu[i1] %*% M[ind3];
       for(j in 1:n) {
           kk2 <- kk1[j] - 1;
           if(kk2 > 0)
              M[ii] <- M[ii] + S[i1,j] %*% kk2 %*% M[ind3-cp1[j]];
           ind4 <- begind[j] + cp[j,] %*% t(kk1-1) - cp[j,j]*kk2 + 1;
           M[ii] <- M[ii] + S[i1,j] %*% (a[j]^kk2*G[ind4] -b[j]^kk2*H[ind4]);
        }
      }
  }      
  return(M) 
}
#--------------------------------------------
# multivariate truncated t distribution
#
dmtrunct <- function(x, mean, S, df, lower, upper, log= FALSE, ...) {
  if(df == Inf)  return(dmtruncnorm(x, mean, S, log = log))
  d <- if(is.matrix(S)) ncol(S) else 1
  x <- if (is.vector(x))  t(matrix(x))  else data.matrix(x)
  if(ncol(x) != d) stop("mismatch of dimensions of 'x' and 'S'")
  if(d > 20) return(rep(NA, NROW(x))) # stop("maximal dimension is 20")
  if(is.matrix(mean)) {
    if((nrow(x) != nrow(mean)) || (ncol(mean) != d)) 
      stop("mismatch of dimensions of 'x' and 'mean'")}
  if(missing(lower)) lower <- rep(-Inf,d)
  if(missing(upper)) upper <- rep(Inf,d)
  if(length(lower) != d | length(upper) != d) stop("dimension mismatch")            
  if(!all(lower < upper)) stop("lower<upper is required")
  ok <- apply((t(x)-lower)>0 & (upper-t(x))>0, 2, all)
  pdf <- rep(0, NROW(x))
  if(sum(ok) > 0) {
    prob <- sadmvt(df, lower, upper, mean, S, ...)
    tmp <- dmt(x[ok,], mean, S, df, log=log)
    pdf[ok] <- if(log) tmp - log(prob) else tmp/prob
    }
  return(pdf)
}

pmtrunct <- function(x, mean, S, df, lower, upper, ...) {
  if(df == Inf)  return(pmtruncnorm(x, mean, S, log = log))
  d <- if(is.matrix(S)) ncol(S) else 1
  x <- if (is.vector(x))  t(matrix(x))  else data.matrix(x)
  if(d > 20) return(rep(NA, NROW(x))) # stop("maximal dimension is 20")
  if (ncol(x) != d) stop("mismatch of dimensions of 'x' and 'S'")
  if (is.matrix(mean)) {
        if ((nrow(x) != nrow(mean)) || (ncol(mean) != d)) 
            stop("mismatch of dimensions of 'x' and 'mean'") }
  if(missing(lower)) lower <- rep(-Inf,d)
  if(missing(upper)) upper <- rep(Inf,d)
  if(length(lower) != d | length(upper) != d) stop("dimension mismatch")            
  if(!all(lower < upper)) stop("lower<upper is required")
  n <- NROW(x)  
  p <- numeric(n)
  for(i in 1:n)  p[i] <- if(any(x[i,] < lower)) 0 else 
     sadmvt(df, lower, pmin(x[i,], upper), mean, S, ...)
  return(p/sadmvt(df,lower, upper, mean, S, ...))
}  

