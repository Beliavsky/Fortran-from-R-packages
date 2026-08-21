########################
#### MLE of a bivariate poisson distribution
#### 3/2015
#### mtsagris@yahoo.gr
#### References: Kazutomo Kawamura (1984)
#### Direct calculation of maximum likelihood
#### estimator for the bivariate poisson distribution
#### Kodai mathematical journal
################################

bp.mle <- function(x1, x2 = NULL) {
  if ( is.null(x2) ) {
    x2 <- x1[, 2]
    x1 <- x1[, 1]
  }
  x1 <- as.numeric(x1)   ;    x2 <- as.numeric(x2)

  ## x1 and x2 are the two variables
  n <- length(x1)  ## sample size
  sx1 <- sum(x1)  ;  sx2 <- sum(x2)
  m1 <- sx1 / n   ;  m2 <- sx2 / n
  ##  m1 and m2 estimates of lambda1* and lambda2* respectively
  ##  funa is the function to be maximised over lambda3
  ind <- Rfast::rowMins( cbind(x1, x2), value = TRUE )
  max1 <- max(x1)  ;    max2 <- max(x2)
  mm <- max( max1, max2 )   ;   mn <- min(max1, max2)
  omn <- 0:mn
  fac <- factorial( omn )
  #ch <- matrix(numeric( (mm + 1)^2 ), nrow = mm + 1, ncol = mm + 1 )
  #for ( i in 1:c(mm + 1) ) {
  #  for ( j in c(i - 1):c(mm + 1) ) {
  #    ch[i, j] <- choose(j, i - 1)
  #  }
  #}
  i <- j <- 1:c(mm + 1)
  ch <- choose( Rfast::rep_row(j, mm + 1), i - 1 )
  rownames(ch) <- colnames(ch) <- 0:mm
  sly1 <- sum( lgamma(x1 + 1) )
  sly2 <- sum( lgamma(x2 + 1) )

  f2a <- list()
  for (j in 1:n) {
    a <-  1:c(ind[j] + 1)
    f2a[[ j ]] <- ch[ a, x1[j] ] * ch[ a, x2[j] ] * fac[ a ]
  }

  funa <- function(l3, f2a, n) {
    f2 <- numeric(n)
    con <-  - m1 - m2 + l3
    expo <- ( l3/( (m1 - l3) * (m2 - l3) ) )^omn
    l1 <- log(m1 - l3)
    l2 <- log(m2 - l3)
    f1 <- sx1 * l1 - sly1 + sx2 * l2 - sly2
    for (j in 1:n) {
      f2[j] <- log( sum( f2a[[ j ]] * expo[ 1:c(ind[j] + 1) ] ) )
    }
    n * con + f1 + sum( f2[abs(f2) < Inf] )
  }

  bar <- optim( cov(x1, x2), funa, f2a = f2a, n = n, control = list(fnscale = -1),
                method = "L-BFGS-B", lower = 0, upper = min(m1, m2) - 0.05, hessian = TRUE )
  l1 <- bar$value  ## maximum of the log-likelihood
  l3 <- bar$par  ## lambda3 estimate
  rho <- l3 / sqrt(m1 * m2)  ## correlation coefficient
  names(rho) <- "correlation"
  l0 <- funa(0, f2a, n)  ## log-likelihood with lam3=0, independence
  stat <- 2 * (l1 - l0)  ## log-likelihood ratio test
  pval1 <- pchisq(stat, 1, lower.tail = FALSE)
  ma <- mm + 20
  f1 <- f2 <- matrix(nrow = ma, ncol = ma)
  con <-  - m1 - m2 + l3

  for ( r in 1:ma ) {
    for ( s in 1:ma ) {
      i <- 0:min(r, s)
      comon <- factorial(i) * ( l3/( (m1 - l3) * (m2 - l3) ) )^i
      f1[r, s] <-  con + (r - 1) * log(m1 - l3) - lgamma(r) +
                   (s - 1) * log(m2 - l3) - lgamma(s) + 
                   log( sum( choose(r - 1, i) * choose(s - 1, i) * comon ) )
      f2[r, s] <- con + r * log(m1 - l3) - lgamma(r + 1) +
                  s * log(m2 - l3) - lgamma(s + 1) + 
                  log( sum( choose(r, i) * choose(s, i) * comon ) )
    }
  }

  tau <- sum( exp(f1)^2 / exp(f2) )
  d2 <-  -m1 + l3 - m2 + l3 + (m1 * m2 - l3^2) * (tau - 1)
  s <- matrix( c(m1, l3, l3, l3, m2, l3, l3, l3,
                ( (m1 - l3) * (m2 - l3) + l3 * (m1 - l3 + m2 - l3) *
                  (l3 * (tau - 1) - 1) ) / d2 ), ncol = 3 ) / n
  v1 <-  -1/bar$hessian ## variance of l3 using the observed information matrix
  v2 <- s[3, 3] ## variance of l3 using the asymptotic covariance matrix
  t1 <- l3 / sqrt(v1)  ## Wald test 1
  t2 <- l3 / sqrt(v2)  ## wald test 2
  pval2 <- pnorm(-t1)
  pval3 <- pnorm(-t2)
  pvalue <- c(pval1, pval2, pval3)

  names(pvalue) <- c('LLR', 'Wald 1', 'Wald 2')
  ci <- rbind( c( l3 - 1.959964 * sqrt(v1), l3 + 1.959964 * sqrt(v1) ),
               c( l3 - 1.959964 * sqrt(v2), l3 + 1.959964 * sqrt(v2) ) )
  colnames(ci) <- c('2.5%', '97.5%')
  rownames(ci) <- c('Observed I', 'Asymptotic I')

  loglik <- c(l0, l1)
  names(loglik) <- c('loglik0', 'loglik1')
  lambda <- c(m1 - l3, m2 - l3, l3)
  names(lambda) <- c('lambda1', 'lambda2', 'lambda3')
  list(lambda = lambda, rho = rho, ci = ci, loglik = loglik, pvalue = pvalue)
}
