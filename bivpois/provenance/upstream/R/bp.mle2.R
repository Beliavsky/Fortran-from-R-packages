########################
#### MLE of a bivariate poisson distribution
#### 3/2015
#### mtsagris@yahoo.gr
#### References: Kazutomo Kawamura (1984)
#### Direct calculation of maximum likelihood
#### estimator for the bivariate poisson distribution
#### Kodai mathematical journal
################################

bp.mle2 <- function(x1, x2 = NULL) {
  if ( is.null(x2) ) {
    x2 <- x1[, 2]
    x1 <- x1[, 1]
  }
  x1 <- as.numeric(x1)   ;    x2 <- as.numeric(x2)

  ## x1 and x2 are the two variables
  n <- length(x1)  ## sample size
  sx1 <- sum(x1)   ;  sx2 <- sum(x2)
  m1 <- sx1 / n    ;    m2 <- sx2 / n
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

  bar <- optimize( funa, c(0, min(m1, m2) - 0.05), f2a = f2a, n = n, tol = 1e-6, maximum = TRUE)

  l3 <- bar$maximum  ## maximum of the log-likelihood
  lambda <- c(m1 - l3, m2 - l3, l3)
  names(lambda) <- c('lambda1', 'lambda2', 'lambda3')
  list(lambda = lambda, loglik = bar$objective)
}
