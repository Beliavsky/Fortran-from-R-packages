dbp <- function(x1, x2 = NULL, lambda, logged = TRUE) {
  if ( is.null(x2) ) {
    x2 <- x1[, 2]
    x1 <- x1[, 1]
  }
  x1 <- as.numeric(x1)   ;    x2 <- as.numeric(x2)

  ## x1 and x2 are the two variables
  n <- length(x1)  ## sample size

  l1 <- lambda[1]  ;  l2 <- lambda[2]  ;  l3 <- lambda[3]
  m1 <- l1 + l3    ;  m2 <- l2 + l3

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
  ly1 <- lgamma(x1 + 1)
  ly2 <- lgamma(x2 + 1)

  f1 <- f2 <- numeric(n)
  con <-  - m1 - m2 + l3
  expo <- ( l3/( (m1 - l3) * (m2 - l3) ) )^omn
  l1 <- log(m1 - l3)
  l2 <- log(m2 - l3)
  for (j in 1:n) {
    f1[j] <- x1[j] * l1 - ly1[j] + x2[j] * l2 - ly2[j]
    f2[j] <- log( sum( ch[ 1:c(ind[j] + 1), x1[j] ] * ch[ 1:c(ind[j] + 1), x2[j] ] *
                       fac[ 1:c(ind[j] + 1) ] * expo[ 1:c(ind[j] + 1) ] ) )
  }

  ep <- which( is.infinite(f2) )
  if ( length(ep) > 0 )  f2[ep] <- 0
  den <- con + f1 + f2
  if ( !logged )  den <- exp(den)
  den
}
