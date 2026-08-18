draw.correlated.binary<- function (no.row, d, prop.vec, corr.mat) 
{
  if ((no.row < 1) | (floor(no.row) != no.row)) {
    stop("Number of replicates must be an integer whose value is at least 1!\n")
  }
  if ((d < 2) | (floor(d) != d)) {
    stop("Dimension must be an integer whose value is at least 2!\n")
  }
  if ((max(prop.vec) >= 1) | (min(prop.vec) <= 0)) {
    stop("Expectations should be greater than 0 and less than 1!\n")
  }
  if (length(prop.vec) != d) {
    stop("Prop vector is misspecified, dimension is wrong!\n")
  }
  if ((ncol(corr.mat) != d) | (nrow(corr.mat) != d)) {
    stop("Correlation matrix is misspecified, dimension is wrong!\n")
  }
  if (sum(corr.mat != t(corr.mat)) > 0) {
    stop("Correlation matrix is not symmetric!\n")
  }
  if (sum(diag(corr.mat) != rep(1, d)) > 0) {
    stop("Not all diagonal elements of correlation matrix are 1!\n")
  }
  if ((max(corr.mat) > 1) | (min(corr.mat) < 0)) {
    stop("Correlations should be greater than or equal to 0 and less than or equal to 1!\n")
  }
  alpha = matrix(0, d, d)
  cor.limit = matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      cor.limit[i, j] = min(sqrt((prop.vec[j] * (1 - prop.vec[i]))/(prop.vec[i] * 
                                                                      (1 - prop.vec[j]))), sqrt((prop.vec[i] * (1 - 
                                                                                                                  prop.vec[j]))/(prop.vec[j] * (1 - prop.vec[i]))))
    }
  }
  if (sum(cor.limit >= corr.mat) < d^2) {
    stop("Correlations are beyond their upper limits imposed by expectations")
  }
  for (i in 1:d) {
    for (j in 1:d) {
      alpha[i, j] = log(1 + corr.mat[i, j] * sqrt((1 - 
                                                     prop.vec[i]) * (1 - prop.vec[j])/(prop.vec[i] * 
                                                                                         prop.vec[j])))
    }
  }
  
  if (d==2) {
    
    x = matrix(0, no.row, d)
    y = matrix(0, no.row, d)
    pois = numeric(d)
    sump = numeric(d)
    for (k in 1:no.row) {
      
      pois[1] = rpois(1, alpha[1,1] - alpha[1,2])  
      pois[2] = rpois(1, alpha[2,2] - alpha[1,2]) 
      pois[3] = rpois(1, alpha[1,2])
      
      sump[1] = pois[1] + pois[3]
      sump[2] = pois[2] + pois[3]
      
      
      x[k, ] = sump
    }
    y[x == 0] = 1
    y[x != 0] = 0
    y
    
  } else {
    
    beta = matrix(0, d, d * d)
    summ = 1
    counter = 0
    while ( summ > 0 ) {
      counter = counter + 1
      minloc = loc.min(alpha, d)
      w = matrix(1, d, d)
      mat.min = apply(alpha, 2, min)
      pos = c(1:d)
      zero.pos = which(mat.min==0)
      if ( length(zero.pos) == 0 ) {
        nonzero.pos = pos
      } else {
        nonzero.pos = pos[-zero.pos]
      }
      
      my.min = apply(matrix(alpha[, -minloc], d, d - length(unique(minloc))), 
                     2, min)
      if (length(my.min) == 1) {
        w[, -minloc][my.min == 0] = 0
        w[-minloc, ][my.min == 0] = 0
      }
      if (length(my.min) > 1) {
        w[, -minloc][, my.min == 0] = 0
        w[-minloc, ][my.min == 0, ] = 0
        w[alpha == 0] = 0
      }
      for (i in 1:d) {
        # if ( sum( mat.min != 0) =  d ) {
        #   T = rep(1,d)
        # }
        beta[i, counter] = alpha[minloc[1], minloc[2]] * 
          1 * ((minloc[1] == i) | (minloc[2] == i) | (i %in% nonzero.pos))
        
      }
      
      if ( 0 %in% diag(alpha)[minloc] ) {
        stop("The method won't work in this parameter setting\n")
      }
      
      alpha = alpha - alpha[minloc[1], minloc[2]] * w
      summ = sum(alpha)
      
    }
    
    tbeta = t(beta)
    w = (tbeta != 0)
    x = matrix(0, no.row, d)
    y = matrix(0, no.row, d)
    pois = numeric(nrow(tbeta))
    sump = numeric(d)
    for (k in 1:no.row) {
      for (j in 1:nrow(tbeta)) {
        pois[j] = rpois(1, max(tbeta[j, ]))
      }
      for (i in 1:d) {
        sump[i] = sum(pois * w[, i])
      }
      x[k, ] = sump
    }
    y[x == 0] = 1
    y[x != 0] = 0
    y
  }
}
