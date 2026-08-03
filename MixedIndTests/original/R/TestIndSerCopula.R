#'@title Statistics and P-values for a test of randomness for a univariate time series
#'
#'@description This function computes Cramer-von Mises statistics from the multilinear copula and their combination for tests of randomness of p consecutives values X(1), ..., X(p). The p-values are computed using Gaussian multipliers.
#'
#'@param x            Time series
#'@param p            Number of consecutive observations
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#'@param B            Number of multipliers samples (default = 1000)
#'@param par          Set to TRUE if one prefers paraller computing (slower)
#'@param ncores       Number of cores for parallel computing (default = 2)
#'@param graph        Set to TRUE if one wants the dependogram of P-values for the Moebius statistics
#'
#'
#'@return \item{stat}{List of Cramer-von Mises statistics cvm, Sn, and test combinations Tn and Tn2 (only pairs)}
#'@return \item{pvalue}{Approximated P-values for the tests using Gaussian multipliers}
#'@return \item{card}{Cardinaly of the subsets for the Moebius statistics}
#'@return \item{subsets}{Subsets for the Moebius statistics}
#'
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'
#'@importFrom foreach %dopar%
#'@import doParallel

#'@examples
#' X <- SimAR1Poisson(c(5,0.2),100)
#' out <- TestIndSerCopula(X,5,3)



TestIndSerCopula =function(x,p,trunc.level=2,B=1000,par=FALSE,ncores=2,graph=FALSE){


 # start_time <- Sys.time()

  v = c(1:(trunc.level-1))
  cA = sum(choose((p-1),v))
  m = 2^(p-1)-1
  n = length(x)
  y = sort(unique(x))
  m0 = length(y)


  if(m0==2){
    out0 <- .C("stats_serial_bin",
               as.double(x),
               as.integer(n),
               as.integer(p),
               p0 = double(1),
               ZA    = double(m),
               J = double(m*m),
               Asets = double(p*m),
               cardA0  = double(m),
               PACKAGE = "MixedIndTests"
    )
    #library(survey)
    Asets0=matrix(out0$Asets,ncol=p)
    Sn = Sn_serial(x,p)$Sn
    J = matrix(out0$J,ncol=m)
    ll = eigen(J)$values
    df = c(rep(1,m))
    chi2 = rep(0,cA)
    cvm  = rep(0,cA)
    cardA = rep(0,cA)
    p0 = out0$p0
    cte = p0*(1-p0)/3.0
    ZA = out0$ZA
    cardA0=out0$cardA0
    Asets=matrix(0,nrow=cA,ncol=p)
    l=1
    for(k in 1:m)
    {
      if(cardA0[k]<=trunc.level)
      {
        cardA[l]=cardA0[k]
        chi2[l] = ZA[k]*ZA[k]
        cvm[l] = chi2[l]*cte^(cardA0[k])
        Asets[l,] = Asets0[k,]
        l = l+1;
      }
    }
    cvm0sim=0*c(1:B)
    for(it in 1:B)
    {
      xi = rnorm(m)
      cvm0sim[it] = t(xi)%*%J%*%xi
    }
    pvalSn = 100*mean(cvm0sim>Sn)
    #pvalSn = 100*survey::pchisqsum(Sn, df, ll, lower.tail = F, method = "saddlepoint")
    pvalcvm=c(rep(0,cA))
    for(k in 1:cA){
      pvalcvm[k]=100*pchisq(chi2[k], 1, lower.tail = F)
    }
  }
  else{

  out0 <- .C("stats_serial",
            as.double(x),
            as.integer(n),
            as.integer(p),
            as.integer(trunc.level),
            stats = double(cA),
            cardA = double(cA),
            M = double(n*n*cA),
            Asets = double(p*cA),
            Sn    = double(1),
            J = double(n*n),
            PACKAGE = "MixedIndTests"
  )

  #end_time <- Sys.time()
  #print(end_time - start_time)
  Asets = matrix(out0$Asets,ncol=p)
  cardA=out0$cardA
  cvm = out0$stats
  Sn = out0$Sn
  Mvec = out0$M
  Jvec = out0$J
  cA = length(out0$card)
  n = length(x)
  # Bootstrapping
  pvalcvm = 0*c(1:cA)

  if(par){
    ## Parallel bootstrapping
    #ncores <- max(2,parallel::detectCores()-2);
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    fun <- c('bootstrap')

 #   start_time <- Sys.time()

    result <- foreach::foreach(i=1:B, .export=fun, .packages = "MixedIndTests") %dopar% bootstrap(Mvec,Jvec,n)

    parallel::stopCluster(cl)
    cvm0sim = 0*c(1:B)
    cvm_sim <- matrix(0,nrow=B,ncol=cA)
    for (i in 1:B){
      cvm_sim[i,] <- result[[i]]$cvm
      cvm0sim[i]  <- result[[i]]$Sn

    }
    for(k in 1:cA){
      pvalcvm[k]=100*mean(cvm_sim[,k]>cvm[k])
    }
    pvalSn = 100*mean(cvm0sim>Sn)


    # end_time <- Sys.time()
    # print(end_time - start_time)
  }else{
  z = rnorm(n*B)
  J  = matrix(Jvec,nrow=n)
  sim=0*c(1:B)
  cvm0sim=0*c(1:B)
  n2 = n*n
  M=array(Mvec,c(cA,n,n))
  for (k in 1:cA)
  {
    M0 = M[k,,]
    for(it in 1:B)
      {
         xi = z[(n*(it-1)+1):(n*it)]
         xic = xi-mean(xi)
         sim[it] = t(xic)%*%M0%*%xic/n
      }
    pvalcvm[k]=100*mean(sim>cvm[k])
  }
  rm(M0)


  for(it in 1:B)
    {
    xi = z[(n*(it-1)+1):(n*it)]
    xic = xi-mean(xi)
    cvm0sim[it] = t(xic)%*%J%*%xic/n
    }
  pvalSn = 100*mean(cvm0sim>Sn)
  #rm(J)
}
  # Much slower!!!
  # start_time <- Sys.time()
  #  for(it in 1:B){
  #
  #    out1=bootstrap(M,n,p,trunc.level)
  #
  #  }
  # end_time <- Sys.time()
  # print(end_time - start_time)

  # rm(Mvec)
  # rm(M)
  # rm(M0)
  # rm(Jvec)

  # Must compute the simulated statistics before

}
#######################################################
  ind = sort(cardA,index.return=TRUE)
  card = ind$x
  Asets = out0$Asets
  AA = matrix(Asets,ncol=p)
  if(cA==1)
  {
    A=AA
  }else{
    A=AA[ind$ix,]
  }

  stats.cvm = cvm[ind$ix]
  pval.cvm = pvalcvm[ind$ix]

  pval = pval.cvm/100+1e-20
  Tn = -2*sum(log(pval))
  Tn2 = -2*sum(log(pval[card==2]))
  pvalTn = 100*(1-pchisq(Tn,2*cA));
  pvalTn2 = 100*(1-pchisq(Tn2,2*(p-1)));
  pvalue=list(cvm=pval.cvm,Sn=pvalSn,Tn=pvalTn,Tn2=pvalTn2)
  stat=list(cvm=stats.cvm,Sn=Sn,Tn=Tn,Tn2=Tn2)
  lagsets = vector(mode="character",length = cA)
  for(i in 1:cA){
    lagsets[i]=as.character(1)
  }
  for ( i in 1:cA)
  {
    for(j in 2:p)
    {
      if(A[i,j])
      {
        lagsets[i] = paste(lagsets[i],as.character(j))
      }
    }

  }


  out0 = list(stat=stat, pvalue=pvalue, card=card,subsets=lagsets)

  #end_time <- Sys.time()
  #print(end_time - start_time)
  if(graph){
    Dependogram(out0)

  }

   return(out0)
}






