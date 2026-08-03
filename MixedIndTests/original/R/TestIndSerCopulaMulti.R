#'@title Statistics and P-values for a test of randomness for a multivariate time series
#'
#'@description This function computes Cramer-von Mises statistics from the multilinear copula and their combination for a tests of randomness for p consecutives values of random vectors X(1), ..., X(p). The p-values are computed using Gaussian multipliers.
#'
#'@param x            Time series matrix
#'@param p            Number of consecutive vectors
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#'@param B            Number of multipliers samples (default = 1000)
#'@param graph        Set to TRUE if one wants the dependogram of P-values for the Moebius statistics
#'
#'
#'@return \item{stat}{List of Cramer-von Mises statistics cvm, tilde Sn, and test combinations tilde Tn and tilde Tn2 (only pairs), as defined in Nasri(2022).}
#'@return \item{pvalue}{Approximated P-values for the tests using Gaussian multipliers}
#'
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'
#'@examples
#'\donttest{
#'data(Y)
#'out <- TestIndSerCopulaMulti(Y,5,5)
#'}



TestIndSerCopulaMulti =function(x,p,trunc.level=2,B=1000,graph=FALSE){


 # start_time <- Sys.time()

  v = c(1:(trunc.level-1))
  cA = sum(choose((p-1),v))
  m = 2^(p-1)-1
  dd = dim(x)
  n = dd[1]
  d = dd[2]

  Stat1 = 0.0;

  Sn1   = 0.0;

  J1    = 0.0;

  M1 = array(0,c(cA,n,n))
   for(k in 1:d)
   {
     out0=Sn_Aserial(x[,k],p,trunc.level)
     Mvec = out0$M
     MM1=array(Mvec,c(cA,n,n))
     M1 = M1+MM1
     Jvec = out0$J
     J = matrix(Jvec,nrow=n)
     J1=J1+J;

     Stat1 = Stat1+out0$stats
     Sn1   = Sn1+ out0$Sn
   }
  for(k in 1:(d-1))
  {
     for(j in (k+1):d)
     {
       out0=Sn_AserialVec(x[,c(k,j)], p, trunc.level)

       Mvec = out0$M
       MM1=array(Mvec,c(cA,n,n))
       M1 = M1+MM1
       Jvec = out0$J
       J = matrix(Jvec,nrow=n)
       J1=J1+J;

       Stat1 = Stat1+out0$stats
       Sn1   = Sn1+ out0$Sn

     }
   }





  cardA=out0$cardA
  cvm = Stat1
  Sn  = Sn1


  # Bootstrapping
  pvalcvm = 0*c(1:cA)


  z = rnorm(n*B)




  sim=0*c(1:B)
  cvm0sim=0*c(1:B)


  for (k in 1:cA)
  {
    M0 = M1[k,,]
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
    cvm0sim[it] = t(xic)%*%J1%*%xic/n
    }
  pvalSn = 100*mean(cvm0sim>Sn)



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


  if(graph){
    Dependogram(out0)

  }

   return(out0)
}






