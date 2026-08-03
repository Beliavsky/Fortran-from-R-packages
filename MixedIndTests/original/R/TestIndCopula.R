#'@title Statistics and P-values for a test of independence between random variables
#'
#'@description This function computes Cramer-von Mises statistics and their combination for a tests of independence between random variables with arbitrary distributions. The P-values are computed using Gaussian multipliers.
#'
#'@param x            Data matrix
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#'@param B            Number of multipliers samples (default = 1000)
#'@param par          Set to TRUE if one prefers paraller computing (slower)
#'@param ncores       Number of cores for parallel computing (default is 2)
#'@param graph        Set to TRUE if one wants the dependogram of P-values for the Moebius statistics
#'
#'@return \item{stat}{List of Cramer-von Mises statistics cvm, Sn  from the multilinear copula, and test combinations Tn and Tn2 (only pairs)}
#'@return \item{pvalue}{Approximated P-values for the tests using Gaussian multipliers}
#'@return \item{card}{Cardinaly of the subsets for the Moebius statistics}
#'@return \item{subsets}{Subsets for the Moebius statistics}
#'
#'@references Genest, Neslehova, Remillard & Murphy (2019). Testing for independence in arbitrary distributions
#'
#'@importFrom foreach %dopar%
#'@import doParallel
#'@examples
#' x <- matrix(rnorm(250),ncol=5)
#' out <-TestIndCopula(x)




TestIndCopula =function(x,trunc.level=2,B=1000,par=FALSE,ncores=2,graph=FALSE){

 if(is.data.frame(x)){x=as.matrix(x)}
 # start_time <- Sys.time()

  dim0 <- dim(x)
  n <- dim0[1]
  d <- dim0[2]
  v = c(2:(trunc.level))

  cA = sum(choose(d,v))

  out0 = .C("stats_nonserial",
            as.double(x),
            as.integer(n),
            as.integer(d),
            as.integer(trunc.level),
            stats = double(cA),
            cardA = double(cA),
            M = double(n*n*cA),
            Asets = double(d*cA),
            Sn    = double(1),
            J = double(n*n),
            PACKAGE = "MixedIndTests"
  )

  #end_time <- Sys.time()
  #print(end_time - start_time)


  cvm = out0$stats
  Sn = out0$Sn
  Mvec = out0$M
  Jvec = out0$J
  m = length(out0$card)

  # Bootstrapping
  pvalcvm = 0*c(1:m)

  if(par){
    #    ncores <- max(2,parallel::detectCores()-2);
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)


    fun <- c('bootstrap')

 #   start_time <- Sys.time()

    result <- foreach::foreach(i=1:B, .export=fun, .packages = "MixedIndTests") %dopar% bootstrap(Mvec,Jvec,n)

    parallel::stopCluster(cl)

    cvm0sim = 0*c(1:B)
    cvm_sim <- matrix(0,nrow=B,ncol=m)
    for (i in 1:B){
      cvm_sim[i,] <- result[[i]]$cvm
      cvm0sim[i]  <- result[[i]]$Sn

    }
    for(k in 1:m){
      pvalcvm[k]=100*mean(cvm_sim[,k]>=cvm[k])
    }
    pvalSn = 100*mean(cvm0sim>=Sn)

    # end_time <- Sys.time()
    # print(end_time - start_time)
  }else{
  z = rnorm(n*B)
  J  = matrix(Jvec,nrow=n)
  sim=0*c(1:B)
  cvm0sim=0*c(1:B)
  n2 = n*n
  M=array(Mvec,c(m,n,n))
  for (k in 1:m)
  {
    M0 = M[k,,]
    for(it in 1:B)
      {
         xi = z[(n*(it-1)+1):(n*it)]
         xic = xi-mean(xi)
         sim[it] = t(xic)%*%M0%*%xic/n
      }
    pvalcvm[k]=100*mean(sim>=cvm[k])
  }
  #rm(M0)


  for(it in 1:B)
    {
    xi = z[(n*(it-1)+1):(n*it)]
    xic = xi-mean(xi)
    cvm0sim[it] = t(xic)%*%J%*%xic/n
    }
  pvalSn = 100*mean(cvm0sim>=Sn)
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

  # rm(M)
  #
  # rm(Jvec)

  # Must compute the simulated statistics before



  ind = sort(out0$cardA,index.return=TRUE)
  card = ind$x
  Asets = out0$Asets
  AA = matrix(Asets,ncol=d)
  if(m==1)
  {
    A=AA
  }else{
    A=AA[ind$ix,]
    }

  stats.cvm = cvm[ind$ix]
  pval.cvm = pvalcvm[ind$ix]
  lim = 1e-20*rep(1,m)
  pval = pmax(pval.cvm,lim)/100
  Tn = -2*sum(log(pval))
  Tn2 = -2*sum(log(pval[card==2]))
  pvalTn = 100*(1-pchisq(Tn,2*m));
  pvalTn2 = 100*(1-pchisq(Tn2,d*(d-1)));
  pvalue=list(cvm=pval.cvm,Sn=pvalSn,Tn=pvalTn,Tn2=pvalTn2)
  stat=list(cvm=stats.cvm,Sn=Sn,Tn=Tn,Tn2=Tn2)
  subsets = vector(mode="character",length = m)

  for ( i in 1:m)
  {
    for(j in 1:d)
    {
      if(A[i,j])
      {
        subsets[i] = paste(subsets[i],as.character(j))
      }
    }

  }


  out0 = list(stat=stat, pvalue=pvalue, card=card,subsets=subsets)

  #end_time <- Sys.time()
  #print(end_time - start_time)
  if(graph){
    Dependogram(out0)

  }

   return(out0)
}






