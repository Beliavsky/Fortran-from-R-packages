#'@title Dependence measures and statistics for test of randomness for a univariate time series
#'
#'@description This function computes copula-based dependence measures for Moebius versions of Spearmans's rho, van der Waerden's coefficient, and Savage's coefficient, as well as their combination for tests of randomness for p consecutive values Y(1), ..., Y(p).
#'
#'@param y            Time series
#'@param p            Number of consecutive observations
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#'@param graph        Set to TRUE if one wants the dependogram of P-values for the Moebius statistics
#'
#'
#'@return \item{stat}{List of statistics (spearman, vdw, savage) and test combinations Ln and Ln2 (only pairs)}
#'@return \item{pvalue}{P-values for the tests}
#'@return \item{card}{Cardinaly of the subsets for the Moebius statistics}
#'@return \item{subsets}{Subsets for the Moebius statistics}
#'
#'@references B.R Nasri & B.N. Remillard (2023). Tests of independence and randomness for arbitrary data using copula-based covariances
#'
#'@examples
#' y<- SimAR1Poisson(c(5,0.2),100)
#' out <- EstDepSerialMoebius(y,4,4)



EstDepSerialMoebius =function(y,p,trunc.level=2,graph=FALSE){


 # start_time <- Sys.time()

  v = c(1:(trunc.level-1))
  cA = sum(choose((p-1),v))
  m = 2^(p-1)-1
  n = length(y)

 # void Stat_A_serial(double *y, int *n, int *d, int *trunc, double *statS, double *statG, double *statE,  double *cardA, double *Asets)

  out0 <- .C("Stat_A_serial",
            as.double(y),
            as.integer(n),
            as.integer(p),
            as.integer(trunc.level),
            statS = double(cA),
            statG = double(cA),
            statE = double(cA),
            cardA = double(cA),
            Asets = double(p*cA),
            PACKAGE = "MixedIndTests"
  )

  #end_time <- Sys.time()
  #print(end_time - start_time)
  Asets = matrix(out0$Asets,ncol=p)
  cardA=out0$cardA

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

  cA = length(out0$card)

  stat=list();
  stat$spearman = out0$statS[ind$ix];
  stat$vdw      = out0$statG[ind$ix];
  stat$savage   = out0$statE[ind$ix];


  Ln=list();
  Ln2=list();


  Ln$spearman = n*sum(stat$spearman^2);
  Ln$vdw = n*sum(stat$vdw^2);
  Ln$savage = n*sum(stat$savage^2);
  Ln2$spearman = n*sum(stat$spearman[card==2]^2);
  Ln2$vdw = n*sum(stat$vdw[card==2]^2);
  Ln2$savage = n*sum(stat$savage[card==2]^2);

  stat$Ln  = Ln;
  stat$Ln2 = Ln2;

  pvalue = list();

  pvalue$Ln$spearman  = 100*(1-pchisq(Ln$spearman,cA));
  pvalue$Ln2$spearman = 100*(1-pchisq(Ln2$spearman,p-1));
  pvalue$Ln$vdw  = 100*(1-pchisq(Ln$vdw,cA));
  pvalue$Ln2$vdw = 100*(1-pchisq(Ln2$vdw,p-1));
  pvalue$Ln$savage  = 100*(1-pchisq(Ln$savage,cA));
  pvalue$Ln2$savage = 100*(1-pchisq(Ln2$savage,p-1));


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
    DependogramZ(out0,n)

  }

   return(out0)
}






