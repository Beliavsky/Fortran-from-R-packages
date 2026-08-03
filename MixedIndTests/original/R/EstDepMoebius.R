#'@title Dependence measures and statistics for test of independence  between random variables
#'@description This function computes copula-based dependence measures for Moebius versions of Spearmans's rho, van der Waerden's coefficient, and Savage's coefficient, as well as their combination for tests of independence between random variables.
#'
#'@param x            Data matrix
#'@param trunc.level  Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.
#'@param graph        Set to TRUE if one wants the dependogram of P-values for the Moebius statistics
#'
#'
#'@return \item{stat}{List of statistics (spearman, vdw, savage) and test combinations Ln and Ln2 (only pairs)}
#'@return \item{pvalue}{P-values for the tests}
#'@return \item{cardA}{Cardinaly of the subsets for the Moebius statistics}
#'@return \item{subsets}{Subsets for the Moebius statistics}
#'
#'@references B.R Nasri & B.N. Remillard (2023). Tests of independence and randomness for arbitrary data using copula-based covariances
#'
#'@examples
#' x <- matrix(rnorm(250),ncol=5)
#' out <-EstDepMoebius(x,3)



EstDepMoebius =function(x,trunc.level=2,graph=FALSE){


  dim0 <- dim(x)
  n <- dim0[1]
  d <- dim0[2]
  v = c(2:(trunc.level))

  cA = sum(choose(d,v))



  out0 <- .C("Stat_A",
            as.double(x),
            as.integer(n),
            as.integer(d),
            as.integer(trunc.level),
            statS = double(cA),
            statG = double(cA),
            statE = double(cA),
            cardA = double(cA),
            Asets = double(d*cA),
            PACKAGE = "MixedIndTests"
  )

  m = length(out0$cardA)

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


  cA = length(out0$card)

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

  df1 = d*(d-1)/2;
  pvalue$Ln$spearman  = 100*(1-pchisq(Ln$spearman,cA));
  pvalue$Ln2$spearman = 100*(1-pchisq(Ln2$spearman,df1));
  pvalue$Ln$vdw  = 100*(1-pchisq(Ln$vdw,cA));
  pvalue$Ln2$vdw = 100*(1-pchisq(Ln2$vdw,df1));
  pvalue$Ln$savage  = 100*(1-pchisq(Ln$savage,cA));
  pvalue$Ln2$savage = 100*(1-pchisq(Ln2$savage,df1));





  out0 = list(stat=stat, pvalue=pvalue, card=card,subsets=subsets)


  if(graph){
    DependogramZ(out0,n)

  }

   return(out0)
}






