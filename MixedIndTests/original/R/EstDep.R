#'@title Kendall's tau and Spearman's rho statistics  for testing independence between random variables
#'
#'@description This function computes the matrix of pairs of Kendall's tau and Spearman's rho statistics between random variables with arbitrary distributions.
#'
#'@param x           Data matrix
#'@param graph       Set to TRUE for a dependogram for all pairs of Kendall's taus and Spearman's rhos.
#
#'
#'@return \item{stat}{List of Kendall's tau and Spearman's rho statistics from multilinear copula, and test combinations LB}
#'@return \item{pvalue}{P-values for the tests statistics}
#'
#'@references Genest, Neslehova, Remillard & Murphy (2018). Testing for independence in arbitrary distributions
#'
#'@examples
#' x <- matrix(rnorm(500),ncol=10)
#' out <-EstDep(x)

EstDep = function(x,graph=FALSE)
{

dd = dim(x)
  n = dd[1]
  d = dd[2]

  tau = matrix(0,ncol=d,nrow=d)
  rho = matrix(0,ncol=d,nrow=d)
  pvaltau = matrix(0,ncol=d,nrow=d)
  pvalrho = matrix(0,ncol=d,nrow=d)

  LBrho=0
  LBtau=0
  tau[d,d]=1
  rho[d,d]=1
  for(j in 1:(d-1))
  {
    for(i in (j+1):d){
    out0 = stat_dep(x[,i],x[,j]);
    tau[i,j] = out0$tau;
    tau[j,i] = tau[i,j]
    rho[i,j] = out0$rho;
    rho[j,i] = rho[i,j]
    varH = out0$s2;
    zrho= sqrt(n)*rho[i,j];
    z = sqrt(n)*tau[i,j];
    ztau =  0.5*z/varH;
    pvaltau[i,j] = 200*pnorm(-abs(ztau))
    pvalrho[i,j] = 200*pnorm(-abs(zrho))
    pvalrho[j,i]=pvalrho[i,j]
    pvaltau[j,i]=pvaltau[i,j]
    LBrho = LBrho+zrho^2
    LBtau = LBtau+ztau^2
    }
    rho[j,j]=1
    tau[j,j]=1
    pvalrho[j,j]=100;
    pvaltau[j,j]=100;
  }


  df = d*(d-1)/2
  pvalLBrho = 100*(1-pchisq(LBrho,df));
  pvalLBtau = 100*(1-pchisq(LBtau,df));


  out0 = list(tau=tau,rho=rho,
             LBtau=LBtau,LBrho=LBrho,pvaltau=pvaltau,
             pvalrho=pvalrho,pvalLBtau=pvalLBtau,pvalLBrho=pvalLBrho)

  if(graph)
  {
    m = d*(d-1)/2


    subsets = vector(mode="character",length = m)
    pvalrho0 = vector(length=m)
    pvaltau0 = vector(length=m)
    i=1

    for(j in 1:d)
    {

      for(k in 1:d)
      {
        if(k>j)
        {
          subsets[i] = paste(as.character(j),as.character(k))
          pvalrho0[i] = pvalrho[j,k]
          pvaltau0[i] = pvaltau[j,k]
          i=i+1
        }

      }
    }

    outtau=list(subsets=subsets,pvalues=pvaltau0)
    Dependogram(outtau,stat="Kendall's tau")
    outrho=list(subsets=subsets,pvalues=pvalrho0)
    Dependogram(outrho,stat="Spearman's rho")

  }
  return(out0)
}
