#'@title Cross-dependences for testing independence between the innovations of 2 series of same length
#'@description This function computes the cross-dependence  between x(t) and y(t-l), for Spearman, van der Waerden and Savage dependence measures, for l=-lag,.., lag, and also the combination (Wald's type)  of these statistics.
#'
#'@param x           Pseudo-observations (or residuals) of first series
#'@param y           Pseudo-observations (or residuals) of second series
#'@param lag         Maximum number of lags around 0
#'@param graph       Set to TRUE for a correlogram for all possible lags.
#
#'
#'@return \item{stat}{Cross-dependences for all lags}
#'@return \item{H}{Sum of squares of cross-dependences}
#'@return \item{pvalue}{P-value of H}
#'@return \item{subsets}{c(-lag:lag)}
#'@return \item{n}{length of the time series}
#'
#'@references Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'@references Nasri & Remillard  (2024). Tests of independence and randomness for arbitrary data using copula-based covariances. JMVA, vol. 201, 105273.
#'@examples
#'data(gas)
#'outr <-crossdep_2series(gas$xres,gas$yres,3)
#'
#'@export
#'


crossdep_2series = function(x,y,lag,graph=TRUE)
{


m = 2*lag+1;
R12S = vector(mode = "double", length=m); #Spearman
R12G = vector(mode = "double", length=m); #van der Waerden
R12E = vector(mode = "double", length=m); #Savage
R120E=R12E
R120G=R12G
R120S=R12S

n = length(x);

y0 = c(y,y)
x0 = c(x,x)
ind = c(1:n)

k=1
for(l in 0:lag){
  y1 = y0[l+ind]
  out0=MixedIndTests::EstDepMoebius(cbind(x,y1))
  R120S[k]= out0$stat$spearman
  R120G[k]= out0$stat$vdw
  R120E[k]= out0$stat$savage
  k=k+1
}
for(l in 1:lag){
  x1 = x0[l+ind]
  out0=MixedIndTests::EstDepMoebius(cbind(y,x1))
  R120S[k]= out0$stat$spearman
  R120G[k]= out0$stat$vdw
  R120E[k]= out0$stat$savage
  k=k+1
}

iid = vector(mode= "integer",length = m);

sum3=0;
      for(l in 0:lag)
            {
              sum3 = sum3+1;
              iid[sum3] = lag+l+1;
            }

     for(l in 1:lag)
         {

          sum3 = sum3+1;
              iid[sum3] = lag-l+1;

         }


for (i in 1:m)
 {
   R12S[iid[i]] = R120S[i];
   R12G[iid[i]] = R120G[i];
   R12E[iid[i]] = R120E[i];
 }
H12=n*sum(R12G^2)
outG = list(stat=R12G, H=H12, pvalue = 1-pchisq(H12,m),subsets=c(-lag:lag),n=n)
H12=n*sum(R12E^2)
outE = list(stat=R12E, H=H12, pvalue = 1-pchisq(H12,m),subsets=c(-lag:lag),n=n)
H12=n*sum(R12S^2)
outS = list(stat=R12S, H=H12, pvalue = 1-pchisq(H12,m),subsets=c(-lag:lag),n=n)

if(graph)
{
  CrossCorrelogram(outS,comb="Spearman:(x,y)")
  CrossCorrelogram(outG,comb="van der Waerden: (x,y)")
  CrossCorrelogram(outE,comb="Savage: (x,y)")
}



out = list(spearman=outS,vdw=outG,savage=outE)
out
}







