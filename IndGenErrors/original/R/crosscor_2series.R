#'@title Cross-correlations for testing independence between the innovations of 2 series of same length
#'@description This function computes the cross-correlations  between x(t) and y(t-l), for l=-lag,.., lag, and also the combination (Wald's type)  of these statistics.
#'
#'@param x           Pseudo-observations (or residuals) of first series
#'@param y           Pseudo-observations (or residuals) of second series
#'@param lag         Maximum number of lags around 0
#'@param graph       Set to TRUE for a correlogram for all possible lags.
#
#'
#'@return \item{stat}{Cross-correlations for all lags}
#'@return \item{LB}{Sum of squares of cross-correlations}
#'@return \item{pvalue}{P-value of LB}
#'@return \item{subsets}{c(-lag:lag)}
#'@return \item{n}{length of the time series}
#'
#'@references Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'
#'@examples
#'data(gas)
#'outr <-crosscor_2series(gas$xres,gas$yres,3)
#'
#'@export
#'

crosscor_2series = function(x,y,lag,graph=TRUE)
{


m = 2*lag+1;
R12 = vector(mode = "double", length=m);


n = length(x);


out0 = .C("crosscor2d",as.double(x),
                  as.double(y),
                  as.integer(n),
                  as.integer(lag),
                  R120=double(m),
                  H12= double(1))

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
   R12[iid[i]] = out0$R120[i];
 }

out = list(stat=R12,LB=out0$H12,pvalue = 1-pchisq(out0$H12,m),subsets=c(-lag:lag),n=n)

if(graph)
{
  CrossCorrelogram(out,comb="(x,y)")
}



out
}







