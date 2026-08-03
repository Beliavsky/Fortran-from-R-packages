#'@title Cramer-von Mises Moebius statistics for testing independence between the innovations of 2 series of same length
#'@description This function computes the Cramer-von Mises statistics between x(t) and y(t-l), for l=-lag,.., lag, and also the combinations of the p-values of these statistics.
#'
#'@param x           Pseudo-observations (or residuals) of first series
#'@param y           Pseudo-observations (or residuals) of second series
#'@param lag         Maximum number of lags around 0
#'@param graph       Set to TRUE for a dependogram for all possible lags.
#
#'
#'@return \item{cvm}{Cramer-von Mises statistics for all lags}
#'@return \item{Wstat}{Sum of (unbiased) Cramer-von Mises statistics}
#'@return \item{Fstat}{Combination of p-values of the Cramer-von Mises statistics}
#'@return \item{pvalue}{List of p-values for the cvm, Wstat, and Fstat}
#'@references Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'
#'@examples
#'data(gas)
#'out <-cvm_2series(gas$xres,gas$yres,3)
#'
#'@export



cvm_2series = function(x,y,lag,graph=TRUE)
{


m = 2*lag+1;
T12 = vector(mode = "double", length=m);
PV12 = vector(mode = "double", length=m);

n = length(x);


out0 = .C("cvm2d",as.double(x),
                  as.double(y),
                  as.integer(n),
                  as.integer(lag),
                  T120=double(m),
                  PV120=double(m),
                  W12= double(1),
                  PVW12 = double(1),
                  F12= double(1),
                  PVF12 = double(1))

iid = vector(mode="double",length=m);

sum3 = 0;

iid = vector(mode= "integer",length = m);


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

for (i in 1:m)  {
     T12[iid[i]] = out0$T120[i];
     PV12[iid[i]] = out0$PV120[i];
      }


pvalue=list(cvm = PV12,Wstat = out0$PVW12,Fstat =out0$PVF12)
out = list(cvm=T12,Wstat=out0$W12, Fstat=out0$F12,pvalue=pvalue,subsets=c(-lag:lag))
if(graph)
{
  dependogram(out,stat="cvm")
}



out
}




