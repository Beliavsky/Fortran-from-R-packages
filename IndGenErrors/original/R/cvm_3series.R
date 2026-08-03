#'@title Cramer-von Mises Moebius statistics for testing independence between the innovations of 3 series of same length
#'@description This function computes the Cramer-von Mises statistics between x(t), y(t-l2), z(t-l3), for l2=-lag2,.., lag2, l3=-lag3,.., lag3,and also the combinations of these statistics.
#'
#'@param x           Pseudo-observations (or residuals) of first series.
#'@param y           Pseudo-observations (or residuals) of second series.
#'@param z           Pseudo-observations (or residuals) of third series.
#'@param lag2        Maximum number of lags around 0 for pairs of series.
#'@param lag3        Maximum number of lags around 0 for the three series.
#
#'
#'@return \item{cvm}{Cramer-von Mises statistics for all lags and for all subsets}
#'@return \item{Wstat}{Sum of (unbiased) Cramer-von Mises statistics for all subsets}
#'@return \item{Fstat}{Combination of p-values of the Cramer-von Mises statistics}
#'@return \item{pvalue}{List of p-values for the cvm, Wstat, and Fstat}
#'@references Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'
#'@examples
#'set.seed(1)
#'x0 = rnorm(100); y = rnorm(100); z = rnorm(100);
# x = abs(x0)*sign(y*z);
# out = cvm_3series(x,y,z,5,2)
#'
#'@export

cvm_3series <-function(x,y,z,lag2,lag3)
{

n = length(x);

m1 = 2*lag2+1;
m2 = (2*lag3+1)*(2*lag3+1);

T12 = vector(mode = "double", length=m1);
PV12 = vector(mode = "double", length=m1);

T13 = vector(mode = "double", length=m1);
PV13 = vector(mode = "double", length=m1);

T23 = vector(mode = "double", length=m1);
PV23 = vector(mode = "double", length=m1);

T123 = vector(mode = "double", length=m2);
PV123 = vector(mode = "double", length=m2);


out0 = .C("cvm3d",
                  as.double(x),
                  as.double(y),
                  as.double(z),
                  as.integer(n),
                  as.integer(lag2),
                  as.integer(lag3),
                  T120=double(m1),
                  T130=double(m1),
                  T230=double(m1),
                  T1230=double(m2),
                  PV120=double(m1),
                  PV130=double(m1),
                  PV230=double(m1),
                  PV1230=double(m2),
                  W12=double(1),
                  W13=double(1),
                  W23=double(1),
                  W123=double(1),
                  W=double(1),
                  PVW12=double(1),
                  PVW13=double(1),
                  PVW23=double(1),
                  PVW123=double(1),
                  PVW=double(1),
                  F12=double(1),
                  F13=double(1),
                  F23=double(1),
                  F123=double(1),
                  Fstat=double(1),
                  PVF12=double(1),
                  PVF13=double(1),
                  PVF23=double(1),
                  PVF123=double(1),
                  PVF=double(1))

################# Reordering of the T stats

sum3 = 0;

iid = vector(mode= "integer",length = m1);


      for(l in 0:lag2)
            {
              sum3 = sum3+1;
              iid[sum3] = lag2+l+1;
            }

     for(l in 1:lag2)
         {

          sum3 = sum3+1;
              iid[sum3] = lag2-l+1;

         }


m = 2*lag3+1;

 sum3 = 0;

ind = vector(mode= "integer",length = m2);


    for(lag in 0:lag3)
        {
          for(l in 0:lag3)
            {
              sum3 = sum3+1;
              ind[sum3] = (lag+lag3)*m +lag3+l+1;
         }
      }

      for(lag  in 0:lag3)
      {
         for(l in 1:lag3)
         {

          sum3 = sum3+1;
              ind[sum3] = (lag+lag3)*m +lag3-l+1;

         }
      }

      for(lag in 1:lag3)
      {
         for(l in 0:lag3)
         {
            sum3 = sum3+1;
              ind[sum3] = (-lag+lag3)*m +lag3+l+1;


         }
      }

      for(lag in 1:lag3)
      {
         for(l in 1:lag)
         {
           sum3 = sum3+1;
              ind[sum3] = (-lag+lag3)*m +lag3-l+1;

         }
       }

       for(lag in 1:(lag3-1))
        {
         for(l in (lag+1):lag3)
         {
            sum3 = sum3+1;
              ind[sum3] = (-lag+lag3)*m +lag3-l+1;
         }
      }



out2 = out0[c(15:34)];

for (i in 1:m1)  {
     T12[iid[i]] = out0$T120[i];
     PV12[iid[i]] = out0$PV120[i];
     T13[iid[i]] = out0$T130[i];
     PV13[iid[i]] = out0$PV13[i];
     T23[iid[i]] = out0$T23[i];
     PV23[iid[i]] = out0$PV23[i];

      }

 for (i in 1:m2)
  {  T123[ind[i]] =out0$T1230[i];
     PV123[ind[i]] =out0$PV1230[i];
  }






subsets = vector(mode="character",length = m2);

sum3=0;

 for ( i in -lag3:lag3)
   {
      for( j in -lag3:lag3)
          {
            sum3 = sum3+1;
            subsets[sum3] = paste('{',i,',',j,'}');
          }
   }

pvalue12=list(cvm = PV12,Wstat = out0$PVW12,Fstat =out0$PVF12)
pvalue13=list(cvm = PV13,Wstat = out0$PVW13,Fstat =out0$PVF13)
pvalue23=list(cvm = PV23,Wstat = out0$PVW23,Fstat =out0$PVF23)
pvalue123=list(cvm = PV123,Wstat = out0$PVW123,Fstat =out0$PVF123)
out12 = list(cvm=T12,Wstat=out0$W12, Fstat=out0$F12,pvalue=pvalue12,subsets=c(-lag2:lag2))
out13 = list(cvm=T13,Wstat=out0$W13, Fstat=out0$F13,pvalue=pvalue13,subsets=c(-lag2:lag2))
out23 = list(cvm=T23,Wstat=out0$W23, Fstat=out0$F23,pvalue=pvalue23,subsets=c(-lag2:lag2))
out123 = list(cvm=T123,Wstat=out0$W123, Fstat=out0$F123,pvalue=pvalue123,subsets=subsets)
out00 = list(Wstat=out0$W,Fstat=out0$Fstat,pvalue.Wstat=out0$PVW,pvalue.Fstat=out0$PVF)
out = list(out12=out12,out13=out13,out23=out23,out123=out123,outall=out00);





out
}





