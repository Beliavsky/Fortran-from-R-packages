#'@title Cross-correlations statistics for testing independence between the innovations of 3 series of same length
#'@description This function computes the cross-correlations for all lags = -lag2, .. lag2, for all pairs, and for pair of lags = (-lag3,-lag3),...(lag3,lag3) for the three series.
#'
#'@param x           Pseudo-observations (or residuals) of first series.
#'@param y           Pseudo-observations (or residuals) of second series.
#'@param z           Pseudo-observations (or residuals) of third series.
#'@param lag2        Maximum number of lags around 0 for pairs of series.
#'@param lag3        Maximum number of lags around 0 for the three series.
#'
#'@return \item{stat}{Cross-correlations for all lags  and for all subsets}
#'@return \item{H}{Sum of squares of cross-correlations for all subsets}
#'@return \item{pvalue}{P-value of stat for all subsets and H}
#'@return \item{n}{length of the time series}
#'
#'@references Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'
#'@examples
#' # Romano-Siegel's example #
#'data(romano_ex)
#'outr = crosscor_3series(romano_ex$x,romano_ex$y,romano_ex$z,5,2)
#'
#'@export





crosscor_3series <-function(x,y,z,lag2,lag3)
{

n = length(x);

m1 = 2*lag2+1;
m2 = (2*lag3+1)*(2*lag3+1);

R12 = vector(mode = "double", length=m1);
R13 = vector(mode = "double", length=m1);
R23 = vector(mode = "double", length=m1);
R123 = vector(mode = "double", length=m2);



out0 = .C("crosscor3d",
                  as.double(x),
                  as.double(y),
                  as.double(z),
                  as.integer(n),
                  as.integer(lag2),
                  as.integer(lag3),
                  R120=double(m1),
                  R130=double(m1),
                  R230=double(m1),
                  R1230=double(m2),
                  H12=double(1),
                  H13=double(1),
                  H23=double(1),
                  H123=double(1),
                  H=double(1))





m = 2*lag3+1;

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



out2 = out0[c(11:15)];

for (i in 1:m1)  {
     R12[iid[i]] = out0$R120[i];
     R13[iid[i]] = out0$R130[i];
     R23[iid[i]] = out0$R230[i];

        }



 for (i in 1:m2)  {  R123[ind[i]] =out0$R1230[i];        }






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

PVH12 = 1-pchisq(out2$H12,m1);
PVH13 = 1-pchisq(out2$H13,m1);
PVH23 = 1-pchisq(out2$H23,m1);
PVH123 = 1-pchisq(out2$H123,m2);
PVH = 1-pchisq(out2$H,(3*m1+m2));
out12 = list(stat=R12,H12=out2$H12,pvalue=PVH12,subsets=c(-lag2:lag2),n=n)
out13 = list(stat=R13,H13=out2$H13,pvalue=PVH13,subsets=c(-lag2:lag2),n=n)
out23 = list(stat=R23,H23=out2$H23,pvalue=PVH23,subsets=c(-lag2:lag2),n=n)
out123 = list(stat=R123,H123=out2$H123,pvalue=PVH123,subsets=subsets,n=n)
out=list(out12=out12,out13=out13,out23=out23,out123=out123,H=out0$H,pvalue.H=PVH)



#
# if(graph)
# {
#   CrossCorrelogram(out12,"{1,2}")
#   CrossCorrelogram(out13,"{1,3}")
#   CrossCorrelogram(out23,"{2,3}")
#   CrossCorrelogram(out123,"{1,2,3}",rot=90)
#
#
# }


out
}





