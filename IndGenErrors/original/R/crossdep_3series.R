#'@title Cross-dependence statistics for testing independence between the innovations of 3 series of same length
#'@description This function computes the cross-dependence for Spearman, van der Waerden and Savage dependence measures, for all lags = -lag2, .. lag2, for all pairs, and for pair of lags = (-lag3,-lag3),...(lag3,lag3) for the three series.
#'
#'@param x           Pseudo-observations (or residuals) of first series.
#'@param y           Pseudo-observations (or residuals) of second series.
#'@param z           Pseudo-observations (or residuals) of third series.
#'@param lag2        Maximum number of lags around 0 for pairs of series.
#'@param lag3        Maximum number of lags around 0 for the three series.
#'
#'@return \item{stat}{Cross-dependences for all lags  and for all subsets}
#'@return \item{H}{Sum of squares of cross-correlations for all subsets}
#'@return \item{pvalue}{P-value of LB for all subsets and H}
#'@return \item{n}{length of the time series}
#'
#'@references Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'@references Nasri & Remillard  (2024). Tests of independence and randomness for arbitrary data using copula-based covariances. JMVA, vol. 201, 105273.
#'@examples
#' #Romano-Siegel's example #
#'data(romano_ex)
#'outr = crossdep_3series(romano_ex$x,romano_ex$y,romano_ex$z,5,2)
#'CrossCorrelogram(outr$spearman$out123,"Savage for {1,2,3}",rot=90)
#'
#'@export

crossdep_3series <-function(x,y,z,lag2,lag3)
{

n = length(x);
y0 = c(y,y)
x0 = c(x,x)
z0 = c(z,z)

m1 = 2*lag2+1;
m2 = (2*lag3+1)*(2*lag3+1);

R12G = vector(mode = "double", length=m1);
R13G = vector(mode = "double", length=m1);
R23G = vector(mode = "double", length=m1);
R123G = vector(mode = "double", length=m2);
R12S = vector(mode = "double", length=m1);
R13S = vector(mode = "double", length=m1);
R23S = vector(mode = "double", length=m1);
R123S = vector(mode = "double", length=m2);
R12E = vector(mode = "double", length=m1);
R13E = vector(mode = "double", length=m1);
R23E = vector(mode = "double", length=m1);
R123E = vector(mode = "double", length=m2);

R120G = vector(mode = "double", length=m1);
R130G = vector(mode = "double", length=m1);
R230G = vector(mode = "double", length=m1);
R1230G = vector(mode = "double", length=m2);
R120S = vector(mode = "double", length=m1);
R130S = vector(mode = "double", length=m1);
R230S = vector(mode = "double", length=m1);
R1230S = vector(mode = "double", length=m2);
R120E = vector(mode = "double", length=m1);
R130E = vector(mode = "double", length=m1);
R230E = vector(mode = "double", length=m1);
R1230E = vector(mode = "double", length=m2);


ind0=c(1:n)


k=1
for(l in 0:lag2)
  {
  # A = {1,2}
  y1 = y0[l+ind0]
  z1 = z0[l+ind0]
  out0=MixedIndTests::EstDepMoebius(cbind(x,y1))
  R120S[k]= out0$stat$spearman
  R120G[k]= out0$stat$vdw
  R120E[k]= out0$stat$savage
  
  #  A = {1,3}  #
  out0=MixedIndTests::EstDepMoebius(cbind(x,z1))
  R130S[k]= out0$stat$spearman
  R130G[k]= out0$stat$vdw
  R130E[k]= out0$stat$savage
  
  #   A = {2,3}  #
  out0=MixedIndTests::EstDepMoebius(cbind(y,z1))
  R230S[k]= out0$stat$spearman
  R230G[k]= out0$stat$vdw
  R230E[k]= out0$stat$savage
  k=k+1
}

for(l in 1:lag2){
  x1 = x0[l+ind0]
  out0=MixedIndTests::EstDepMoebius(cbind(y,x1))
  R120S[k]= out0$stat$spearman
  R120G[k]= out0$stat$vdw
  R120E[k]= out0$stat$savage
  out0=MixedIndTests::EstDepMoebius(cbind(z,x1))
  R130S[k]= out0$stat$spearman
  R130G[k]= out0$stat$vdw
  R130E[k]= out0$stat$savage
  out0=MixedIndTests::EstDepMoebius(cbind(z,y1))
  R230S[k]= out0$stat$spearman
  R230G[k]= out0$stat$vdw
  R230E[k]= out0$stat$savage
  k=k+1
}

#  A = {1,2,3}  #

k=1

for(lag in 0:lag3)
 {
  for( l in 0:lag3)
  {
  y1 = y0[lag+ind0]
  z1 = z0[l+ind0]
  out0=MixedIndTests::EstDepMoebius(cbind(x,y1,z1),trunc.level=3)
  R1230S[k]= out0$stat$spearman[out0$card==3]
  R1230G[k]= out0$stat$vdw[out0$card==3]
  R1230E[k]= out0$stat$savage[out0$card==3]
  k=k+1
}
}

for(lag in 0:lag3)
  {
  for( l in 1:lag3)
  {
    y1 = y0[l+lag+ind0]
    x1 = x0[l+ind0]
    out0=MixedIndTests::EstDepMoebius(cbind(z,x1,y1),trunc.level=3)
    R1230S[k]= out0$stat$spearman[out0$card==3]
    R1230G[k]= out0$stat$vdw[out0$card==3]
    R1230E[k]= out0$stat$savage[out0$card==3]
    k=k+1
  }
}

for(lag in 1:lag3)
{
  for( l in 0:lag3)
  {
    x1 = x0[lag+ind0]
    z1 = z0[l+lag+ind0]
    out0=MixedIndTests::EstDepMoebius(cbind(y,x1,z1),trunc.level=3)
    R1230S[k]= out0$stat$spearman[out0$card==3]
    R1230G[k]= out0$stat$vdw[out0$card==3]
    R1230E[k]= out0$stat$savage[out0$card==3]
    k=k+1
  }
}

for(lag in 1:lag3)
{
  for( l in 1:lag)
  {
    x1 = x0[lag+ind0]
    z1 = z0[lag-l +ind0]
    out0=MixedIndTests::EstDepMoebius(cbind(y,x1,z1),trunc.level=3)
    R1230S[k]= out0$stat$spearman[out0$card==3]
    R1230G[k]= out0$stat$vdw[out0$card==3]
    R1230E[k]= out0$stat$savage[out0$card==3]
    k=k+1
  }
}

for(lag in 1:lag3)
{
  for( l in (lag+1):lag3)
  {
    x1 = x0[l+ind0]
    y1 = y0[l-lag +ind0]
    out0=MixedIndTests::EstDepMoebius(cbind(z,x1,y1),trunc.level=3)
    R1230S[k]= out0$stat$spearman[out0$card==3]
    R1230G[k]= out0$stat$vdw[out0$card==3]
    R1230E[k]= out0$stat$savage[out0$card==3]
    k=k+1
  }
}

#############################################




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




#############################################

for (i in 1:m1)  
    {
     R12G[iid[i]] = R120G[i];
     R13G[iid[i]] = R130G[i];
     R23G[iid[i]] = R230G[i];
     R12S[iid[i]] = R120S[i];
     R13S[iid[i]] = R130S[i];
     R23S[iid[i]] = R230S[i];
     R12E[iid[i]] = R120E[i];
     R13E[iid[i]] = R130E[i];
     R23E[iid[i]] = R230E[i];
    }



 for (i in 1:m2)  
 {  
    R123G[ind[i]] = R1230G[i];  
    R123S[ind[i]] = R1230S[i]; 
    R123E[ind[i]] = R1230E[i]; 
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

H12  = n*sum(R120S^2)
H13  = n*sum(R130S^2)
H23  = n*sum(R230S^2)
H123 = n*sum(R1230S^2)
H    = H12+H13+H23+H123
PVH12 = 1-pchisq(H12,m1);
PVH13 = 1-pchisq(H13,m1);
PVH23 = 1-pchisq(H23,m1);
PVH123 = 1-pchisq(H123,m2);
PVH = 1-pchisq(H,(3*m1+m2));

out12S  = list(stat=R12S,H12=H12,pvalue=PVH12,subsets=c(-lag2:lag2),n=n)
out13S  = list(stat=R13S,H13=H13,pvalue=PVH13,subsets=c(-lag2:lag2),n=n)
out23S  = list(stat=R23S,H23=H23,pvalue=PVH23,subsets=c(-lag2:lag2),n=n)
out123S = list(stat=R123S,H123=H123,pvalue=PVH123,subsets=subsets,n=n)

H12  = n*sum(R120G^2)
H13  = n*sum(R130G^2)
H23  = n*sum(R230G^2)
H123 = n*sum(R1230G^2)
H    = H12+H13+H23+H123
PVH12 = 1-pchisq(H12,m1);
PVH13 = 1-pchisq(H13,m1);
PVH23 = 1-pchisq(H23,m1);
PVH123 = 1-pchisq(H123,m2);
PVH = 1-pchisq(H,(3*m1+m2));
out12G  = list(stat=R12G,H12=H12,pvalue=PVH12,subsets=c(-lag2:lag2),n=n)
out13G  = list(stat=R13G,H13=H13,pvalue=PVH13,subsets=c(-lag2:lag2),n=n)
out23G  = list(stat=R23G,H23=H23,pvalue=PVH23,subsets=c(-lag2:lag2),n=n)
out123G = list(stat=R123G,H123=H123,pvalue=PVH123,subsets=subsets,n=n)

H12  = n*sum(R120E^2)
H13  = n*sum(R130E^2)
H23  = n*sum(R230E^2)
H123 = n*sum(R1230E^2)
H    = H12+H13+H23+H123
PVH12 = 1-pchisq(H12,m1);
PVH13 = 1-pchisq(H13,m1);
PVH23 = 1-pchisq(H23,m1);
PVH123 = 1-pchisq(H123,m2);
PVH = 1-pchisq(H,(3*m1+m2));
out12E  = list(stat=R12E,H12=H12,pvalue=PVH12,subsets=c(-lag2:lag2),n=n)
out13E  = list(stat=R13E,H13=H13,pvalue=PVH13,subsets=c(-lag2:lag2),n=n)
out23E  = list(stat=R23E,H23=H23,pvalue=PVH23,subsets=c(-lag2:lag2),n=n)
out123E = list(stat=R123E,H123=H123,pvalue=PVH123,subsets=subsets,n=n)

outS    = list(out12=out12S,out13=out13S,out23=out23S,out123=out123S,H=H,pvalue.H=PVH)
outG    = list(out12=out12G,out13=out13G,out23=out23G,out123=out123G,H=H,pvalue.H=PVH)
outE    = list(out12=out12E,out13=out13E,out23=out23E,out123=out123E,H=H,pvalue.H=PVH)

out     = list(spearman=outS, vdw = outG, savage = outE)

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





