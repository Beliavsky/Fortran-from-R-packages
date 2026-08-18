SbivgeomRoy <- Vectorize(function(x,y,theta1,theta2,theta3)
{
if(theta3 < (theta1+theta2-1)/(theta1*theta2) | theta1<=0 | theta1>=1 | theta2 <=0 | theta3 <=0 | theta3 >1) stop("parameters' values are not feasible")
a <- ceiling(max(x,0))
b <- ceiling(max(y,0))
theta1^a*theta2^b*theta3^(a*b)
}
)
