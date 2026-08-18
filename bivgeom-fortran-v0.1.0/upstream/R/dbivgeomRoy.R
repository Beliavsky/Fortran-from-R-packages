dbivgeomRoy <- Vectorize(function(x,y,theta1,theta2,theta3)
{
if(theta3 < (theta1+theta2-1)/(theta1*theta2) | theta1<=0 | theta1>=1 | theta2 <=0 | theta3 <=0 | theta3 >1) stop("parameters' values are not feasible")
if(ceiling(x)!=floor(x) | ceiling(y)!=floor(y) | x<0 | y<0)
return(0)
else
theta1^x*theta2^y*theta3^(x*y)*(1-theta1*theta3^y-theta2*theta3^x+(theta1*theta2*theta3)*theta3^(x+y))
}
)
