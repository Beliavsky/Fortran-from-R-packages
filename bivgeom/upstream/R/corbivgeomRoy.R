corbivgeomRoy <-
function(theta1,theta2,theta3)
{
# checking the parameters' feasibility
if(theta3 < (theta1+theta2-1)/(theta1*theta2) | theta1<=0 | theta1>=1 | theta2 <=0 | theta3 <=0 | theta3 >1) stop("parameters' values are not feasible")
# if theta3=1 we have independence and hence uncorrelation
if(theta3==1)
{
return(0)
}
else
{
# the correlation can be obtained as the sum of an infinite series
# we truncate this sum considering as the last term
# the maximum between twice the 0.99999 quantiles of X and Y
x <- 0:(max(2*qgeom(0.99999,1-theta1),2*qgeom(0.99999,1-theta2)))
# Covariance of X and Y
EXY <- theta2*sum((theta1*theta3)^(x+1)/(1-theta2*theta3^(x+1)))
# first moment and variance of X
mx <- theta1/(1-theta1)
vx <- theta1/(1-theta1)^2
# first moment and variance of Y
my <- theta2/(1-theta2)
vy <- theta2/(1-theta2)^2
(EXY-mx*my)/sqrt(vx*vy)
}
}
