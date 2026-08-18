RelbivgeomRoy <- Vectorize(function(theta1, theta2, theta3)
{
if(theta3 < (theta1+theta2-1)/(theta1*theta2) | theta1<=0 | theta1>=1 | theta2 <=0 | theta3 <=0 | theta3 >1) stop("parameters' values are not feasible")
alpha <- 1e-5
# the reliabilty parameter can be obtained passing through the sum of an infinite series
# we truncate this sum considering as the last term
# the maximum between twice the (1-alpha) quantiles of X and Y
x <- 0:(max(2*qgeom(1-alpha,1-theta1),2*qgeom(1-alpha,1-theta2)))
#x<-0
#R <- 1-(1-theta2)*theta1*sum((theta1*theta2)^x*(theta3)^((x+1)^2))
R <- 1-sum(theta2^x*(theta1*theta3^x)^(x+1)*(1-theta3^(x+1)*theta2))
}
)