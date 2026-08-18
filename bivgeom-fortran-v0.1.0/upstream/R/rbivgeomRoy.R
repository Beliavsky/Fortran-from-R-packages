rbivgeomRoy <-
function(n, theta1, theta2, theta3)
{
# checking the parameters' feasibility
if(theta3 < (theta1+theta2-1)/(theta1*theta2) | theta1<=0 | theta1>=1 | theta2 <=0 | theta3 <=0 | theta3 >1) stop("parameters' values are not feasible")
# threshold used for root search
eps <- 0.00001
u <- runif(n)
v <- runif(n)
# sampling x (unconditionally, from a geometric distribution with parameter 1-theta1)
x <- qgeom(u, 1-theta1)
xv <- cbind(x,v)
# numerically find the real root of the equation F_{y|x} - u = 0
res <- apply(xv, 1,function(z) uniroot(function(y) FyxbivgeomRoy(y, theta1=theta1, theta2=theta2, theta3=theta3, x=z[1])-z[2], lower=-1+eps,upper=2*qgeom(max(1-eps,z[2]),1-theta2))$root)
# sampling y, conditionally on x
y <- ceiling(res) # take the smallest integer >= to the real number obtained
cbind(x,y)
}
