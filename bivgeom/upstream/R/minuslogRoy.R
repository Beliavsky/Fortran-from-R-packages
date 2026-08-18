minuslogRoy <-
function(x,y,theta1=0.5,theta2=0.5,theta3=1)
{
if(theta3 < (theta1+theta2-1)/(theta1*theta2) | theta1*theta2*theta3==0 | theta3>1)
{s <- length(x)*99}
else
{
s <- -sum(log(dbivgeomRoy(x,y,theta1,theta2,theta3)))
}
s
}
