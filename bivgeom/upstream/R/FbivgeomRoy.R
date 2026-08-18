FbivgeomRoy <- Vectorize(function(x,y,theta1,theta2,theta3)
{
if(x<0 | y<0) return(0)
else
1-theta1^(floor(x)+1)-theta2^(floor(y)+1)+theta1^(floor(x)+1)*theta2^(floor(y)+1)*theta3^((floor(x)+1)*(floor(y)+1))
}
)
