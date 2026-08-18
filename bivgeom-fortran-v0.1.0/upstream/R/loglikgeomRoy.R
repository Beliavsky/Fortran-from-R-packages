loglikgeomRoy <-
function(par,x,y)
{
theta1<-par[1]
theta2<-par[2]
theta3<-par[3]
-sum(log(dbivgeomRoy(x,y,theta1,theta2,theta3)))
}
