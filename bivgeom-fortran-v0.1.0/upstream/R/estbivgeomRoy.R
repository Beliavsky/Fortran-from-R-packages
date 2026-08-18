estbivgeomRoy <-
function(x, y, method="LS")
{
est <- numeric(3)
# Least-Squares 
if (method=="LS")
{
d <- cbind(x,y)
z <- x*y
mod <- lm(log(S.n(d,d))~ -1 + x + y + z)
est <- exp(mod$coef)
}
# Maximum Likelihood
else if(method=="ML")
{
thetamin <- 1e-3 # theta1 and theta2 should lie in (0,1)! theta3 in (0,1]...
res <- mle2(minuslogRoy, method="L-BFGS-B", data=list(x=x,y=y), start=list(theta1=.5, theta2=.5,theta3=1), lower=c(theta1=thetamin,theta2=thetamin,theta3=thetamin),upper=c(theta1=1-thetamin,theta2=1-thetamin,theta3=1))
est <- res@fullcoef
}
# Method of Moments and Proportion
else if(method=="MMP")
{
hatt1 <- mean(x)/(1+mean(x))
hatt2 <- mean(y)/(1+mean(y))
hatt3 <- mean(x>=1&y>=1)/hatt1/hatt2
est <- c(hatt1, hatt2, hatt3)
}
# Method of Moments 1
else if(method=="MM1")
{
hatt1 <- mean(x)/(1+mean(x))
hatt2 <- mean(y)/(1+mean(y))
hatt3 <- (mean(x[x>0&y>0])-1)/hatt1/mean(x[x>0&y>0])
est <- c(hatt1, hatt2, hatt3)
}
# Method of Moments 2
else if(method=="MM2")
{
hatt1 <- mean(x)/(1+mean(x))
hatt2 <- mean(y)/(1+mean(y))
hatt3 <- (mean(y[y>0&x>0])-1)/hatt2/mean(y[y>0&x>0])
est <- c(hatt1, hatt2, hatt3)
}
# Method of Moments 3
else if(method=="MM3")
{
hatt1 <- mean(x)/(1+mean(x))
hatt2 <- mean(y)/(1+mean(y))
### solution of 2nd degree equation
hatt3 <- (mean(x[x>0&y>0])+mean(y[y>0&x>0])-2)/(hatt1*mean(x[x>0&y>0])+hatt2*mean(y[x>0&y>0]))
est <- c(hatt1, hatt2, hatt3)
}
# Method of Moments 4
else if (method=="MM4")
{
e10 <- mean(x[x>0 & y>0])
e01 <- mean(y[x>0 & y>0])
s <- e10 + e01
hatt1 <- mean(x)/(1+mean(x))
hatt2 <- mean(y)/(1+mean(y))
hatt3<- ((hatt1+hatt2)*(s-1) - sqrt(s*(s-2)*(hatt1-hatt2)^2+(hatt1+hatt2)^2))/(2*s*hatt1*hatt2)
est <- c(hatt1, hatt2, hatt3)
}
names(est) <- c("theta1","theta2","theta3")
return(est)
}
