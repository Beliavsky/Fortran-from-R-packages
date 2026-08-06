# Jacobian function for the computation of standard error of the coefficients of the conditional variance
ConJacob.stvar <-
function(a,lag,l,v,c)
{
Jv <- function(a)
{
S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S

S11 <- S[1:l,1:l] ; S12 <- S[1:l,(l+1):(l*(lag+1))] ; S22 <- S[(l+1):(l*(lag+1)),(l+1):(l*(lag+1))]
SS <- solve(S22)  ; Q <- SS/v ; B1<-SS%*%t(S12) ; s2 <- S11-S12%*%SS%*%t(S12) 

	c.var <- c(s2[1,1]*vech(SS)/(v+l*lag-2))
	var.coef <- c(v*s2[1,1]/(v+l*lag-2),c.var)

var.coef
}

return(list(Jv=jacobian(Jv,a)))

}
