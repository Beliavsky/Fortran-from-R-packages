#This computes the parameters of the conditional distribution using the paramters of the joint distribution
Par.stvar <-
function(a,S,MU,l,lag,v,T,M,c)
{

S11 <- S[1:l,1:l] ; S12 <- S[1:l,(l+1):(l*(lag+1))] ; S22 <- S[(l+1):(l*(lag+1)),(l+1):(l*(lag+1))]
SS <- solve(S22)  ; Q <- SS/v ; B1 <- SS%*%t(S12) ; s2 <- S11-S12%*%SS%*%t(S12) 

Delta <- matrix(nrow=(T-lag),ncol=l)
for(i in 1:(T-lag))
{
if (c!=0) {Delta[i,] <- M[1:l,i] - t(B1)%*%M[(l+1):nrow(M),i]}
}

if (c!=0) {Delta0 <- MU - t(B1)%*%(kronecker(diag(ones(lag)),MU)) }

if (c==0) {Delta <- 0 ; Delta0 <- 0}

c.var <- (s2[1,1]*vech(SS)/(v+l*lag-2))
var.coef <- c(v*s2[1,1]/(v+l*lag-2),c.var)

return(list(Delta=Delta,Delta0=Delta0,Q=Q,s2=s2,S=S,B1=B1,var.coef=var.coef))

}
