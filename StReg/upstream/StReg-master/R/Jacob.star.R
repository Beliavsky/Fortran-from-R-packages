Jacob.star <-
function(a,lag,l,v,c)
{
J <- function(a)
{
S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S

S11 <- S[1:l,1:l] ; S12 <- S[1:l,(l+1):(l*(lag+1))] ; S22 <- S[(l+1):(l*(lag+1)),(l+1):(l*(lag+1))]
SS <- solve(S22)  ; Q <- SS/v ; B1 <- SS%*%(S12) ; s2 <- S11-S12%*%SS%*%(S12) 

MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)

Delta0 <- MU - t(B1)%*%(kronecker(diag(ones(lag)),MU)) 

      Cc <- c(Delta0,B1,vech(s2))

Cc
}

return(list(J=jacobian(J,a)))

}
