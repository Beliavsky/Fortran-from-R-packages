Jacob.dlrm <-
function(a,lag,l,v,c)
{
J <- function(a)
{
S <- BlockTop(a[1:((lag+3)*(l+1)*l/2)],l,lag)$S

S11 <- S[1:1,1:1] ; S12 <- S[1:1,(1+1):ncol(S)] ; S22 <- S[(1+1):ncol(S),(1+1):ncol(S)]
SS <- solve(S22)  ; Q <- SS/v ; B1 <- SS%*%(S12) ; s2 <- S11-S12%*%SS%*%(S12) 

if(c!=0) {MU <- matrix(c(a[(((lag+3)*(l+1)*l/2)+1):(((lag+3)*(l+1)*l/2)+l*c)]),l,c)}

if(c!=0) {Delta0 <- MU[1:1,] - t(B1)%*%MU[c(2:l,rep(1:l,lag)),]}
if(c==0) {Delta0 <- 0}

      Cc <- c(Delta0,B1,vech(s2))

Cc
}

return(list(J=jacobian(J,a)))

}


