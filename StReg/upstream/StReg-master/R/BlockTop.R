## This function produces matrix in block Toeplitz form using the parameters/initial values supplied as "a" vector

BlockTop <-
function(a,l,lag)
{
C <- array(dim=c(l,l,lag+3))
for(x in 1:(lag+3))
{
C[,,x] <- c(xpnd(a[((l*(l+1)/2)*(x-1)+1):(x*(l+1)*l/2)]))
}

f <- c(C)
A <- array(matrix(c(f,f[1:(l*l)],f),nrow=l,ncol=(2*(lag+4)-1)*l),dim=c(l,l,2*(lag+4)-1))
MM <- matrix(nrow=(lag+4)*(l),ncol=(lag+4)*(l))

for(i in 1:(lag+4))
{
for(j in 1:(lag+4))
{
MM[(l*i-(l-1)):(i*l),(l*j-(l-1)):(j*l)] <- A[,,j+(i-1)]
}
}

S <- MM[1:((1+lag)*l),]%*%t(MM[1:((1+lag)*l),])
return(list(S=S))
}
