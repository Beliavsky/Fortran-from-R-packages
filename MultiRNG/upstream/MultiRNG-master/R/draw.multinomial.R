draw.multinomial=function(no.row,d,theta,N){
  if((no.row<1)|floor(no.row)!=no.row){
    stop("Number of replicate samples must be integer whose value is at least 1!\n")
  }
  if((d<1)|floor(d)!=d){
    stop("Dimension must be integer whose value is at least 2!\n")
  }
  if (length(theta)!=d){
    stop("Length of the parameter vector does not match the dimension!\n")
  }
  if (min(theta)<0){
    stop("Parameter vector contains negative values!\n")
  }
  if (sum(theta)!=1){
    stop("Sum of probabilities must be 1!\n")
  }
  if((N<2)|floor(N)!=N){
    stop("Size must be an integer whose value is at least 2!\n")
  }
  mult=matrix(0,no.row,d)
  mytheta=theta
  for (r in 1:no.row){
    theta=mytheta
    size=N
    mult[r,1]=rbinom(1,size,theta[1])
    for (j in 2:(d-1)){
      size=N-sum(mult[r,1:(j-1)])
      theta[j]=theta[j]/sum(theta[j:d])
      mult[r,j]=rbinom(1,size,theta[j])
    }
    mult[r,d]=N-sum(mult[r,1:(d-1)])
  }
  mult
}