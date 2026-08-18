draw.dirichlet.multinomial=function(no.row,d,alpha,beta,N){
  if((no.row<1)|(floor(no.row)!=no.row)){
    stop("Number of replicates must be an integer whose value is at least 1!\n")
  }
  if((d<2)|(floor(d)!=d)){
    stop("Dimension must be an integer whose value is at least 2!\n")
  }
  if(length(alpha)!=d){
    stop("Shape vector is misspecified, dimension is wrong!\n")
  }
  if(min(alpha)<=0){
    stop("Shape vector cannot contain non-positive numbers!\n")
  }
  if(beta<=0){
    stop("Common scale parameter must be positive!\n")
  }
  if((N<2)|floor(N)!=N){
    stop("Size must be an integer whose value is at least 2!\n")
  }
  dirichlet=apply(draw.dirichlet(no.row,d,alpha,beta),2,mean)
  if(sum(dirichlet)!=1){
    dirichlet[d]=1-sum(dirichlet[1:d-1])
  }
  draws=draw.multinomial(no.row,d,dirichlet,N)
  draws
}