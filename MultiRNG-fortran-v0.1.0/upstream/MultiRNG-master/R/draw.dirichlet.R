draw.dirichlet=function(no.row,d,alpha,beta){
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
  mygamma=matrix(rgamma(no.row*d,alpha,beta),no.row,d,byrow=T)
  mybeta=matrix(0,no.row,d)
  for (i in 1:no.row){
    mybeta[i,]=mygamma[i,]/sum(mygamma[i,])
  }
  mybeta
}