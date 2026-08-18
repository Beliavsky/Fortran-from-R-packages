draw.multivariate.hypergeometric=function(no.row,d,mean.vec,k){
  if((no.row<1)|(floor(no.row)!=no.row)){
    stop("Number of replicates must be an integer whose value is at least 1!\n")
  }
  if((d<2)|(floor(d)!=d)){
    stop("Dimension must be an integer whose value is at least 2!\n")
  }
  if(length(mean.vec)!=d){
    stop("Number of items are misspecified, dimension is wrong!\n")
  }
  if(min(mean.vec)<=0){
    stop("Number of items vector cannot contain non-positive numbers!\n")
  }
  if(sum(floor(mean.vec)!=mean.vec)>0){
    stop("Number of items vector cannot contain non-integer numbers!\n")
  }
  if((k<=0)|(floor(k)!=k)){
    stop("Number of items to be sampled must be a positive integer!\n")
  }
  if(k>sum(mean.vec)){
    stop("Number of items to be sampled cannot be greater than the total items!\n")
  }
  x=matrix(0,no.row,d)
  tot.m=sum(mean.vec)
  myk=k
  for (i in 1:no.row){
    summ=tot.m
    k=myk
    for (j in 1:(d-1)){
      x[i,j]=rhyper(1,mean.vec[j],summ-mean.vec[j],k)
      k=k-x[i,j]
      summ=summ-mean.vec[j]}
    x[i,d]=k
  }
  x
}