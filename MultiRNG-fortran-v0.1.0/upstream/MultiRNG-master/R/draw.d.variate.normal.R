draw.d.variate.normal=function(no.row,d,mean.vec,cov.mat){
  if ((no.row<1)|(floor(no.row)!=no.row)) {
    stop("Number of replicates must be an integer whose value is at least 1!\n")
  }
  if((d<2)|(floor(d)!=d)){
    stop("Dimension must be an integer whose value is at least 2!\n")
  } 
  if(length(mean.vec)!=d){
    stop("Mean vector is misspecified, dimension is wrong!\n")
  } 
  if((ncol(cov.mat)!=d)|(nrow(cov.mat)!=d)){
    stop("Variance-covariance matrix is misspecified, dimension is wrong!\n")
  }
  if(min(eigen(cov.mat)$values)<0){
    stop("Variance-covariance matrix must be symmetric and positive definite!\n")
  }
  z=matrix(rnorm(no.row*d),no.row,d)
  x=z%*%chol(cov.mat)+t(matrix (rep(mean.vec,no.row),nrow=d))
  x
}