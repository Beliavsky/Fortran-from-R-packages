draw.multivariate.laplace=function(no.row,d,gamma,mu,Sigma){
  if ((d<2)|(floor(d)!=d)){
    stop("Dimension must be an integer whose value is at least 2!\n")
  }
  if(gamma<=0){
    stop("Shape parameter must be positive!\n")
  }
  if((no.row<2)|(floor(no.row)!=no.row)){
    stop("Number of replicates must be an integer whose value is at least 2!\n")
  }
  if(d<2){
    stop("Dimension must be at least 2!\n")
  }
  if(length(mu)!=d){
    stop("Mean vector is misspecified, dimension is wrong!\n")
  }
  if((ncol(Sigma)!=d)|(nrow(Sigma)!=d)){
    stop("Variance-covariance matrix is misspecified, dimension is wrong!\n")
  }
  if(sum(Sigma!=t(Sigma))+min(eigen(Sigma)$values<=0)){
    stop("Variance-covariance matrix must be symmetric and positive definite!\n")
  }
  mul.laplace=matrix(0,no.row,d)
  for (i in 1:no.row){
    s=generate.point.in.sphere(1,d)
    mul.laplace[i,]=(rgamma(1,d,1)^(1/gamma))*t(chol(Sigma))%*%t(s)+mu
  }
  mul.laplace
}