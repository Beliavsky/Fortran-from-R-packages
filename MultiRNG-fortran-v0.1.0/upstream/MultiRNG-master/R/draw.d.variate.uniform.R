draw.d.variate.uniform=function(no.row,d,cov.mat){
  if((no.row<1)|(floor(no.row)!=no.row)){
    stop("Number of subjects must be an integer whose value is at least 1!\n")
  }
  if(d<2){
    stop("Number of variables must be at least 2!\n")
  } 
  if((ncol(cov.mat)!=d)|(nrow(cov.mat)!=d)){
    stop("Variance-covariance matrix is misspecified, dimension is wrong!\n")
  } 
  if(sum(cov.mat!=t(cov.mat))+min(eigen(cov.mat)$values<=0)){
    stop("Variance-covariance matrix must be symmetric and positive definite!\n")
  }
  draw=draw.d.variate.normal(no.row,d,mean.vec=rep(0,d),cov.mat)
  x=pnorm(draw)
  x
}