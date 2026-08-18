draw.wishart<-function(no.row,d,nu,sigma){
  if((no.row<1)|(floor(no.row)!=no.row)){
    stop("Number of replicates must be an integer whose value is at least 1!\n")
  }
  if((d<2)|(floor(d)!=d)){
    stop("Dimension must be an integer whose value is at least 2!\n")
  }
  if(nu<d){
    stop("Distribution is not proper !\n")
    stop("Degrees of freedom should be greater than or equal to the dimension!\n")
  }
  if(floor(nu)!=nu){
    stop("Degrees of freedom should be an integer!\n")
  }
  if((ncol(sigma)!=d)|(nrow(sigma)!=d)){
    stop("Scale matrix is misspecified, dimension is wrong!\n")
  }
  if(min(eigen(sigma)$values)<0){
    stop("Scale matrix must be symmetric and positive definite!\n")
  }
  wishart<-matrix(0,no.row,d^2)
  for (i in 1:no.row){
    alpha.i<-draw.d.variate.normal(nu,d,rep(0,d),sigma)
    wishart[i,]<-t(alpha.i)%*%alpha.i 
  }
  # This function generates Wishart deviates in the form of rows.
  # To obtain the Wishart matrix, convert each row to a matrix where
  # rows are filled first.
  wishart
}