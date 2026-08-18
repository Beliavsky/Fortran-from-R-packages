draw.inv.wishart=function(no.row,d,nu,inv.sigma){
  if((no.row<1)|(floor(no.row)!=no.row)){
    stop("Number of replicates must be an integer whose value is at least 1!\n")
  }
  if((d<2)|(floor(d)!=d)){
    stop("Dimension must be an integer whose value is at least 2!\n")
  }
  if(nu<d){
    stop("Distribution is degenerate!\n")
    stop("Degrees of freedom should be greater than or equal to the dimension!\n")
  }
  if(nu==d+1){
    warning("Expectation does not exist!\n")
  }
  if(floor(nu)!=nu){
    stop("Degrees of freedom should be an integer!\n")
  }
  if((ncol(inv.sigma)!=d)|(nrow(inv.sigma)!=d)){
    stop("Inverse scale matrix is misspecified, dimension is wrong!\n")
  }
  if(min(eigen(inv.sigma)$values)<0){
    stop("Inverse scale matrix must be symmetric and positive definite!\n")
  }
  inv.wishart=draw.wishart(no.row,d,nu,solve(inv.sigma))
  # This function generates Wishart deviates in the form of rows.
  # To obtain the Inverted-Wishart matrix, convert each row to a matrix
  # where rows are filled first.
  inv.wishart
}