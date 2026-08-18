loc.min=function(my.mat,d){
  w=is.matrix(my.mat)
  if (w==F){
    stop("This is not a matrix!\n")
  }
  if (nrow(my.mat)!=ncol(my.mat)){
    stop("This is not a square matrix!\n")
  }
  n=nrow(my.mat)
  my.vec=as.vector(t(my.mat))
  my.vec[my.vec==0]=999
  my.index=min((1:length(my.vec))[my.vec==min(my.vec)])
  row.index=floor((my.index-1)/n)+1
  col.index=my.index-d*floor((my.index-1)/n)
  c(row.index,col.index)
}