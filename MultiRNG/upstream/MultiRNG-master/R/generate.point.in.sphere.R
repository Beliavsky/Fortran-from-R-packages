generate.point.in.sphere=function(no.row,d){
  if ((d<2)|(floor(d)!=d)){
    stop("Dimension must be an integer whose value is at least 2!\n")
  }
  my.mat=matrix(0,no.row,d)
  if ((d==2)|(d>4)){
    for (i in 1:no.row){
      index=0
      while (index<1){
        u.mat=runif(d)*sample(c(-1,1), d, replace = TRUE)
        summ=sum(u.mat^2)
        my.mat[i,]=u.mat/sqrt(summ)
        index=1*(summ<=1)
      }
    }
  }
  if (d==3){
    for (i in 1: no.row){
      index=0
      while (index<1){
        u1=runif(1,-1,1)
        u2=runif(1,-1,1)
        s1=u1^2+u2^2
        w=(s1<=1)
        index=1*w
        my.mat[i,1][w]=2*u1[w]*sqrt(1-s1[w])
        my.mat[i,2][w]=2*u2[w]*sqrt(1-s1[w])
        my.mat[i,3][w]=1-2*s1[w]
      }
    }
  }
  if (d==4){
    for (i in 1: no.row){
      index=0
      while (index<1){
        u1=runif(1,-1,1)
        u2=runif(1,-1,1)
        u3=runif(1,-1,1)
        u4=runif(1,-1,1)
        s1=u1^2+u2^2
        s2=u3^2+u4^2
        w1=(s1<=1)
        w2=(s2<=1)
        index=1*(w1&w2)
        my.mat[i,1][w1&w2]=u1[w1&w2]
        my.mat[i,2][w1&w2]=u2[w1&w2]
        my.mat[i,3][w1&w2]=u3[w1&w2]*sqrt((1-s1[w1&w2])/s2[w1&w2])
        my.mat[i,4][w1&w2]=u4[w1&w2]*sqrt((1-s1[w1&w2])/s2[w1&w2])
      }
    }
  }
  my.mat
}
