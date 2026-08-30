tileArea <- function(x,y,rw) {
   n <- length(x)
   eps <- sqrt(.Machine$double.eps)
   area <- 0
   for(i in 1:n) {
       ip <- if(i==n) 1 else i+1
       tmp <- .Fortran(
                  "stoke",
                  x1=as.double(x[i]),
                  y1=as.double(y[i]),
                  x2=as.double(x[ip]),
                  y2=as.double(y[ip]),
                  rw=as.double(rw),
                  area=double(1),
                  s1=double(1),
                  eps=as.double(eps),
                  PACKAGE="deldir"
              )
        area <- area+tmp[["area"]]*tmp[["s1"]]
    }
    area
}
