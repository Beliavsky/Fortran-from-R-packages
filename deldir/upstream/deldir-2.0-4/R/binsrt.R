binsrtR <- function(x,y,rw) {
    n    <- length(x)
    ind  <- rep(0,n)
    rslt <- .Fortran("binsrt",
                     x=as.double(x),
                     y=as.double(y),
                     rw=as.double(rw),
                     n=as.integer(n),
                     ind=as.integer(ind),
                     rind=as.integer(ind),
                     tx=double(n),
                     ty=double(n),
                     ilst=integer(n),
                     PACKAGE="deldir"
                )
    list(x=rslt$tx,y=rslt$ty,ind=rslt$ind,rind=rslt$rind)
}
