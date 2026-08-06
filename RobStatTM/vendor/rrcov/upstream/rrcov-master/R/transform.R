transform_ilr <- function(x)
{
    x.ilr <- x[, -ncol(x), drop=FALSE]
    colnames(x.ilr) <- paste("X", 1:ncol(x.ilr), sep="")
    rownames(x.ilr) <- NULL
    D <- ncol(x)

    for (i in 1:ncol(x.ilr)){
        x.ilr[,i]=sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop=FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
        ##x.ilr[,i]=sqrt((D-i)/(D-i+1))*log(apply(as.matrix(x[,(i+1):D]), 1, gm)/(x[,i]))	
    }

	return(-x.ilr)
}
