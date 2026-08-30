findNewInOld <- function(xnew,xold,ynew,yold,tolerance=sqrt(.Machine$double.eps)) {
#
# Check that x and y lengths match.
if(length(xnew) != length(ynew) | length(xold) != length(yold))
    stop("Mismatch in lengths of x and y sequences.\n")

# If either the old or new sequence of points is empty, there
# is nothing to find.
if(!(length(xnew) & length(xold))) return(numeric(0))

# Set holder for found indices. 
ind  <- numeric(length(xnew))

# Scale up tolerance if the "old" points are large in modulus.
sfac <- max(1,mean(sqrt(xold^2+yold^2)))
tol  <- tolerance*sfac

# Search for matches.
for(i in seq_along(xnew)) {
    for(j in seq_along(xold)) {
        xok <- abs(xnew[i] - xold[j]) < tol
        yok <- abs(ynew[i] - yold[j]) < tol
        if(xok & yok) {
            ind[i] <- j
            break
        }
    }
}
ind
}
