divchain.default <- function (x,y,z,...) {
#
    if(missing(z)) {
        if(inherits(x,"ppp")) z <- x$marks
        else stop("Argument \"z\" was not supplied .\n")
    }
    z  <- factor(z)
    dd <- deldir(x,y,z=z,...)
    divchain(dd)
}
