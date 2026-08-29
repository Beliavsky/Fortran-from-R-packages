# HEADER ####################################################
# This is file spam/R/summary.R.                            #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################


summary.spam <- function(object,...) {
    nz <- length(object@entries)
    dens <- nz/prod(object@dimension)*100
    SS <- .format.spam(object, validate = FALSE)
    ans <- list(nnz=nz, density=dens, dim=object@dimension, format=SS)
    class(ans) <- "summary.spam"
    ans
}

summary.spam.chol.NgPeyton <- function(object,...) {
  nrow <- object@dimension[1]
  nnzR <- object@rowpointers[nrow+1]-1
  dens <- nnzR/(nrow^2)
  nnzc <- length(object@colindices)
  nnzA <- object@nnzA
  fill <- nnzR/((nnzA+nrow)/2)
  ans <-  list(nnzR=nnzR, nnzcolindices=nnzc,density=dens,fillin=fill, dim=c(nrow,nrow), nzzA=nnzA)
  class( ans) <- "summary.spam.chol.NgPeyton"
  ans
}

setMethod("summary","spam",summary.spam)
setMethod("summary", "spam.chol.NgPeyton", summary.spam.chol.NgPeyton)

print.summary.spam <- function(x, ...) {
    cat("Matrix object of class 'spam' of dimension ",x$dim[1],"x",
        x$dim[2],",\n",sep="")
    cat("    with ",x$nnz," (row-wise) nonzero elements.\n",sep="")
    cat("    Density of the matrix is ",signif(x$dens,3),"%.\n",sep="")
    cat(paste0("Class 'spam' (", x$format$name, ")\n"))
    invisible(x)
}

print.summary.spam.chol.NgPeyton <- function(x,...) {
  cat("(Upper) Cholesky factor of class 'spam.chol.NgPeyton' of dimension ", x$dim[1],
                "x", x$dim[1], " with ",x$nnzR," (row-wise) nonzero elements.", sep = "", fill=TRUE)
  cat("    Density of the factor is ", signif(x$dens * 100, 3),"%.\n", sep = "")
  cat("    Fill-in ratio is ", signif(x$fill, 3),"\n", sep = "")
  cat("    (Optimal argument for 'chol' is 'memory=list(nnzR=",x$nnzR,
            ifelse(x$nnzA<x$nnzc,paste(",nnzcolindices=",x$nnzc, sep = ""),""),")'.)\n", sep = "")
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(x)
}
