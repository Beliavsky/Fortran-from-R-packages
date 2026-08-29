# HEADER ####################################################
# This is file spam/R/spam_solve.R.                         #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



########################################################################
########################################################################
#
# Contains routines linked to solving spd linear systems. Namely:
#    chol, solve, backsolve, forwardsolve
#    determinant                (the later because it is based on chol)
#
# As well as associated S4 elements.
#
# The key element is a new class: "spam.chol.NgPeyton", the output
# of 'chol'
#
########################################################################
########################################################################


# lindx=  colindices
# xlindx= colpointers
# xlnz=   rowpointers
# snode=snmember
# xsuper=supernodes
# c(... nnztmp,cachesize)= memory
#### Summary of memory usage in cholstepwise:
# m
# nsuper
# nnzlmax
#
# Internally, we additionally assign several arrays (of length):
#  m: adj, colcnt
#  7*m+3: iwork

# iwork: is often used, anytime:


########################################################################

print.spam.chol.NgPeyton <- function(x,...) {
  nrow <- x@dimension[1]
  nnzR <- x@rowpointers[nrow+1]-1
  cat("(Upper) Cholesky factor of dimension ", nrow,
                "x", nrow, " with ",nnzR," (row-wise) nonzero elements (",printSize(x),").", sep = "", fill=TRUE)
  cat("    (The object is supposed to be used with: 'as.spam', 'backsolve', 'forwardsolve', etc.)\n",
      fill=TRUE)
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(NULL)
}

setMethod("show","spam.chol.NgPeyton", function(object) {
  nrow <- object@dimension[1]
  nnzR <- object@rowpointers[nrow+1]-1
  cat("(Upper) Cholesky factor of dimension ", nrow,
                "x", nrow, " with ",nnzR," (row-wise) nonzero elements (",printSize(object),").", sep = "", fill=TRUE)
  cat("    (The object is supposed to be used with: 'as.spam', 'backsolve', 'forwardsolve', etc.)\n",
      fill=TRUE)
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(NULL)
        })


"diag.of.spam.chol.NgPeyton" <- function(x, nrow, ncol)
  return( x@entries[x@rowpointers[-(x@dimension[1]+1)]])


setMethod("diag",    "spam.chol.NgPeyton", diag.of.spam.chol.NgPeyton)
setMethod("diag<-",  "spam.chol.NgPeyton", function(x) stop("diagonal cannot be changed on 'spam.chol.NgPeyton' object"))

setMethod("print",   "spam.chol.NgPeyton", print.spam.chol.NgPeyton)
setMethod("length",  "spam.chol.NgPeyton",function(x) x@rowpointers[x@dimension[1]+1]-1)
setMethod("length<-","spam.chol.NgPeyton",function(x,value) stop("length cannot be changed on 'spam.chol.NgPeyton' object") )
setMethod("dim",     "spam.chol.NgPeyton",function(x) x@dimension)
setMethod("dim<-",   "spam.chol.NgPeyton",function(x,value) stop("dimension cannot be altered on 'spam.chol.NgPeyton' object") )

# summary.spam.chol.NgPeyton is in file `summary.R`

setMethod("c","spam.chol.NgPeyton", function(x,...){
    nrow <- x@dimension[1]
    nnzR <- x@rowpointers[nrow+1]-1
    nsuper <- as.integer( length(x@supernodes)-1)
    if( getOption("spam.force64") || .format.spam(x)$package != "spam")
        SS <- .format64()
    else
        SS <- .format32

    xcolindices <- .C64('calcja',
                        SIGNATURE = rep(SS$signature, 7),

                        nrow,
                        nsuper,
                        x@supernodes,
                        x@colindices,

                        x@colpointers,
                        x@rowpointers,
                        xja = vector_dc( SS$type, nnzR),

                        INTENT=c("r", "r", "r", "r",
                                 "r", "r", "w"),
                        NAOK = getOption("spam.NAOK"),
                        PACKAGE = SS$package)$xja
    cx <- .C64("spamcsrdns",
               SIGNATURE = c(SS$signature, "double" , SS$signature, SS$signature, "double"),

               nrow = nrow,
               entries = x@entries,
               colindices = xcolindices,
               rowpointers = x@rowpointers,

               res = vector_dc( "double", nrow*nrow),

               INTENT = c("r", "r", "r", "r",
                          "w"),
               NAOK=getOption("spam.NAOK"),
               PACKAGE = SS$package)$res
  if (length( list(...)) < 1)
    return( cx)
  else
    c( cx,c(...))
})

as.spam.chol.NgPeyton <- function(x, eps = getOption("spam.eps")) {
    if( getOption("spam.force64") || .format.spam(x)$package != "spam")
        SS <- .format64()
    else
        SS <- .format32

  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  nrow <- x@dimension[1]
  nnzR <- x@rowpointers[nrow+1]-1
  nsuper <- length(x@supernodes)-1

  colindices <- .C64('calcja',
                     SIGNATURE=rep(SS$signature, 7),

                     nrow,
                     nsuper,
                     x@supernodes,
                     x@colindices,

                     x@colpointers,
                     x@rowpointers,
                     xja=vector(SS$type, nnzR), #!1#

                     INTENT=c("r", "r", "r", "r", "r", "r", "w"),
                     NAOK = getOption("spam.NAOK"),
                     PACKAGE = SS$package)$xja

  return(.newSpam(
    entries=x@entries,
    colindices=colindices,
    rowpointers=x@rowpointers,
    dimension=x@dimension
  ))
}



"as.matrix.spam.chol.NgPeyton" <- function(x,...){
    ## print("as.matrix.spam.chol.NgPeyton")
    nrow <- x@dimension[1]
    nnzR <- x@rowpointers[nrow+1]-1L
    ## newx <- new("spam")
    nsuper <-  length(x@supernodes)-1L
    ## xcolindices <- .Fortran('calcja',
    ##                         as.integer(nrow), as.integer(nsuper), as.integer(x@supernodes), as.integer(x@colindices), as.integer(x@colpointers), as.integer(x@rowpointers),
    ##                         xja=vector("integer",nnzR),
    ##                         NAOK = getOption("spam.NAOK"),PACKAGE = "spam")$xja
    if( getOption("spam.force64") || .format.spam(x)$package != "spam" )
        SS <- .format64()
    else
        SS <- .format32

    xcolindices <- .C64('calcja',
                        SIGNATURE = c(rep(SS$signature,7)),
                        nrow,
                        nsuper,
                        x@supernodes,
                        x@colindices,

                        x@colpointers,
                        x@rowpointers,
                        xja = vector_dc( SS$type, nnzR),

                        INTENT=c("r", "r", "r", "r",
                                 "r", "r", "w"),
                        NAOK = getOption("spam.NAOK"),
                        PACKAGE = SS$package)$xja
    return(array(.C64("spamcsrdns",
                      SIGNATURE = c(SS$signature, "double" , SS$signature, SS$signature,
                             "double"),
                 nrow = nrow,
                 entries = x@entries,
                 colindices = xcolindices,
                 rowpointers = x@rowpointers,

                 res = vector_dc( "double", nrow*nrow),  # dotCall takes care of size!

                 INTENT = c("r", "r", "r", "r",
                            "w"),
                 NAOK=getOption("spam.NAOK"),
                 PACKAGE = SS$package)$res,
               c(nrow,nrow))      # we preserve dimensions
         )
}

setMethod("as.spam","spam.chol.NgPeyton", as.spam.chol.NgPeyton)
setMethod("as.matrix","spam.chol.NgPeyton",as.matrix.spam.chol.NgPeyton)
setMethod("as.vector","spam.chol.NgPeyton",
          function(x){
            as.vector.spam(as.spam.chol.NgPeyton(x))
          })

########################################################################


setGeneric("backsolve", def = function(r, x, ...) standardGeneric("backsolve"),
           useAsDefault= function(r, x, ...) base::backsolve(r, x, ...))

setGeneric("forwardsolve", def = function(l, x, ...) standardGeneric("forwardsolve"),
           useAsDefault= function(l, x, ...) base::forwardsolve(l, x, ...))

# adapted from methods
#setGeneric("forwardsolve", function(l, x, k, upper.tri = FALSE, transpose = FALSE, ...)
#           standardGeneric("forwardsolve"),
#           useAsDefault = function(l, x, k = ncol(l), upper.tri = FALSE, transpose = FALSE, ...)
#                  base::forwardsolve(l, x, k = k, upper.tri = upper.tri, transpose = transpose, ... ),
#           signature = c("l", "x"))#, where = where)
##### setGenericImplicit("forwardsolve")#, restore=FALSE)


setMethod("chol","spam", chol.spam)
setMethod("solve","spam",solve.spam)
setMethod("chol2inv","spam", chol2inv.spam)
setMethod("chol2inv","spam.chol.NgPeyton", chol2inv.spam)

setMethod("backsolve","spam",              #signature(r="spam",x='ANY'),
          backsolve.spam)
setMethod("backsolve","spam.chol.NgPeyton",#signature(r="spam.chol.NgPeyton",x='ANY'),
          backsolve.spam, sealed=TRUE)
setMethod("forwardsolve","spam",               forwardsolve.spam)
setMethod("forwardsolve","spam.chol.NgPeyton", forwardsolve.spam)

# ?? sealed=TRUE ??  from help(setMethod):
#  'sealed' prevents the method being redefined, but should never be needed
#    when the method is defined in the source code of a package.

######################################################################


"ordering.default" <- function(x, inv=FALSE) stop('Operation not defined form this class')

setGeneric("ordering", function(x, inv=FALSE) standardGeneric("ordering"))

setMethod("ordering","spam.chol.NgPeyton",function(x,inv=FALSE)
          {
            if (inv) return(x@invpivot) else return(x@pivot) })

setMethod("ordering","matrix",function(x,inv=FALSE)
          {
            if (dim(x)[1]!=dim(x)[2])
              stop("ordering is defined for square matrices only")
            if(inv)return(dim(x)[1]:1) else return(1:dim(x)[1]) })

setMethod("ordering","spam",function(x,inv=FALSE)
          {
            if (dim(x)[1]!=dim(x)[2])
              stop("ordering is defined for square matrices only")
            if(inv) return(dim(x)[1]:1) else return(1:dim(x)[1]) })



########################################################################
#  force to spam matrices. Would not be required with inheritance

setMethod("image","spam.chol.NgPeyton",
          function(x,cex=NULL,...){
            image.spam(as.spam.chol.NgPeyton(x),cex=cex,...)
          })


setMethod("display","spam.chol.NgPeyton",
          function(x,...){
            display.spam(as.spam.chol.NgPeyton(x),...)
          })

setMethod("t","spam.chol.NgPeyton",
          function(x){
            t.spam(as.spam.chol.NgPeyton(x))
          })

setMethod("chol","spam.chol.NgPeyton",
          function(x){
           x
          })

