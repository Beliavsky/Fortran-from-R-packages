# HEADER ####################################################
# This is file spam/R/definitions.R.                        #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################


# We check that
# 1.  #colindices equals #entries
# 2.	class(entries) is "numeric"
# 3.	rowpointers[1] = 1 and rowpointers[MAX] = #entries + 1
# 4.	dimension[1] ==  #rowpointers - 1
#
# --- Test in C ---
# c1.  rowpointers must be increasing
# c2.  No NaN/NA/Inf in rowpointer, colindices or entries
# c3.  colindices are in domain [1, dimension[2]]
#
# Warning: if class(colindices) != class(rowpointers)
validate_spam <- function(object) {
  # 1
  if(length(object@entries) != length(object@colindices))
    return("Length of slot 'entries' does not equal to length of slot 'colindices'.")

  # 2
  if(!is.double(object@entries))
    return("Slot 'entries' must be of type double (not integer).")

  # 3
  if(object@rowpointers[1] != 1)
    return("Slot 'rowpointers' must be have 1 on the first position.")
  if(object@rowpointers[length(object@rowpointers)] != length(object@entries) + 1)
    return("Last value of rowpointers is incorrect.")

  # 4
  if(length(object@rowpointers) != object@dimension[1] + 1)
    return("Slot 'rowpointers' has not expected length.")


  # --- For performance issues, the following tests should be implemented in C ---

  # c1
  result <- any(head(object@rowpointers, n=-1) > object@rowpointers[-1], na.rm=TRUE)
  if(result == TRUE)
    return("Slot 'rowpointers' is not increasing.")

  # c2
  if(any(!is.finite(object@rowpointers)))
    return("Slot 'rowpointers' contains non-finite elements.")
  if(any(!is.finite(object@colindices)))
    return("Slot 'colindices' contains non-finite elements.")

  # raise a warning if any non-finite element is in entries:
  if(!getOption("spam.NAOK") && any(!is.finite(object@entries)))
    warning("Slot 'entries' contains non-finite elements.")

  if(any(object@colindices < 1) | any(object@colindices > object@dimension[2]))
    return("Slot 'colindices' contains element outside the domain [1, ncol].")

  return(TRUE)
}



# Start with the two classes. Remaining part contains `spam` only material.
setClass(
  Class = "spam",
  slots = c(entries="numeric",
            colindices="numeric",
            rowpointers="numeric",
            dimension="numeric"
    ),
  prototype = prototype(entries=0,
                        colindices=1,
                        rowpointers=c(1,2),
                        dimension=c(1,1)
    ),
  validity = validate_spam
)

setClass("spam.chol.NgPeyton",
         representation(entries="numeric",      colindices="integer",
                        colpointers="integer",  rowpointers="integer",
                        dimension="integer",
                        pivot="integer",        invpivot="integer",
                        supernodes="integer",   snmember="integer",
                        memory="integer",       nnzA="integer")
         )

# lindx=  colindices
# xlindx= colpointers
# xlnz=   rowpointers
# snode=snmember
# xsuper=supernodes
# c(... nnztmp,cachesize)= memory
########################################################################



# Internal function to create new spam or spam64 object, depending on the size.
.newSpam <- function(entries=0,
                     colindices=1,
                     rowpointers=NULL,
                     dimension=c(1,1),
                     force64=getOption("spam.force64")
                     )
{
    if(is.null(rowpointers)) {
        rowpointers <- c(1,rep_len64(2,dimension[1]))
    }

    if(force64 || length(colindices) > 2147483647 ||
       length(colindices) > 2147483647 || max(dimension) > 2147483647 ){
        newx <- new("spam")
        slot(newx,"entries",check=FALSE) <- entries
        entries <- NULL
        slot(newx,"colindices",check=FALSE) <- as.numeric(colindices)
        colindices <- NULL
        slot(newx,"rowpointers",check=FALSE) <- as.numeric(rowpointers)
        rowpointers <- NULL
        slot(newx,"dimension",check=FALSE) <- as.numeric(dimension)
        return(newx)
    } else {
        newx <- new("spam")
        slot(newx,"entries",check=FALSE) <- entries
        entries <- NULL
        slot(newx,"colindices",check=FALSE) <- as.integer(colindices)
        colindices <- NULL
        slot(newx,"rowpointers",check=FALSE) <- as.integer(rowpointers)
        rowpointers <- NULL
        slot(newx,"dimension",check=FALSE) <- as.integer(dimension)
        return(newx)
    }
}


# Create a string containing the matrix size with appropriate unit.
# not planned to be exported.
# If length(x)==1, this represents the number of doubles (or other size).
printSize <- function(x, size=8, digits=2L) {
    mat <- if (length(x)==1) 8*as.numeric(x) else  as.numeric( object.size(x))
    unitpower <- max(min( floor( log(mat, 1024)),3),1)
    paste(round(mat/1024^unitpower, digits=digits), c("Kb","Mb","Gb")[unitpower])
}


.print.spam <- function(x, digits = 2, rowpointer = FALSE, minimal = TRUE, # show row pointers in left margin
              zerosymbol = "'", ...){

    if (minimal) {
      print(x@entries, ...)
      return(invisible(NULL))
    }

    # never true if used via print.spam
    if(is.matrix(x)) {
      x <- as.spam(x) }
    digits <- round(digits)

    stopifnot(is.spam(x), is.numeric(digits), length(digits)==1,
              !rowpointer || (rowpointer && (is.spam(x) || is(x, "spam.chol.NgPeyton"))), # rowpointes are only available for spam
              is.character(zerosymbol), nchar(zerosymbol) == 1)

    if (length(x@entries) == 0 || ((length(x@entries) == 1) && (x@entries[1] == 0))) {
        cat("zero-matrix. class: spam (", .format.spam(x)$name, "). dim: ",  x@dimension[1],
              "x", x@dimension[2], ".", sep = "", fill = TRUE)

        return(invisible(NULL))
    }

    tmp <- round(x@entries, digits = digits)
    x@entries <- ifelse(tmp == 0 & x@entries < 0, 0, tmp) # because sprintf() prints "-0"

    if(digits == 0 || any(grep("\\.", format(x@entries, scientific = FALSE)))) {
        digits <- min(digits, max(nchar(gsub("(.*)(\\.)|([0]*$)","", x@entries))))
    } else {
        digits <- 0 }

    width <- max(nchar(as.character(round(abs(x@entries))))) #left side of comma
    if(any(round(x@entries,digits = digits) < 0)) { # minus sign
        width <- width+1 }
    if(digits > 0) {
        width <- width+1+digits }
    if(width < 3 && ncol(x) >= 10) { #minus sign or two digits for col number
        width <- width+1 }

    if(nrow(x) > getOption("spam.printsize") || ncol(x)*(width+1) > getOption("width")[1]) {
        cat("class: spam (", .format.spam(x)$name, "). dim: ", x@dimension[1], "x",
            x@dimension[2], ". density: ", signif(100*length(x@entries)/prod(x@dimension),3),
            "%.\nrow-wise nonzero elements:",   sep = "", fill = TRUE)
        if(length(x@entries)>10)
            cat(round(x@entries[1:10], digits), "...\n")
        else
            cat(round(x@entries, digits), "\n")

        return(invisible(NULL))
    }

    headerstring <- paste("Matrix of dimension ",x@dimension[1],"x",
                          x@dimension[2], " with (row-wise) nonzero elements", sep ="")
    if (rowpointer)
      headerstring <- paste(headerstring, ",\n number left of the vertical bars denote the respective rowpointer:\n", sep = "")
    else
      headerstring <- paste(headerstring, ":\n", sep = "")

    cat(headerstring)
    cat("    ", formatC(1:ncol(x), width = width), "\n", sep = " ")
    cat("    ", rep("-", ncol(x)*(width+1)), "\n" , sep = "")
    zerostring <- paste(c(rep(" ",width-1), zerosymbol), collapse="")
    rowtmp <- rep(zerostring, x@dimension[2])
    for(i in 1:x@dimension[1]){
        if(rowpointer)
            cat(formatC(x@rowpointers[i], width = 3), "| ", sep = "")
        else
            cat(formatC(i, width = 3), "| ", sep = "")
        if(x@rowpointers[i]==x@rowpointers[i+1]){
            cat(rowtmp, "\n")
            next
        }
        row <- rowtmp
        coli <- x@colindices[x@rowpointers[i]:(x@rowpointers[i+1]-1)]
        row[coli] <- sprintf(paste0("%", width, ".", digits,"f", collapse=""),
                             x@entries[x@rowpointers[i]:(x@rowpointers[i+1]-1)],
                             width = digits+2)
        cat(row, "\n")
    }
    cat(paste0("class: spam (", .format.spam(x)$name, ")\n"))
}


print_nnzpos <- function(x, ...){
  stopifnot(is.spam(x) || is(x, "spam.chol.NgPeyton"))
  if(prod(x@dimension) < getOption("spam.printsize"))
    stop("dim too large for nnzpos=TRUE")
  x@entries <- seq_along(x@entries)
  spam::print.spam(x, digits = 0, minimal = FALSE, ...)
}


print.spam <- function(x, ...) {
  SS <- .format.spam(x, validate = FALSE)
  if (prod(x@dimension) < getOption("spam.printsize")) {
    print(as.matrix.spam(x), ...)
  } else {
    if ( (length(x@entries)==1) && (x@entries[1]==0)) {
      cat("Zero matrix of dimension ",x@dimension[1],"x",
                x@dimension[2],".\n",sep="", fill=TRUE)
    } else {
      .print.spam(x, ...)
    }
  }

  cat(paste0("Class 'spam' (", SS$name, ")\n"))
  invisible(NULL)
}


setMethod("show","spam",  function(object) print.spam(x = object))
setMethod("print","spam", print.spam)

setMethod("length<-","spam",function(x,value) stop("operation not allowed on 'spam' object") )
setMethod("length","spam",function(x) length(x@entries) )
## equivalent to x@rowpointers[x@dimension[1]+1]-1
## however this resulted in a strange error:
## R version 3.2.1 (2015-06-18) -- "World-Famous Astronaut"
## Platform: x86_64-unknown-linux-gnu (64-bit)
## > require(spam)
## Spam version 1.4-0 (2016-08-29) is loaded.
## > gctorture(TRUE)
## > length(spam(1))
## Error: unimplemented type 'integer' in 'coerceToInteger'


setMethod("c","spam", function(x,...){
    if(getOption("spam.force64") || .format.spam(x)$package != "spam")
        SS <- .format64()
    else
        SS <- .format32

    dimx <- x@dimension

    cx <- .C64("spamcsrdns",
               SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature,
                             "double"),

               nrow = x@dimension,
               entries = x@entries,
               colindices = x@colindices,
               rowpointers = x@rowpointers,

               res = vector_dc("double", prod(dimx)),

               INTENT = c("r", "r", "r", "r",
                   "rw"),
               NAOK = getOption("spam.NAOK"),
               PACKAGE = SS$package
               )$res

    if (length( list(...)) < 1){
        return( cx)
    }else{
        return( c( cx, c(...)))
    }
})

########################################################################
# diag and derivatives
"diag.spam" <- function(x=1, nrow, ncol)  {
  ## structure of function according an older version of base:::diag
  if (is.spam(x)) return( diag.of.spam( x, nrow, ncol))

  if (is.array(x) && length(dim(x)) != 1)
    stop("first argument is array, but not matrix.")

  if (missing(x))
    n <- as.integer(nrow)
  else if (length(x) == 1 && missing(nrow) && missing(ncol)) {
    n <- as.integer(x)
    x <- 1
  }
  else n <- length(x)
  if (!missing(nrow))
    n <- as.integer(nrow)
  if(missing(ncol))
    ncol <- n
  p <- as.integer(ncol)

  m <- min(n, p)

  ###CheckMe: are the following two lines efficient?
  entries <- vector("double", m)
  entries[1:m] <- x
  return(.newSpam(
      entries=entries,
        colindices=1:m,
        rowpointers=c(1:m,rep(m+1,n+1-m)),
        dimension=c(n,p)
    ))
}

spam_diag <- function(x=1, nrow, ncol)  diag.spam(x=x, nrow, ncol)

"diag<-.spam" <-  function(x, value) {
    if(getOption("spam.force64") || .format.spam(x)$package != "spam")
        SS <- .format64()
    else
        SS <- .format32

    nrow <- x@dimension[1]
    minrc <- min( x@dimension)
    if (length(value)==1)
        value <- rep(value,minrc)
    else if (length(value)!=minrc)
        stop("replacement diagonal has wrong length")

                                        #!7#
    if(2147483647 < minrc + length(x@entries))
        SS <- .format64()

    z <- .C64("setdiagmat",
              SIGNATURE=c(SS$signature, SS$signature, "double", SS$signature,
                          SS$signature, "double", SS$signature),

              nrow = nrow,
              n = minrc,
              a = c(x@entries, vector("double", minrc)),
              ja = c(x@colindices, vector(SS$type, minrc)),
              ia = x@rowpointers,
              diag = value,
              iw = vector_dc(SS$type, nrow), #!8#

              INTENT=c("r", "r", "rw", "rw", "rw", "r", "r"),
              NAOK= getOption("spam.NAOK"),
              PACKAGE = SS$package)


    nz <- z$ia[nrow+1]-1
    entries <- z$a
    colindices <- z$ja
    rowpointers <- z$ia
    rm(z)

    length(entries) <- nz
    length(colindices) <- nz

    return(.newSpam(
        entries=entries,
        colindices=colindices,
        rowpointers=rowpointers,
        dimension=x@dimension
    ))
}

"diag.spam<-" <- get("diag<-.spam")

diag.of.spam <- function(x, nrow, ncol) {
    if(getOption("spam.force64") || .format.spam(x)$package != "spam")
        SS <- .format64()
    else
        SS <- .format32

    len <- min(x@dimension)

    result <- .C64("getdiag",
                   SIGNATURE=c("double", SS$signature, SS$signature, SS$signature, "double"),

                   a = x@entries,
                   colindices =  x@colindices,
                   rowpointers = x@rowpointers,
                   len = len,
                   diag = vector_dc("double",len),

                   INTENT=c("r", "r", "r", "r", "w"),
                   NAOK = getOption("spam.NAOK"),
                   PACKAGE = SS$package
                   )
    return(result$diag)
}


setMethod("diag","spam",diag.of.spam)
setMethod("diag<-","spam",get("diag<-.spam"))


########################################################################

t.spam <- function(x){
# TODO: Check if we need to switch to 64-bit (in case of too many rows/cols)
    if(getOption("spam.force64") || .format.spam(x)$package != "spam")
        SS <- .format64()
    else
        SS <- .format32

    dimx <- x@dimension
    nz <- x@rowpointers[dimx[1]+1]-1
    z <- .C64("transpose",
              SIGNATURE=c(SS$signature, SS$signature,
                          "double", SS$signature, SS$signature,
                          "double", SS$signature, SS$signature),
              n = dimx[1],
              m = dimx[2],

              a = x@entries,
              ja = x@colindices,
              ia = x@rowpointers,

              entries = vector_dc("double",nz),
              colindices = vector_dc(SS$type,nz),
              rowpointers = vector_dc(SS$type,dimx[2]+1),

              INTENT = c("r", "r",
                         "r", "r", "r",
                         "w", "w", "rw"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)


    return(.newSpam(
        entries=z$entries,
        colindices=z$colindices,
        rowpointers=z$rowpointers,
        dimension=dimx[2:1]
    ))
}

setMethod("t","spam",t.spam)


########################################################################


"is.spam" <- function(x) is(x,"spam")

"as.spam" <- function(x,  eps = getOption("spam.eps"))
    stop("coercion not defined form this class")

"spam" <- function(x, nrow = 1, ncol = 1, eps = getOption("spam.eps"))
    stop("argument 'x' should be of mode 'numeric' (or 'spam')")

as.spam.spam <- function(x, eps = getOption("spam.eps"))  {
    force64 <- getOption("spam.force64")

    if(force64)
        SS <- .format64()
    else
        SS <- .format.spam(x)

    if (eps < .Machine$double.eps)
        stop("'eps' should not be smaller than machine precision", call.=FALSE)
    dimx <- x@dimension

    z <- .C64("cleanspam",
              SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature, "double"),

              nrow=dimx[1],
              entries=x@entries,
              colindices=x@colindices,
              rowpointers=x@rowpointers,
              eps=eps,

              INTENT=c("r", "rw", "rw", "rw", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package
              )
    nz <- z$rowpointers[dimx[1]+1]-1
    if(nz==0) {
        return(.newSpam(
            dimension=dimx,
            force64 = force64
        ))
    }

    return(.newSpam(
        entries=z$entries[1:nz],
        colindices=z$colindices[1:nz],
        rowpointers=z$rowpointers[1:(dimx[1]+1)],
        dimension=dimx,
        force64=force64
    ))
}

"cleanup" <- function(x, eps = getOption("spam.eps")) {
  if (is.spam(x)) as.spam.spam(x,eps) else x
}

as.spam.matrix <- function(x, eps = getOption("spam.eps")) {
    force64 <- getOption("spam.force64")
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)

    dimx <- dim(x)
    nonz <- abs(x) > eps
    nz <- sum(nonz, na.rm = TRUE) + sum(!is.finite(nonz))

    #only zero entries
    if(nz==0) {
        return(.newSpam(
            dimension=dimx,
            force64=force64
        ))
    }

    # EXPLICIT-STORAGE-FORMAT: Depending on the length of x, use 32 or 64-bit:
    SS <- if( force64 || nz > 2147483647 || max(dimx) > 2147483647 ){
              .format64()
          }else{
              .format32
          }

    z <- .C64("spamdnscsr",
              SIGNATURE=c(SS$signature, SS$signature, "double", SS$signature,
                          "double", SS$signature, SS$signature, "double"),

              nrow = dimx[1],
              ncol = dimx[2],
              x = x,
              dimx[1],

              entries = vector_dc("double",nz),
              colindices = vector_dc(SS$type,nz),
              rowpointers = vector_dc(SS$type,dimx[1]+1),
              eps = eps,

              INTENT=c("r", "r", "r", "r",
                       "w", "w", "w", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package
              )

    return(.newSpam( entries=z$entries,
                    colindices=z$colindices,
                    rowpointers=z$rowpointers,
                    dimension=dimx,
                    force64 = force64))
}

"as.spam.numeric" <- function(x, eps = getOption("spam.eps")) {
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)

    force64 <- getOption("spam.force64")

    poselements <- (abs(x)>eps)
    if (any(!is.finite(x))) {
        poselements[!is.finite(x)] <- TRUE
    }
    lx <- length(x)
    nz <- sum(as.numeric(poselements))

    # empty matrix
    if (identical(nz,0)){
        return(
            .newSpam(
                # rowpointers = c(1,rep_len64(2,lx)),
                dimension = c(lx,1),
                force64 = force64
            )
        )
    }

    return(
        .newSpam(
            entries = as.double(x[poselements]),
            colindices = rep_len64(1, nz),
            rowpointers = cumsum(c(1, poselements)),
            dimension = c(lx,1),
            force64 = force64
        )
    )
}

as.spam.dist <- function(x, eps = getOption("spam.eps")) {
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
    if (any(!is.finite(x))) {
        warning("'NA/NaN/Inf' coerced to zero")
        x[!is.finite(x)] <- 0
    }
    dimx <- attr(x,"Size")
    size <- dimx*(dimx-1)/2

    force64 <- getOption("spam.force64")

    if(force64 || dimx^2 > 2147483647)
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("disttospam",
              SIGNATURE = c(SS$signature, "double",
                            "double", SS$signature, SS$signature,
                            "double"),

              nrow=dimx, #r
              x=as.vector(x,"double"), #r

              entries=vector_dc("double",size), #w
              colindices=vector_dc("integer",size), #w
              rowpointers=vector_dc("integer",dimx+1), #w

              eps=eps, #r

              INTENT = c("r", "r",
                         "w", "w", "w",
                         "r"),
              NAOK=getOption("spam.NAOK"),
              PACKAGE = SS$package
              )
    nz <- z$rowpointers[dimx+1]-1
    if(nz==0) return( .newSpam(
                              # rowpointers = c( 1, rep( 2, dimx)),
                              dimension = c( dimx, dimx),
                              force64 = force64))


    return(.newSpam(
        entries=z$entries[1:nz],
        colindices=z$colindices[1:nz],
        rowpointers=z$rowpointers[1:(dimx+1)],
        dimension=c(dimx,dimx),
        force64 = force64
    ))
}

"as.spam.list" <-  function(x, eps = getOption("spam.eps"))
    spam.list(x,eps)

spam.numeric <- function(x, nrow = 1, ncol = 1, eps = getOption("spam.eps"))  {
    force64 <- getOption("spam.force64")
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
    if (any(!is.finite(x))) {
        warning("'NA/NaN/Inf' coerced to zero")
        x[!is.finite(x)] <- 0
    }
    lenx <- length( x)
    if (missing(nrow))
        nrow <- ceiling( lenx/ncol)
    else if (missing(ncol))
        ncol <- ceiling( lenx/nrow)
    dimx <- c(nrow, ncol)
    if (lenx != prod(nrow,  ncol)) {
        if(lenx==1 && abs(x)<eps) {
            return(.newSpam(
                dimension=dimx,
                force64=force64
            ))
        }
        else if(prod(nrow,ncol)%%lenx!=0)
            warning("ncol*nrow indivisable by length(x)")

        x <- rep_len64(x, prod(nrow,ncol))
    }


    dimx <- c(nrow, ncol)

    nz <- sum(as.numeric(abs(x) > eps))

    if(is.na(nz) | is.nan(nz))
        stop("NA or NaN in matrix.")

    if(nz==0) {
        return(.newSpam(
            dimension=dimx,
            force64=force64
        ))
    }

    # EXPLICIT-STORAGE-FORMAT: Depending on the length of x, use 32 or 64-bit:
    if(force64 || nz > 2147483647 || max(dimx) > 2147483647 )
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("spamdnscsr",
              SIGNATURE=c(SS$signature, SS$signature, "double", SS$signature,
                          "double", SS$signature, SS$signature, "double"),

              nrow = dimx[1],
              ncol = dimx[2],
              x = x,
              dimx[1],

              entries = vector_dc("double",nz),
              colindices = vector_dc(SS$type,nz),
              rowpointers = vector_dc(SS$type,dimx[1]+1),
              eps = eps,

              INTENT=c("r", "r", "r", "r",
                       "w", "w", "w", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package
              )

    entries <- z$entries
    colindices <- z$colindices
    rowpointers <- z$rowpointers
    z <- NULL

    return( .newSpam( entries = entries,
                     colindices = colindices,
                     rowpointers = rowpointers,
                     dimension = dimx,
                     force64 = force64))
}

setOldClass(c("dist", "numeric"))


setGeneric("as.spam")
setMethod("as.spam","spam",   as.spam.spam)
setMethod("as.spam","matrix", as.spam.matrix)
setMethod("as.spam","numeric",as.spam.numeric)
setMethod("as.spam","dist",   as.spam.dist)

setGeneric("spam")
setMethod("spam","numeric",spam.numeric)
setMethod("spam","spam",as.spam)



triplet <- function(x, tri=FALSE){
# inverse of spam.list
  dimx <- dim(x)
  if (length(dimx)!=2) stop("'x' should be a matrix like object of dim 2")
  if (is.spam(x)) {
    return(c({ if (tri) list(i=rep(1:dimx[1],diff(x@rowpointers)),
                  j= x@colindices) else list(indices=cbind(rep(1:dimx[1],diff(x@rowpointers)),
                  x@colindices) ) }, list(values=x@entries)
             )
           )
  } else {
    return(c({ if (tri) list(i=rep.int(1:dimx[1],dimx[2]),
                             j=rep.int(1:dimx[2],rep.int(dimx[1],dimx[2]))) else
              list(indices=cbind(rep.int(1:dimx[1],dimx[2]),
                  rep.int(1:dimx[2],rep.int(dimx[1],dimx[2]))))
             } , list(values=c(x))
             )
           )
  }
}

########################################################################

as.matrix.spam <- function(x, ...) {
    if( getOption("spam.force64") || .format.spam(x)$package == "spam64")
        SS <- .format64()
    else
        SS <- .format32

    dimx <- x@dimension
    m <- .C64("spamcsrdns",
              SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature,
                            "double"),

              nrow = dimx[1],
              entries = x@entries,
              colindices = x@colindices,
              rowpointers = x@rowpointers,

              res = vector_dc("double", prod(dimx)),

              INTENT = c("r", "r", "r", "r",
                         "rw"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)$res

    dim(m) <- dimx
    return(m)
}


setMethod("as.matrix","spam",as.matrix.spam)

"as.vector.spam" <- function(x, mode="any"){
    if(getOption("spam.structurebased")) {
        return( as.vector( x@entries, mode=mode))
    } else {
        inefficiencywarning( gettextf("This `as.vector` coercion operation may be inefficient"), prod(dim(x)))
        return( as.vector( as.matrix(x), mode=mode))
    }
}

setMethod("as.vector","spam",as.vector.spam)


########################################################################


"complement.spam" <- function(x){
# this function returns the structure of the zeros of the spam object x.
    force64 <- getOption("spam.force64")

    nrow <- x@dimension[1]
    ncol <- x@dimension[2]
    nnz <- x@rowpointers[nrow+1]-1
    nz <- prod(nrow,ncol) - nnz

                                        # we work through special cases
    if(nz == 0){
        ## print("1")
        return( spam(0, nrow, ncol))
    }
    if(nnz == 1 && x@entries == 0){
        ## print("2")
        return( spam(1, nrow, ncol))
    }

                                        # normal case, proceed to efficient function
    ## print("3")
    if( force64 || .format.spam(x)$package == "spam64" || nz > 2147483647 || nrow > 2147483646)
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("notzero",
              SIGNATURE=c(SS$signature, SS$signature,
                          SS$signature, SS$signature, SS$signature, SS$signature,
                          SS$signature, SS$signature),

              ja = x@colindices,
              ia = x@rowpointers,

              nrow = nrow,
              ncol = ncol,
              nnz = nnz,
              nz = nz,

              colindices = vector_dc( SS$type, nz),
              rowpointers = vector_dc( SS$type, nrow+1),

              INTENT=c("r", "r",
                       "r", "r", "r", "r",
                       "w", "w"),
              NAOK=getOption("spam.NAOK"),
              PACKAGE= SS$package
              )
    ## newx <- new("spam")
    ## slot(newx,"entries",check=FALSE) <- rep.int(1.0,nz)
    ## slot(newx,"colindices",check=FALSE) <- z$colindices
    ## slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
    ## slot(newx,"dimension",check=FALSE) <- x@dimension
    return(.newSpam(
        entries = rep_len64(1.0, nz),
        colindices = z$colindices,
        rowpointers = z$rowpointers,
        dimension = x@dimension,
        force64 = force64
    ))
}


# ASSIGNING:
##########################################################################################

# as S3subsetting causes problems, we eliminate this...
#"[<-.spam" <- function (x, rw, cl,value) {#cat("qq");
#                                          assign.spam(x,rw,cl,value) }


setReplaceMethod("[", signature(x = "spam",
			 i = "missing", j = "missing", value = "ANY"),
	  function (x, i, j, value) {#cat("mm");
                                     assign.spam(x,1:x@dimension[1],1:x@dimension[2],value)})


setMethod("[<-",signature(x="spam",i="vector",j="vector", value = "spam"),
	  function (x, i, j, value)
          {#### FIXME Highly inefficient!
            inefficiencywarning( "Performing inefficient replacement...", prod(dim(value)))
            as.spam(assign.spam(x,i,j,as.matrix(value)))} )
setMethod("[<-",signature(x="spam",i="missing",j="vector", value = "spam"),
	  function (x, i, j, value)
          {#### FIXME Highly inefficient!
            inefficiencywarning( "Performing inefficient replacement...", prod(dim(value)))
            as.spam(assign.spam(x,1:x@dimension[1],j,as.matrix(value)))} )
setMethod("[<-",signature(x="spam",i="vector",j="missing", value = "spam"),
	  function (x, i, j, value)
          {#### FIXME Highly inefficient!
            inefficiencywarning( "Performing inefficient replacement...", prod(dim(value)))
            as.spam(assign.spam(x,i,1:x@dimension[2],as.matrix(value)))} )


setMethod("[<-",signature(x="spam",i="vector",j="missing", value = "ANY"),
	  function (x, i, j, value) {#cat(i);
            assign.spam(x,i,1:x@dimension[2],value)} )

setMethod("[<-",signature(x="spam",i="vector",j="vector", value = "ANY"),
	  function (x, i, j, value) {#cat("vv");
            assign.spam(x,i,j,value)} )

setMethod("[<-",signature(x="spam",i="missing",j="vector", value = "ANY"),
	  function (x, i, j, value) {#cat(j);
            assign.spam(x,1:x@dimension[1],j,value)} )

setMethod("[<-",signature(x="spam",i="matrix",j="missing", value = "ANY"),
	  function (x, i, j, value) {#cat("Mm");
            assign.spam(x,i,NULL,value) })

setMethod("[<-",signature(x="spam",i="matrix",j="matrix",value = "ANY"),
	  function (x, i, j, value) {#cat("MM");
            assign.spam(x,cbind(c(i),c(j)),NULL,value) })

setMethod("[<-",signature(x="spam",i="spam",j="missing", value = "ANY"),
	  function (x, i, j, value)
{
    dimx <- x@dimension
    nrow <- dimx[1]
    ncol <- dimx[2]
    if ( i@dimension[1]>nrow | i@dimension[2]>ncol)
        stop("subscript out of bounds",call.=FALSE)
    if ( ( (i@rowpointers[i@dimension[1]+1]-1) %%length(value))!= 0)
        stop("number of items to replace is not a multiple of replacement length")

    nzmax <- as.integer(min(prod(nrow,ncol), i@rowpointers[i@dimension[1]+1]+x@rowpointers[nrow+1]-2))

    if (length(value)!=  (i@rowpointers[i@dimension[1]+1]-1) )
        value <- rep(value, (i@rowpointers[i@dimension[1]+1]-1) %/%length(value))

    ## z <- .Fortran("subass",
    ##               as.integer(nrow), as.integer(ncol),
    ##               as.double(x@entries),      as.integer(x@colindices),    as.integer(x@rowpointers),
    ##               b=as.double(value),  bj=as.integer(i@colindices), bi=as.integer(i@rowpointers),
    ##               c=vector("double",nzmax),jc=vector("integer",nzmax),ic=vector("integer",nrow+1),
    ##               nzmax=as.integer(nzmax),
    ##               PACKAGE="spam")
    if(getOption("spam.force64") || .format.spam(x)$package != "spam" || .format.spam(i)$package != "spam")
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("subass",
              ## subroutine subass(nrow,ncol,a,ja,ia,b,jb,ib,c,jc,ic,nzmax)
              SIGNATURE = c(SS$signature, SS$signature,
                            "double", SS$signature, SS$signature,
                            "double", SS$signature, SS$signature,
                            "double", SS$signature, SS$signature,
                            SS$signature ),

              nrow,
              ncol,

              x@entries,
              x@colindices,
              x@rowpointers,

              b = value,
              bj = i@colindices,
              bi = i@rowpointers,

              c = vector_dc( "double", nzmax),
              jc = vector_dc( SS$type, nzmax),
              ic = vector_dc( SS$type, nrow+1),

              nzmax = nzmax,

              INTENT = c("r", "r",
                         "r", "r", "r",
                         "r", "r", "r",
                         "w", "w", "w",
                         "r" ),
              PACKAGE = SS$package)

    cnz <- z$ic[nrow+1]-1

    return(.newSpam(
        entries = z$c[1:cnz],
        colindices = z$jc[1:cnz],
        rowpointers = z$ic[1:(nrow+1)],
        dimension = c(nrow,ncol)))
}      )

setMethod("[<-", signature(x = "spam", i = "ANY", j = "ANY", value = "ANY"),
	  function(x,i,j, value){#  cat(value,class(value))
          stop("invalid or not-yet-implemented 'spam' subassigning")})



"assign.spam" <-
function (x, rw, cl,value)
{
  # we separate into cases where:
  # (A) rw matrix:
  #     1: logical: transformation to spam and extract structure
  #     2: two column matrix: extract (i,j) as given by the lines.
  #     3: all else extract   x[ c( rw)]
  # (B) rw and cl one element: ((i,j)
  # (C) rw and cl vectors:  (i1:i2,j1:j2)               [i1<=i2, j1<=j2]
  #                         (c(i1,...,ii),c(j1,...,jj)) [arbitrary block]

  if (!is.numeric(value)) stop(paste("Assignment of numeric structures only, here",class(value)))

  dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]

  if (is.matrix(rw)) {
    if (is.logical(rw)) {
      return( x[as.spam(rw)] <- value)
    }
    if (dim(rw)[2]==2) {
      ir <- rw[,1]
      jr <- rw[,2]
    } else  {
      ir <- c(rw-1) %% nrow + 1
      jr <- c(rw-1) %/% nrow + 1
      rw <- cbind(ir,jr)
    }
    if ( (min(ir)<1)|(max(ir)>x@dimension[1])|(min(jr)<1)|(max(jr)>x@dimension[2]))
      stop("subscript out of bounds",call.=FALSE)
    if (any(duplicated(cbind(ir,jr))))
      stop("only unique index for subassigning",call.=FALSE)
    nir <- length(ir)
    if ( (nir%%length(value))!= 0)
      stop("number of items to replace is not a multiple of replacement length")
    value <- rep(value, nir%/%length(value))


    ord <- order(ir,jr)
    rw <- rw[ord,,drop=F]
    ## print("1")
    ## bia <- .Fortran("constructia",
    ##                 as.integer(nrow),as.integer(nir),
    ##                 rowpointers=vector("integer",nrow+1),
    ##                 ir=as.integer(c(rw[,1],0)),
    ##                 PACKAGE="spam")$rowpointers
    if( getOption("spam.force64") || nrow+1 > 2147483647 )
        SS <- .format64()
    else
        SS <- .format32

    bia <- .C64("constructia",
                SIGNATURE = c( SS$signature, SS$signature, SS$signature, SS$signature),

                nrow,
                nir,
                rowpointers = vector_dc( SS$type, nrow+1),
                ir = c(rw[,1],0),

                INTENT = c("r", "r", "w", "r"),
                PACKAGE = SS$package)$rowpointers
    nzmax <- as.integer(min(prod(nrow,ncol), nir+x@rowpointers[nrow+1])+2)
    ## z <- .Fortran("subass",
    ##               as.integer(nrow),as.integer(ncol),
    ##               as.double(x@entries), as.integer(x@colindices), as.integer(x@rowpointers),
    ##               b=as.vector(value[ord],"double"),
    ##               bj=as.vector(rw[,2],"integer"),  bi=as.integer(bia),
    ##               entries=vector("double",nzmax),
    ##               colindices=vector("integer",nzmax),
    ##               rowpointers=vector("integer",nrow+1),
    ##               nzmax=as.integer(nzmax),
    ##               NAOK=getOption("spam.NAOK"),PACKAGE="spam")
    if( getOption("spam.force64") ||  .format.spam(x)$package != "spam" || length( bia) >  2147483647 )
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("subass",
                  SIGNATURE = c(SS$signature, SS$signature,
                      "double", SS$signature, SS$signature,
                      "double", SS$signature, SS$signature,
                      "double", SS$signature, SS$signature,
                      SS$signature ),

                  nrow,
                  ncol,

                  x@entries,
                  x@colindices,
                  x@rowpointers,

                  b = as.vector( value[ord], "double" ),
                  bj = as.vector( rw[,2], "integer" ),
                  bi = bia,

                  entries = vector_dc( "double", nzmax),
                  colindices = vector_dc( SS$type, nzmax),
                  rowpointers = vector_dc( SS$type, nrow+1),

                  nzmax = nzmax,

                  INTENT = c("r", "r",
                      "r", "r", "r",
                      "r", "r", "r",
                      "w", "w", "w",
                      "r" ),
                  NAOK=getOption("spam.NAOK"),
                  PACKAGE = SS$package)
    cnz <- z$rowpointers[nrow+1]-1
    if (cnz<0) {
      cat("Negative cnz in subassigning, forced to one. Please report.")
      return( spam(0))
    }

    return(.newSpam(
        entries =  z$entries[1:cnz],
        colindices = z$colindices[1:cnz],
        rowpointers = z$rowpointers,
        dimension = c(nrow,ncol)
        ))

  }
  # negative subsetting:
  if ( max(rw)<0 )    rw <- seq_len( nrow)[rw]
  if ( max(cl)<0 )    cl <- seq_len( ncol)[cl]

  # logical
  if (is.logical(rw))    rw <- seq_len( nrow)[rw]
  if (is.logical(cl))    cl <- seq_len( ncol)[cl]

  # sanity check
  if (length(rw)==0) stop("You should assign at least one element for the rows",call.=FALSE)
  if (length(cl)==0) stop("You should assign at least one element for the columns",call.=FALSE)


  if ( (min(rw)<1)|(max(rw)>x@dimension[1])|(min(cl)<1)|(max(cl)>x@dimension[2]))
    stop("subscript out of bounds",call.=FALSE)

  if (is.vector(rw) && is.vector(cl)) {
    if (any(duplicated(rw))||any(duplicated(cl)))
      stop("only unique index for subassigning",call.=FALSE)

    nrw <- length(rw)   # length returns an integer, so is a product therof
    ncl <- length(cl)
    bnz <- nrw*ncl

    if ( (bnz%%length(value))!= 0)
      stop("number of items to replace is not a multiple of replacement length")

    # we pack the value into a vector _row by row_
    value <- c(t(array(as.double(value),c(nrw,ncl))[order(rw),order(cl)]))

    bia <- vector("integer",nrow)  # bia has size of nrow + 1
    bia[rw] <- ncl        # in each row we have ncl new objects
    bia <- as.integer(c(1,cumsum(bia)+1))

    # we construct now a sparse matrix containing the "value" at positions rw and cl.
    # then we use the subassign function.
    nzmax <- as.integer(min(prod(nrow,ncol), bnz+x@rowpointers[nrow+1])+2)

    ## print("2")
    ## z <- .Fortran("subass",
    ##               as.integer(nrow),as.integer(ncol),
    ##               as.double(x@entries), as.integer(x@colindices) ,as.integer(x@rowpointers),
    ##               b=as.double(value),
    ##               bj=as.integer(rep(sort(as.integer(cl)),nrw)),
    ##               bi=as.integer(bia),
    ##               entries=vector("double",nzmax),colindices=vector("integer",nzmax),
    ##               rowpointers=vector("integer",nrow+1),
    ##               nzmax=as.integer(nzmax),
    ##               NAOK=getOption("spam.NAOK"),PACKAGE="spam")
     if( getOption("spam.force64") ||  .format.spam(x)$package != "spam" || length( bia) >  2147483647 )
        SS <- .format64()
    else
        SS <- .format32
    z <- .C64("subass",
              SIGNATURE = c(SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  SS$signature ),

              nrow,
              ncol,

              x@entries,
              x@colindices,
              x@rowpointers,

              b = value,
              bj = rep_len64(sort(as.integer(cl)),nrw*length(cl)),
              bi = bia,

              entries = vector_dc( "double", nzmax),
              colindices = vector_dc( SS$type, nzmax),
              rowpointers = vector_dc( SS$type, nrow+1),

              nzmax = nzmax,

              INTENT = c("r", "r",
                  "r", "r", "r",
                  "r", "r", "r",
                  "w", "w", "w",
                  "r" ),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)
    cnz <- z$rowpointers[nrow+1]-1

    return(.newSpam(
        entries = z$entries[1:cnz],
        colindices = z$colindices[1:cnz],
        rowpointers = z$rowpointers,
        dimension = c(nrow,ncol)
        ))
  }
  stop("invalid or not-yet-implemented 'spam' subsetting")
}


#!8#
# x:spam, y:matrix, returns matrix
.spam.matmul.mat <- function(x,y)
{
    nrow <- x@dimension[1]
    ncol <- x@dimension[2]
    yrow <- dim(y)[1]
    ycol <- dim(y)[2]

    # EXPLICIT-STORAGE-FORMAT
    SS <- .format.spam(x)
    if( getOption("spam.force64")|| as.numeric(nrow)*as.numeric(ycol) > 2147483647)
        SS <- .format64()

    if(yrow != ncol) stop("not conformable for multiplication")
    z <- .C64("amuxmat",
              SIGNATURE=c(SS$signature, SS$signature, SS$signature, "double",
                  "double", "double", SS$signature, SS$signature),

              nrow,
              yrow,
              ycol,
              y,

              y = vector_dc("double",nrow*ycol),
              x@entries,
              x@colindices,
              x@rowpointers,

              INTENT=c("r", "r", "r", "r",
                       "w", "r", "r", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)$y
    dim(z) <- c(nrow,ycol)
    return(z)
}

.spam.matmul <- function(x,y) {
    # --- CHANGED ---
    # Refactored. -> See .spam.matmul.vector

  if (is.matrix(y)) y <- as.spam(y)
  if (is.matrix(x)) x <- as.spam(x)

  #matrix multiply two sparse spam matrices
  xn <- x@dimension[1]
  xm <- x@dimension[2]
  yl <- y@dimension[2]
  if(xm != y@dimension[1])
    stop("matrices not conformable for multiplication")

  # EXPLICIT-STORAGE-FORMAT
  SS <- .format.spam(x, y)
  if( getOption("spam.force64") || as.numeric(xn)*as.numeric(yl) > 2147483647)
        SS <- .format64()

  z <- .C64("amubdg",
            SIGNATURE = rep(SS$signature, 10),

            xn,
            xm,
            yl,
            x@colindices,

            x@rowpointers,
            y@colindices,
            y@rowpointers,
            unused2=vector_dc(SS$type, xn),

            nz = vector_dc(SS$type,1),
            unused=vector_dc(SS$type, yl),

            INTENT=c("r", "r", "r", "r",
                     "r", "r", "r", "w",
                     "w", "w"),
            NAOK = getOption("spam.NAOK"),
            PACKAGE = SS$package)
  nzmax <- z$nz
  z <- NULL

  # EXPLICIT-STORAGE-FORMAT
  # As we now know the number of entries, we can decide again if we need 64-bit
  SS <- .format.spam(x, y)
  if( getOption("spam.force64") || nzmax > 2147483647)
        SS <- .format64()

  z <- .C64("amub",
            SIGNATURE=c(SS$signature, SS$signature, SS$signature, "double",
                        SS$signature, SS$signature, "double", SS$signature,
                        SS$signature, "double", SS$signature, SS$signature,
                        SS$signature, SS$signature, SS$signature),

            xn,
            yl,
            1L,
            x@entries,

            x@colindices,
            x@rowpointers,
            y@entries,
            y@colindices,

            y@rowpointers,
            entries = vector_dc("double",nzmax),
            colindices = vector_dc(SS$type,nzmax),
            rowpointers = vector_dc(SS$type,xn+1),

            nzmax,
                                        # TODO: Check if we can get rid of this argument and allocate memory in Fortran:
            vector_dc(SS$type, yl),
            ierr = vector_dc(SS$type, 1),

            INTENT=c("r", "r", "r", "r",
                     "r", "r", "r", "r",
                     "r", "w", "w", "w",
                     "r", "r", "w"),
            NAOK = getOption("spam.NAOK"),
            PACKAGE = SS$package)
  nz <- z$rowpointers[xn+1]-1
  if(z$ierr != 0) stop("insufficient space for sparse matrix multiplication")


  if(nz==0){#trap zero matrix
    return(.newSpam(
      dimension=c(xn,yl)
    ))
  }

   entries <- z$entries
   colindices <- z$colindices
   rowpointers <- z$rowpointers
   z <- NULL

   length(entries) <- nz
   length(colindices) <- nz

   z <- .C64("sortrows",
             SIGNATURE=c(SS$signature, "double", SS$signature, SS$signature),

             xn,
             entries=entries,
             colindices=colindices,
             rowpointers=rowpointers,

             INTENT=c("r", "rw", "rw", "rw"),
             NAOK = getOption("spam.NAOK"),
             PACKAGE = SS$package)

  return(.newSpam(
    entries=z$entries,
    colindices=z$colindices,
    rowpointers=z$rowpointers,
    dimension=c(xn,yl)
  ))
}


.spam.matmul.vector <- function(x,y)
{
    # --- CHANGED ---
    # Refactored the function .spam.matmul into the functions
    # .spam.matmul and .spam.matmul.vector

    # if we have x*Y, we calculate t(t(Y)*x)
    # Use "is.spam(y)" instead of "is.vector(x)" because the later might be
    # misleading in case that x has any attributes attached.
    if(is.spam(y)) {
        A <- t(y)
        b <- x
    } else {
        A <- x
        b <- y
    }
    SS <- .format.spam(A)

    nrow <- A@dimension[1]
    ncol <- A@dimension[2]
    if(length(b) != ncol)  stop("not conformable for multiplication")
    z <- .C64("amux",
              NAOK = getOption("spam.NAOK"),
              SIGNATURE=c(SS$signature, "double", "double", "double",
                          SS$signature, SS$signature),

          nrow,
          b,
          y = vector_dc("double",nrow),
          A@entries,

          A@colindices,
          A@rowpointers,

          INTENT=c("r", "r", "w", "r", "r", "r"),
          PACKAGE = SS$package)$y

    if(is.spam(y))
        dim(z) <- c(1,nrow)
    else
        dim(z) <- c(nrow,1)
    return(z)
}

setMethod("%*%",signature(x="spam",y="spam"),    .spam.matmul)
setMethod("%*%",signature(x="spam",y="matrix"),  .spam.matmul.mat)
setMethod("%*%",signature(x="spam",y="numeric"), .spam.matmul.vector)
setMethod("%*%",signature(x="matrix",y="spam"),  .spam.matmul)
setMethod("%*%",signature(x="numeric",y="spam"), .spam.matmul.vector)

upper.tri.spam <- function(x,diag=FALSE)
  {
    dimx <- x@dimension
    nrow <- dimx[1]

    SS <- .format.spam(x)
    z <- .C64("getu",
              SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature,
                  "double", SS$signature, SS$signature),
              nrow,
              x@entries,
              x@colindices,
              x@rowpointers,

              entries = x@entries,
              colindices = x@colindices,
              rowpointers = x@rowpointers,

              INTENT = c("r", "r", "r", "r",
                  "w", "w", "w"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package )

    nz <- z$rowpointers[dimx[1]+1]-1

    if (!diag) {
        z <- .C64("getdia",
                  SIGNATURE = c(SS$signature, SS$signature, SS$signature,
                                    "double", SS$signature, SS$signature,
                      SS$signature, "double", SS$signature, SS$signature),

                  n = nrow,
                  m = nrow,
                  job = 1,

                  entries = z$entries[1:nz],
                  colindices = z$colindices[1:nz],
                  rowpointers = z$rowpointers,

                  len = nrow,
                  diag = vector_dc("double", nrow),
                  idiag = vector_dc(SS$type, nrow),
                  ioff = 0,

                  INTENT = c("r", "r", "r",
                      "rw", "rw", "rw",
                      "w", "w", "w", "r"),
                  NAOK = getOption("spam.NAOK"),
                  PACKAGE = SS$package
                  )
        nz <- z$rowpointers[nrow+1]-1
    }
    if(getOption("spam.trivalues")){
        return(.newSpam(
            entries = z$entries[1:nz],
            colindices = z$colindices[1:nz],
            rowpointers = z$rowpointers,
            dimension = dimx,
            ))
    }else{
        return(.newSpam(
            entries = rep_len64(1,nz),
            colindices=z$colindices[1:nz],
            rowpointers=z$rowpointers,
            dimension=dimx
            ))
    }
  }

lower.tri.spam <- function(x,diag=FALSE)
{
    dimx <- x@dimension
    nrow <- dimx[1]
    SS <- .format.spam(x)
    z <- .C64("getl",
              SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature,
                  "double", SS$signature, SS$signature),
              nrow,
              x@entries,
              x@colindices,
              x@rowpointers,

              entries = x@entries,
              colindices = x@colindices,
              rowpointers = x@rowpointers,

              INTENT = c("r", "r", "r", "r",
                  "w", "w", "w"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package )
    nz <- z$rowpointers[nrow+1]-1

    if (!diag) {
        z <- .C64("getdia",
                  SIGNATURE = c(SS$signature, SS$signature, SS$signature,
                      "double", SS$signature, SS$signature,
                      SS$signature, "double", SS$signature, SS$signature),

                  n = nrow,
                  m = nrow,
                  job = 1,

                  entries = z$entries[1:nz],
                  colindices = z$colindices[1:nz],
                  rowpointers = z$rowpointers,

                  len = nrow,
                  diag = vector_dc("double", nrow),
                  idiag = vector_dc(SS$type, nrow),
                  ioff = 0,

                  INTENT = c("r", "r", "r",
                      "rw", "rw", "rw",
                      "w", "w", "w", "r"),
                  NAOK = getOption("spam.NAOK"),
                  PACKAGE = SS$package
                  )
        nz <- z$rowpointers[nrow+1]-1
    }

    if(getOption("spam.trivalues")){
        return(.newSpam(
            entries = z$entries[1:nz],
            colindices = z$colindices[1:nz],
            rowpointers = z$rowpointers,
            dimension = dimx,
            ))
    }else{
        return(.newSpam(
            entries = rep_len64(1,nz),
            colindices=z$colindices[1:nz],
            rowpointers=z$rowpointers,
            dimension=dimx
            ))
    }
}



setGeneric("upper.tri")
setMethod("upper.tri","spam",upper.tri.spam)
setGeneric("lower.tri")
setMethod("lower.tri","spam",lower.tri.spam)


# fields uses the construct of vector representation for a diagonal matrix.

# Create a special matrix multiply for diagonal matrices.
# Diagonal matrix assumed to be just a vector.
# NOTE: this is not a symmetric operation:
#  when a left vector is given it is a diagonal matrix
#  when a right vector is given it is a vector.
#
.spam.diagmulmat <- function(x,y){
    nrow <- y@dimension[1]
    if(length(x) != nrow)
        stop("not conformable for multiplication")

    SS <- .format.spam(y)

    z <- .C64("diagmua",
              SIGNATURE = c(SS$signature,
                  "double", SS$signature, "double"),

              nrow,

              entries = y@entries,
              y@rowpointers,
              as.vector(x,"double"),

              INTENT = c("r",
                  "rw", "r", "r"),
              NAOK=getOption("spam.NAOK"),
              PACKAGE = SS$package)$entries

    y@entries <- z
    return(y)
}

.spam.diagaddmat <- function(x,y){
      ## subroutine diagaddmat (nrow, n, a, ja, ia, diag, iw)
    nrow <- y@dimension[1]
    minrc <- min( y@dimension)
    if(length(x) != minrc) stop("not conformable for addition")

    SS <- .format.spam(y)

    z <- .C64("diagaddmat",
              SIGNATURE = c( SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  "double", SS$signature),

              nrow = nrow,
              n = minrc,

              a = c(y@entries, rep_len64(0,minrc)),
              ja = c(y@colindices, rep_len64(0,minrc)),
              ia = y@rowpointers,

              diag = x,
              iw = vector_dc(SS$type, nrow),

              INTENT = c("r", "r",
                  "rw", "rw", "rw",
                  "r", "rw"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package )

    nz <- z$ia[nrow+1]-1

    return(.newSpam(
        entries = z$a[1:nz],
        colindices = z$ja[1:nz],
        rowpointers = z$ia,
        dimension = y@dimension))
}

setGeneric("%d*%",function(x,y,...)standardGeneric("%d*%"))

setMethod("%d*%",signature(x="matrix",y="ANY"),       function(x,y){x%*%y} )
setMethod("%d*%",signature(x="numeric",y="matrix"),   function(x,y){x*y} )
setMethod("%d*%",signature(x="numeric",y="numeric"),  function(x,y){cbind(x*y)} )

setMethod("%d*%",signature(x="spam",y="spam"),    .spam.matmul )
setMethod("%d*%",signature(x="spam",y="ANY"),     .spam.matmul )
setMethod("%d*%",signature(x="numeric",y="spam"), .spam.diagmulmat )


setGeneric("%d+%",function(x,y,...)standardGeneric("%d+%"))
setMethod("%d+%",signature(x="matrix",y="ANY"),      function(x,y){ x+y } )
setMethod("%d+%",signature(x="numeric",y="matrix"),  function(x,y){ diag(x)+y} )
setMethod("%d+%",signature(x="numeric",y="numeric"), function(x,y){ diag(x)+y} )

setMethod("%d+%",signature(x="spam",y="spam"),     function(x,y){ x+y})
setMethod("%d+%",signature(x="spam",y="ANY"),      function(x,y){ x+y})
setMethod("%d+%",signature(x="numeric",y="spam"),  .spam.diagaddmat )


all.equal.spam <- function (target, current, tolerance = .Machine$double.eps^0.5,
    scale = NULL, check.attributes = FALSE,...)
{
    if (check.attributes)
        warning("attributes are not supported for 'spam' objects. Ignoring 'check.attributes' argument")
    if (!is.spam(target)) stop("'target' should be of class 'spam'")
    if (!is.spam(current)) {
        return(paste("target is spam, current is ", data.class(current), sep = ""))
    }
    msg <- NULL
    lt <- length(target)
    lc <- length(current)
    if (lt != lc) {
      return(paste("Lengths (", lt, ", ", lc, ") differ", sep = ""))
    }
    dt <- target@dimension
    dc <- current@dimension
    if ( !all( dt == dc ))
      return(paste("Dimensions ([",dt[1],",",dt[2],"], [",
                    dc[1],",",dc[2], "]) differ", sep = ""))
    # --- CHANGED ---
    # Add suppressWarnings
    tmp <- suppressWarnings(sum(target@colindices != current@colindices))
    if ( tmp > 0)
      msg <- c(msg,paste("Column-sparsity structure differ (at least",
                    tmp,"instance(s))"))

    tmp <- suppressWarnings(sum(target@rowpointers != current@rowpointers))
    if ( tmp > 0)
      msg <- c(msg,paste("Row-sparsity structure differ (at least",
                    tmp,"instance(s))"))

    what <- if (is.null(scale)) {
      scale <- mean(abs(target@entries))
      if (is.finite(scale) && scale > tolerance) {
        "relative"
      } else {
        scale <- 1
        "absolute"
      }
    } else {
        "scaled"
    }
    xy <- suppressWarnings(mean(abs(target@entries - current@entries) / scale))
    if (is.na(xy) || xy > tolerance)
        msg <- c(msg,paste("Mean", what, "difference:", format(xy)))
    if (is.null(msg))
        TRUE
    else msg
}


isSymmetric.spam <- function(object, tol = 100 * .Machine$double.eps, ...)
{
  # very similar to is.Symmetric.matrix
  test <-  all.equal.spam(object, t.spam(object), tolerance = tol, ...)

  # Possibility that structure is different but not contents

  if (!isTRUE(test)) {
    object <- as.spam.spam(object)
    test <-  all.equal.spam(object, t.spam(object), tolerance = tol, ...)
  }
  isTRUE(test)

}


setMethod("all.equal",signature(target="spam",current="spam"), all.equal.spam )
setMethod("all.equal",signature(target="matrix",current="spam"),
          function (target, current, tolerance = .Machine$double.eps^0.5,
                    scale = NULL, check.attributes = FALSE,eps = getOption("spam.eps"),...)
{
    if (check.attributes)
      warning("attributes are not supported for 'spam' objects. Ignoring 'check.attributes' argument")
    msg <- NULL

    dimx <- dim(target)
    nz <- length(target)
     ## z <- .Fortran("spamdnscsr", nrow = as.integer(dimx[1]), ncol = as.integer(dimx[2]),
     ##    x = as.double(target), as.integer(dimx[1]), entries = vector("double",
     ##        nz), colindices = vector("integer", nz), rowpointers = vector("integer",
     ##        dimx[1] + 1), eps = as.double(eps),
     ##              NAOK = getOption("spam.NAOK"), PACKAGE = "spam")
    if( getOption("spam.force64") || .format.spam(current)$package != "spam" || prod(dimx) > 2147483647)
        SS <- .format64()
    else
        SS <- .format32
    z <- .C64("spamdnscsr",
              SIGNATURE = c(SS$signature, SS$signature, "double", SS$signature,
                  "double", SS$signature, SS$signature, "double"),

              nrow = dimx[1],
              ncol = dimx[2],
              x = target,
              dimx[1],
              entries = vector_dc( "double", nz),
              colindices = vector_dc( SS$type, nz),
              rowpointers = vector_dc( SS$type, dimx[1] + 1),
              eps = eps,

              INTENT = c("r", "r", "r", "r",
                  "w", "w", "w", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)

    lt <- z$rowpointers[dimx[1] + 1] - 1
    lc <- length(current)

    if (lt != lc) {
        return(paste("Lengths (", lt, ", ", lc, ") differ", sep = ""))
    }
    dt <- dim(target)
    dc <- current@dimension
    if ( !all( dt == dc ))
        return(paste("Dimensions ([",dt[1],",",dt[2],"], [",
                     dc[1],",",dc[2], "]) differ", sep = ""))
    tmp <- sum(z$colindices[1:lt] != current@colindices)
    if ( tmp>0)
        msg <- c(msg,paste("Column-sparsity structure differ (at least",
                           tmp,"instance(s))"))

    tmp <- sum(z$rowpointers != current@rowpointers)
    if ( tmp>0)
        msg <- c(msg,paste("Row-sparsity structure differ (at least",
                           tmp,"instance(s))"))

    xy <- mean(abs(z$entries[1:lt] - current@entries))
    what <- if (is.null(scale)) {
        xn <- mean(abs(z$entries))
        if (is.finite(xn) && xn > tolerance) {
            xy <- xy/xn
            "relative"
        }
        else "absolute"
    }
    else {
        xy <- xy/scale
        "scaled"
    }
    if (is.na(xy) || xy > tolerance)
        msg <- c(msg,paste("Mean", what, "difference:",
            format(xy)))
    if (is.null(msg))
        TRUE
    else msg

}
 )
setMethod("all.equal",signature(target="spam",current="matrix"), function(target, current, ...)
     { all.equal(current, target, ...) })
setMethod("isSymmetric","spam", isSymmetric.spam)
