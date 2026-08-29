# HEADER ####################################################
# This is file spam/R/kronecker.R.                          #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################






kronecker.spam <- function(X,Y,FUN = "*", make.dimnames = FALSE, ...)
{
    ## print("1")
    if(getOption("spam.force64"))
        SS <- .format64()
    else
        SS <- .format32

  if (make.dimnames)
    warning("dimnames not supported within sparse matrices")
  lenx <- length(X)
  leny <- length(Y)
  kronlen <- lenx*leny

  if(is.spam(X)){
    Xdim <- X@dimension
    Xentries <- X@entries
    Xcol <- X@colindices
    Xrow <- X@rowpointers
    if(.format.spam(X)$package == "spam64")
        SS <- .format64()
  }else if (is.vector(X)){
    Xentries <- X
    Xcol <- rep.int(as.integer(1),lenx)
    Xdim <- as.integer(c(lenx,1))
    Xrow <- seq_len(lenx+1)
  } else {
    Xentries <- as.double(t(X))
    Xdim <- dim(X)
    Xcol <- rep.int(as.integer(1:Xdim[2]),Xdim[1])
    Xrow <- seq.int(1, by=Xdim[2], length.out=Xdim[1]+1)
  }

  if(is.spam(Y)){
    Ydim <- Y@dimension
    Yentries <- Y@entries
    Ycol <- Y@colindices
    Yrow <- Y@rowpointers
    if(.format.spam(Y)$package == "spam64")
        SS <- .format64()
  } else if (is.vector(Y)){
    Yentries <- Y
    Ycol <- rep.int(as.integer(1),leny)
    Ydim <- as.integer(c(leny,1))
    Yrow <- seq_len(leny+1)
  } else {
    Yentries <- as.double(t(Y))
    Ydim <- dim(Y)
    Ycol <- rep.int(as.integer(1:Ydim[2]),Ydim[1])
    Yrow <- seq.int(1, by=Ydim[2], length.out=Ydim[1]+1)
  }
  ## kronxy <- new("spam")

    if (FUN=="*") {
      z <- .C64("kroneckermult",
                SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature,
                              SS$signature, SS$signature, "double", SS$signature, SS$signature,
                              "double", SS$signature, SS$signature),

                Xdim[1],
                Xentries,
                Xcol,
                Xrow,

                Ydim[1],
                Ydim[2],
                Yentries,
                Ycol,
                Yrow,

                entries = vector_dc( "double", kronlen),
                colindices = vector_dc( SS$type, kronlen),
                rowpointers = vector_dc( SS$type, Xdim[1]*Ydim[1]+1),

                INTENT = c( "r", "r", "r", "r",
                           "r", "r", "r", "r", "r",
                           "w", "w", "w"),
                NAOK=getOption("spam.NAOK"),
                PACKAGE = SS$package)
        kronxy <- .newSpam( entries = z$entries)
    ## slot(kronxy, "entries", check=FALSE) <-  z$entries
    }else {
        z <- .C64("kroneckerf",
                 SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature,
                               SS$signature, SS$signature, "double", SS$signature, SS$signature,
                               "double", "double", SS$signature, SS$signature),

                 Xdim[1],
                 Xentries,
                 Xcol,
                 Xrow,

                 Ydim[1],
                 Ydim[2],
                 Yentries,
                 Ycol,
                 Yrow,

                 ent1 = vector_dc( "double", kronlen),
                 ent2 = vector_dc( "double", kronlen),
                 colindices = vector_dc( SS$type, kronlen),
                 rowpointers = vector_dc( SS$type, Xdim[1]*Ydim[1]+1),

                 INTENT = c("r", "r", "r", "r",
                            "r", "r", "r", "r", "r",
                            "w", "w", "w", "w"),
                 NAOK=getOption("spam.NAOK"),
                 PACKAGE = SS$package)

        FUN <- match.fun(FUN)
        kronxy <- .newSpam( entries = FUN(z$ent1,z$ent2,...) )
    ## slot(kronxy, "entries", check=FALSE) <-  FUN(z$ent1,z$ent2,...)
    if (z$rowpointers[Xdim[1]*Ydim[1]+1]-1 < prod(Xdim,Ydim))
      warning("Sparseness structure of 'kronecker(X,Y)' preserved when applying 'FUN'.", call.=FALSE)
    }
    kronxy <- .newSpam(
        entries = kronxy@entries,
        colindices = z$colindices,
        rowpointers = z$rowpointers,
        dimension = Xdim*Ydim)
  ## slot(kronxy, "colindices", check=FALSE) <- z$colindices
  ## slot(kronxy, "rowpointers", check=FALSE) <- z$rowpointers
  ## slot(kronxy, "dimension", check=FALSE) <- Xdim*Ydim   # no need to use prod

  return(kronxy)
}


"kronecker.default" <- base::kronecker
setGeneric("kronecker")

setMethod("kronecker","spam", kronecker.spam)
setMethod("kronecker",signature(X="spam",Y="spam"), kronecker.spam)
setMethod("kronecker",signature(X="spam",Y="ANY"), kronecker.spam)
setMethod("kronecker",signature(X="ANY",Y="spam"), kronecker.spam)
