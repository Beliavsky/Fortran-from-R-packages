# HEADER ####################################################
# This is file spam/R/dim.R.                                #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################


# This is the actual dim...

"dim<-.spam" <- function(x, value) {
    if (is.spam(x)) {

        dimx <- x@dimension
        pdim <- prod(dimx)
        vlen <- prod(value)
        if( !identical(pdim,vlen))
            stop( sprintf("dims [product %d] do not match the length of object [%d]. Do you want `pad`",
                          pdim,vlen))

        if (length(value)>2)
            stop("dims should be of length 1 or 2")
        if (identical(length(value),1L))
            return( c(x) )

        if(any(dimx<1))
            stop("the dims contain negative values")

        tmp <- cbind(st=rep(1:dim(x)[1],diff(x@rowpointers)), nd=x@colindices)
        ind <- tmp[,1]+(tmp[,2]-1)*dimx[1] - 1

        slist <- list(i = ind%%value[1]   +1,
                       j = ind%/%value[1] +1,
                       x@entries)

        return( spam.list( slist, nrow=value[1], ncol=value[2],
                          eps = .Machine$double.eps))


    } else  {
        dim(x) <- value
        x
    }
}


########################################################################
# dim and derivatives

"pad<-.spam" <- function(x,value) {

    force64 <- getOption("spam.force64")

    # check if value is valid
    if ( (min(value)<1 ) || any(!is.finite(value)))
        stop("dims should be postive integers.")
    if (!identical( length(value), 2L))
        stop("dims should be of length 2.")

    dimx <- x@dimension
    last <- value[1]+1

  # In three steps:
  #  1) Address col truncation
            # to safe time, we also take into account if we have fewer or equal rows
  #  2) Augment rows
  #  3) if fewer rows and more columns, truncate
  # In any case, dimensions are fixed at the end.


  # If fewer cols required, we run reducedim
  if (dimx[2]>value[2]){

#     subroutine reducedim(a,ja,ia,eps,bnrow,bncol,k,b,jb,ib)

      if( force64 || .format.spam(x)$package == "spam64")
          SS <- .format64()
      else
          SS <- .format32

      z <- .C64("reducedim",
                SIGNATURE=c("double", SS$signature, SS$signature,
                    "double", SS$signature, SS$signature, SS$signature,
                    "double", SS$signature, SS$signature),
                oldra = x@entries,
                oldja = x@colindices,
                oldia = x@rowpointers,

                eps = getOption("spam.eps"),
                min(value[1],dimx[1]),
                value[2],
                nz = 1,

                entries=vector_dc("double",length(x@entries)),
                colindices=vector_dc(SS$type,length(x@entries)),
                rowpointers=vector_dc(SS$type,last),

                INTENT=c("r", "r", "r",
                         "r", "r", "r", "w",
                         "w", "w", "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE = SS$package)

    if (z$nz==1 ){ #was identical( z$nz,1L)
        ## print("2")
            return(
                .newSpam(
                    # entries=x@entries,
                    # colindices=x@colindices,
                    # rowpointers=c(1,rep_len64(2,value[1])),
                    dimension=value,
                    force64=force64
                    )
                )
    }
      nz <- z$nz-1
      x <- .newSpam(
          entries=z$entries[1:nz],
          colindices=z$colindices[1:nz],
          rowpointers=z$rowpointers[1:min(last,dimx[1]+1)],
          dimension=value, #actually here dim 2 = value 2 but dim1 maybe not yet
          force64=force64
          )
  }
    # augment rows
  if  (dimx[1]<value[1]){
      ## print("3")
      x <- .newSpam(
              entries=x@entries,
              colindices=x@colindices,
              rowpointers= c( x@rowpointers,
                  rep_len64( x@rowpointers[length(x@rowpointers)],value[1]-dimx[1])),
              dimension=value,
              force64=force64
              )
  }

    # special case: fewer rows and more columns, truncate
  if((dimx[1]>=value[1])&(dimx[2]<=value[2])) { ## added =, think about it again
      ## print("4")
      lastelement <- (x@rowpointers[last]-1)

      x <- .newSpam(
          entries= x@entries[1:lastelement],
          colindices= x@colindices[1:lastelement],
          rowpointers= x@rowpointers[1:last],
          dimension=value,
          force64=force64
          )
  }
    #before dim x = value x was here with slot option
  return(x)

}



setMethod("dim",   "spam", function(x) x@dimension )
setMethod("dim<-",   "spam", get("dim<-.spam"))

setGeneric("pad<-", function(x, value) standardGeneric("pad<-"))
setMethod("pad<-",   "spam", get("pad<-.spam"))
setMethod("pad<-",   "matrix",
          function(x, value) {
              if (!identical( length(value), 2L)) stop("dims should be of length 2.")
              tmp <- matrix(0, value)
              mr <- 1:min(value[1], dim(x)[1])
              mc <- 1:min(value[2], dim(x)[2])
              tmp[mr,mc] <- x[mr,mc]
              return(tmp)
          })
