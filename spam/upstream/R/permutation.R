# HEADER ####################################################
# This is file spam/R/permutation.R.                        #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################








checkpivot <- function(pivot, len, type="Pivot") {
  if(is.null(pivot))                       return()
  if(!is.vector(pivot))                        stop(paste(type,"is not a vector."))     
  pivot <- as.vector(pivot,"integer")
  if (length(pivot) !=len)           stop(paste(type,"of wrong length."))
  tmp <- sort.int(pivot)
  if(tmp[1]!=1 ||  any(tmp-seq_len(len)!=0))   stop(paste("Invalid",type))
  return()
}



"permutation.spam" <- function(A, P=NULL, Q=NULL, ind=FALSE, check=TRUE){
  # eliminated .Internal calls as this creates a "Note" on CRAN checks.
                                        # Only 1-2% timing loss, see end of the file.
    if(  getOption("spam.force64") )
        SS <- .format64()
    else
        SS <- .format.spam(A)
  
    
    nrow <- A@dimension[1]
    ncol <- A@dimension[2]
    
    if (is.null(P)&is.null(Q))
        stop("At least one permutation should be specified")

    nz <- A@rowpointers[nrow+1]-1


    if (check){
        checkpivot(P,nrow,"Permutation")
        checkpivot(Q,ncol,"Permutation")
    }
    
    if (is.null(Q)) {
                                        #      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm)
                                        #      B = P A
        P <- as.integer(P)
        if(ind) P <- order(P)
                                        #    if(ind) P <- .Internal(order(T,F,P))
        
        
        z <- .C64("rperm",
                  ## subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm)
                  SIGNATURE = c( SS$signature, "double", SS$signature, SS$signature,
                                "double", SS$signature, SS$signature,
                                SS$signature),
                  
                  nrow,
                  A@entries,
                  A@colindices,
                  A@rowpointers,
                  
                  entries = vector_dc( "double", nz), 
                  colindices = vector_dc( SS$type, nz),
                  rowpointers = vector_dc( SS$type, nrow + 1),
                  
                  P,

                  INTENT = c("r", "r", "r", "r",
                             "w", "w", "w",
                             "r"),
                  NAOK = getOption("spam.NAOK"),
                  PACKAGE = SS$package)
    } else {  
        if (is.null(P)){
                                        #      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,iwork) 
                                        #      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), iwork(*)
                                        #      double precision a(*), ao(*) 
                                        #      B = A Q 
            Q <- as.integer(Q)
            if(ind) Q <- order(Q)
                                        #       if(ind) Q <- .Internal(order(T,F,Q))
            ## z <- .Fortran("cperm",
            ##               as.integer(nrow),
            ##               A@entries, as.integer(A@colindices),
            ##               as.integer(A@rowpointers),
            ##               entries = vector("double",nz),
            ##               colindices = vector("integer", nz),
            ##               rowpointers = vector("integer", nrow + 1),
            ##               Q,
            ##               NAOK = getOption("spam.NAOK"), PACKAGE = "spam")
            z <- .C64("cperm",
                      SIGNATURE = c(SS$signature, "double", SS$signature, SS$signature,
                                    "double", SS$signature, SS$signature,
                                    SS$signature),
                      
                      nrow,
                      A@entries,
                      A@colindices,
                      A@rowpointers,
                      
                      entries = vector_dc( "double", nz),
                      colindices = vector_dc( SS$type, nz),
                      rowpointers = vector_dc( SS$type, nrow + 1),
                      
                      Q,

                      INTENT = c("r", "r", "r", "r",
                                 "w", "w", "w",
                                 "r"),
                      NAOK = getOption("spam.NAOK"),
                      PACKAGE = SS$package)

        } else {  
                                        #      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,pperm,qperm,iwork)
                                        #      B = P A Q 
            Q <- as.integer(Q)
                                        #        if(ind) Q <- .Internal(order(T,F,Q))
            if(ind) Q <- order(Q)
            P <- as.integer(P)
                                        #        if(ind) P <- .Internal(order(T,F,P))
            if(ind) P <- order(P)
            ## z <- .Fortran("dperm",
            ##               as.integer(nrow),
            ##               A@entries,as.integer(A@colindices),
            ##               as.integer(A@rowpointers),
            ##               entries = vector("double",nz),
            ##               colindices = vector("integer", nz),
            ##               rowpointers = vector("integer", nrow + 1),
            ##               P,Q,
            ##               NAOK = getOption("spam.NAOK"), PACKAGE = "spam")
            z <- .C64("dperm",
                      ## subroutine dperm (nrow,a,ja,ia,ao,jao,iao,pperm,qperm)
                      SIGNATURE = c( SS$signature, "double", SS$signature, SS$signature,
                                    "double", SS$signature, SS$signature,
                                    SS$signature, SS$signature),
                      
                      nrow,
                      A@entries,
                      A@colindices,
                      A@rowpointers,
                      
                      entries = vector_dc( "double", nz),
                      colindices = vector_dc( SS$type, nz),
                      rowpointers = vector_dc( SS$type, nrow + 1),
                      
                      P,Q,

                      INTENT = c("r", "r", "r", "r",
                                 "w", "w", "w",
                                 "r", "r"),
                      NAOK = getOption("spam.NAOK"),
                      PACKAGE = SS$package)
        }   
        
    }
    ## newx <- new("spam")
    ## slot(newx, "entries", check = FALSE) <- z$entries
    ## slot(newx, "colindices", check = FALSE) <- z$colindices
    ## slot(newx, "rowpointers", check = FALSE) <- z$rowpointers
    ## slot(newx, "dimension", check = FALSE) <- c(nrow,ncol)
    return(.newSpam(
        entries = z$entries,
        colindices = z$colindices,
        rowpointers = z$rowpointers,
        dimension = c(nrow,ncol)
    ))
}

permutation.matrix <- function(A, P=NULL, Q=NULL, ind=FALSE, check=TRUE){
  nrow <- dim(A)[1]
  ncol <- dim(A)[1]
  
  if (is.null(P)&is.null(Q))     stop("At least one permutation should be specified")
  
  if (check){
    checkpivot(P,nrow,"Permutation")
    checkpivot(Q,ncol,"Permutation")
  }

  if (ind) {
    if (is.null(Q))     return(A[P,])
    if (is.null(P))     return(A[,Q])
    return(A[P,Q])
  } else {
    if (is.null(Q))     return(A[order(P),])
    if (is.null(P))     return(A[,order(Q)])
    return(A[order(P),order(Q)])
  }
}

setGeneric("permutation",function(A, P=NULL, Q=NULL, ind=FALSE, check=TRUE)standardGeneric("permutation"))
setMethod("permutation","matrix",permutation.matrix)
setMethod("permutation","spam",permutation.spam)


### ss <- sample(1:100000)
### system.time( for( i in 1:1000) tt<-order(ss))
### system.time( for( i in 1:1000) tt<-.Internal(order(T,F,ss)))
