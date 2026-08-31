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
# Contains common routine for chol.spam and determinant.spam 
# The lattre and Associated S4 elements are in different files.
#
#
########################################################################
########################################################################


# lindx=  colindices       of length nsub   = nnzR
# xlindx= colpointers                nsuper = nnzcolindices
# xlnz=   rowpointers
# snode=  snmember
# xsuper= supernodes
# c(... nnztmp,cachesize)= memory

#### Summary of memory usage in cholstepwise:
# Fortran      R
# m            nrow,ncol
# nsuper                     nnzcolindices
# nnzlmax                    nnzR
#
# Internally, we additionally assign several arrays (of length):
#  m: adj, colcnt
#  7*m+3: iwork
#  iwork: but 7*m+3 is used!

.chol.spam.basis <- function(x,  pivot,  memory,  eps,  verbose){

  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision", call.=FALSE)

  nrow <- x@dimension[1]
  nnzA <- x@rowpointers[nrow+1]-1L
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'", call.=FALSE)

  if (any( diag.of.spam(x, nrow, nrow) < eps))
    stop("Input matrix to 'chol' not positive definite (up to eps)", call.=FALSE)
# base rule:
#  `nrow(x) * .Machine$double.neg.eps * max(diag(x)`.

  if (length(pivot)==1) {
        if (pivot==FALSE) {
            doperm <- 0
         } else if(pivot==TRUE) {
            doperm <- 1
        } else {
            doperm <- as.integer( switch(match.arg(pivot,c("MMD","RCM")), MMD=1, RCM=2))
        }
  } else  if (length(pivot)==nrow) {
        doperm <- 0
        if (getOption("spam.cholpivotcheck")) {
            checkpivot(pivot,nrow)
        }
  } else stop("'pivot' should be 'MMD', 'RCM' or a valid permutation")

  ### IMPROVEME get better parameter values
  nnzcfact <- c(5,1,5)
  nnzRfact <- c(5,1,2)
  # nnzcolindices = length of array holding the colindices
  if(is.null(memory$nnzcolindices))  {
    nnzcolindices <- ifelse((nnzA/nrow < 5), # very sparse matrix
                           max(1000,nnzA*(1.05*nnzA/nrow-3.8)),
                           nnzA)*nnzcfact[doperm+1]
    nnzcolindices <- max(nnzcolindices,nnzA)
    if (verbose) cat("Guessing upper bound of column indicies as ", nnzcolindices, " (",
            printSize(nnzcolindices),").\n", sep='')
  } else {
    nnzcolindices <- max(memory$nnzcolindices,nnzA)
    memory$nnzcolindices <- NULL
  }
  # nnzR = length of array holding the nonzero values of the factor
  if(is.null(memory$nnzR)) {
    nnzR <- min(max(4*nnzA,floor(.2*nnzA^1.3))*nnzRfact[doperm+1],nrow*(nrow+1)/2)
    if (verbose) cat("Guessing upper bound of factor size.\n")
  } else {
    nnzR <- memory$nnzR
    memory$nnzR <- NULL
  }


  if(is.null(memory$cache))     cache <- 512  else {
    cache <- memory$cache
    memory$cache <- NULL
  }

  if (length( memory)>0 )
    warning("The component(s) ", paste("'",names(memory),"'",sep='',collapse=","),
            " of the argument 'memory'\npassed to function 'chol' not meaningful and hence ignored.",call.=FALSE)


  if(getOption("spam.cholsymmetrycheck")) {
    if (verbose) cat("Checking symmetry.")    
    test <- isSymmetric(x, tol = eps*100)
#  from help of isSymmetric:
#     isSymmetric(object, tol = 100 * .Machine$double.eps, ...)
    if (!isTRUE(test))
      stop("Input matrix to 'chol' not symmetric (up to 100*eps)",call.=FALSE)
  }


    if (verbose) {
        cat("Factorizing matrix of dimension ", nrow, " with ", length(x@entries)," entries (",
            printSize(x),  "). \nReserved memory for factor with ", nnzR, " entries (",
            printSize(nnzR),").\n", sep='')
        }

    force64 <- getOption("spam.force64")
    # we elaborate on the sizes and check for 2^31-2 to allow one additional iteration per counter
    if(force64 || 7*nrow+3  > 2147483646 || nnzcolindices  > 2147483646 || nnzR  > 2147483646) {
        SS <- .format64()
    } else {
        SS <- .format32
    }


    if (length(pivot)==1) {
        if (pivot==FALSE) {
            pivot <- as.vector( seq_len(nrow), SS$type)
        } else if(pivot==TRUE) {
            pivot <- vector(SS$type,nrow)
        } else {
           pivot <- vector(SS$type, nrow)
        }
    } else  pivot <- as.vector( pivot, SS$type)



    if (verbose) {
        cat("\nWorking with",SS$name,"spam matrix")
        if (SS$type != class(x@dimension)) cat(", casting original matrix")

        ## note that object.size(x)  is not required as we have reading mode
        m1 <- object.size(x) +
            6*object.size(pivot)+
            (nnzcolindices +     # length lindx
             2*nrow +nnzA + 7*nrow+3)*   # work arrays {colcnt, split}, {adj, adjncy}  iwork
            ifelse(SS$signature=='double',8,4) +      # depending on 32/64-bit
            nnzR*8    # final entries
        cat(".\nApproximate amount of memory needed: ", printSize( m1),".\n",
            "Starting now the heavy calculation (this might take some time ...\n",sep='')
                                        # gestimate for total use is:
#        m2 <- (8+6)*object.size(pivot)+
#            nnzcolindices*ifelse(SS$signature=='double',8,4) + nnzR*8
    }
  cholstep.intent <- c("rw", "rw", "r", "r",  # we are passing back a few of the arguments 
                    "r", "r", "rw", "rw",
                    "rw", "rw", "w", "w",    # 'w' for lindx xlindx  # _FOR SURE_
                    "rw", "rw", "w", "w",   # 'w' for lnz  xlnz
                    "rw", "rw", "rw", "rw")

  cholstep.signature <- rep(SS$signature, 20)
  cholstep.signature[c(3,15)] <- "double"
  z <- .C64("cholstepwise",                # 11 huge objects...
            SIGNATURE=cholstep.signature,

            nrow = as.vector(nrow, mode=SS$type),
            nnzA = as.vector(nnzA, SS$type),
            d =  x@entries,                          #Size _D_nrow
            jd = as.vector(x@colindices, SS$type),   #Size nrow

            id = as.vector(x@rowpointers, SS$type),#5#Size nrow
            doperm = as.vector(doperm,    SS$type),
            invp = vector(SS$type,nrow),             #Size nrow
            perm = pivot,                            #Size nrow

            nnzlindx = vector(SS$type, 1), #9
            nnzcolindices = as.vector(nnzcolindices, SS$type),
            lindx = vector_dc(SS$type, nnzcolindices), #Size nnzcolindices
            xlindx = vector_dc(SS$type, nrow+1),       #Size nrow

            nsuper = vector(SS$type, 1), #13
            nnzR = as.vector(nnzR, SS$type),
            lnz = vector_dc( "double", nnzR),      #Size _D_nnzR
            xlnz = vector_dc( SS$type, nrow+1),    #Size nrow

            snode = vector(SS$type, nrow), #17     #Size nrow
            xsuper = vector(SS$type, nrow+1),      #Size nrow
            cachesize = as.vector(cache, SS$type),
            ierr = vector(SS$type, 1),

            INTENT=cholstep.intent,
            NAOK = getOption("spam.NAOK"),
            PACKAGE = SS$package)

  if (verbose) cat("... . Size of R-Fortran transfered object: ", printSize( z),").\n", sep='')

  if(z$ierr == 1) stop("Singularity problem when calculating the Cholesky factor.")
  if(z$ierr == 6) stop("Inconsitency in the input",call.=FALSE)

  while( z$ierr>1) {
    if (verbose) cat("Not enough memory allocated. (Starting another round ...\n")    
    if(z$ierr == 4) {
      tmp <- ceiling(nnzR*getOption("spam.cholincreasefactor")[1])
      warning("Increased 'nnzR' with 'NgPeyton' method\n",
                    "(currently set to ",tmp," from ",nnzR,")",call.=FALSE)
      nnzR <- tmp
    }
    if(z$ierr == 5) {
      tmp <- ceiling(nnzcolindices*getOption("spam.cholincreasefactor")[2])
      warning("Increased 'nnzcolindices' with 'NgPeyton' method\n",
         "(currently set to ",tmp," from ",nnzcolindices,")",call.=FALSE)
      nnzcolindices <- tmp
    }
    if(force64 || nnzcolindices  > 2147483647) {
        SS <- .format64()
    } else {
        SS <- .format32
    }
    cholstep.signature <- rep(SS$signature, 20)
    cholstep.signature[c(3,15)] <- "double"


    z <- .C64("cholstepwise",
              SIGNATURE=cholstep.signature,

              nrow = as.vector(nrow, mode=SS$type),
              nnzA = as.vector(nnzA, SS$type),
              d =  x@entries,
              jd = as.vector(x@colindices, SS$type),

              id = as.vector(x@rowpointers, SS$type),
              doperm = as.vector(doperm,    SS$type),
              invp = vector(SS$type, nrow),
              perm = pivot,

              nnzlindx = vector(SS$type,1),   # length of colindices
              nnzcolindices = as.vector(nnzcolindices, SS$type),
              lindx = vector(SS$type, nnzcolindices),  # one backup!!
              xlindx = vector(SS$type, nrow+1),   # col pointers. will be shorter

              nsuper = vector(SS$type, 1),        #  number of supernodes
              nnzR = as.vector(nnzR, SS$type),    #  allocation for elements in the factor
              lnz = vector("double",nnzR),        # factor 
              xlnz = vector(SS$type, nrow+1),     # row pointers of factor

              snode = vector(SS$type, nrow),
              xsuper = vector(SS$type, nrow+1),
              cachesize = as.vector(cache, SS$type),
              ierr = vector(SS$type, 1),

              INTENT=cholstep.intent,
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)

    if (verbose) cat("... . Size of R-Fortran transfered object: ", printSize( z),").\n", sep='')
    if(z$ierr == 1) stop("Singularity problem or matrix not positive definite", call.=FALSE)
  }
  z  
}


########################################################################

