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
#    chol.spam, chol.update and triangular solves.
# Associated S4 elements and determinant are in different files.
#
#
########################################################################
########################################################################






# TODO below, check with solve.spam and triangular solves for use of vector_dc
# determinant.spam from "rw" to "w" signature.


chol.spam <- function(x, pivot="MMD", method="NgPeyton", memory=list(),
                      eps=getOption("spam.eps"), Rstruct=NULL, ..., verbose=FALSE){

  if (is(Rstruct,"spam.chol.NgPeyton")) invisible( update.spam.chol.NgPeyton(Rstruct,x,...))

    if (verbose) {
        timeused <- proc.time()
        cat("\nStarting Cholesky factorization with method '",method,"'.\n")
    }
  if (method != "NgPeyton")
      warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'",
                     method), domain = NA)

  z <- .chol.spam.basis(x, pivot=pivot, memory=memory, eps=eps, verbose=verbose)

  nnzRfinal <- z$xlnz[length(z$xlnz)]-1  # reduce to final size

  newx <- new("spam.chol.NgPeyton")
  slot(newx, "entries", check=FALSE) <- z$lnz[1:nnzRfinal]
  slot(newx, "colindices", check=FALSE) <- z$lindx[1:z$nnzlindx]
  slot(newx, "colpointers", check=FALSE) <- z$xlindx[1:(z$nsuper+1)]
  slot(newx, "rowpointers", check=FALSE) <- z$xlnz
  slot(newx, "dimension", check=FALSE) <- c(z$nrow, z$nrow)
  slot(newx, "pivot", check=FALSE) <- z$perm
  slot(newx, "invpivot", check=FALSE) <- z$invp
  slot(newx, "supernodes", check=FALSE) <- z$xsuper[1:(z$nsuper+1)]
  slot(newx, "snmember", check=FALSE) <- z$snode
  slot(newx, "memory", check=FALSE) <- c(z$nnzcolindices, z$nnzR, z$cachesize)
  slot(newx, "nnzA", check=FALSE) <- z$nnzA

  if (verbose) {
      cat("Final size of Cholesky factor (of class 'spam.chol.NgPeyton'): ",
          printSize(newx), ".\n")
      cat("Finished factorization (total time used:", (proc.time()-timeused)[1], "s).\n\n")
  }

  invisible(newx)
}



########################################################################


update.spam.chol.NgPeyton <- function(object,x,...){
    ## print("update.spam.chol.NgPeyton")
  nrow <- object@dimension[1]
  if (!is.spam(x))
    stop("Covariance should be a 'spam' object.")
  if ((x@rowpointers[nrow+1]-1) != object@nnzA)
    stop("Updated covariance entries do not match length of original one.")

    if( getOption("spam.force64") || .format.spam(object)$package != "spam" || .format.spam(x)$package != "spam" )
        SS <- .format64()
    else
        SS <- .format32

    u <- .C64("updatefactor",
     ##          subroutine updatefactor( m,nnzd,
     ## &     d,jd,id, invp,perm,
     ## &                lindx,xlindx, nsuper,lnz,xlnz,
     ## &                snode, xsuper,
     ## &                cachesize,ierr)
              SIGNATURE = c(SS$signature, SS$signature,
                            "double", SS$signature, SS$signature,
                            SS$signature, SS$signature, SS$signature, SS$signature, SS$signature,
                            "double", SS$signature, SS$signature, SS$signature, SS$signature, SS$signature ),
              nrow,
              object@nnzA,

              d =  x@entries,
              jd = x@colindices,
              id = x@rowpointers,

              object@invpivot,
              object@pivot,
              lindx = object@colindices,
              xlindx = object@colpointers,
              nsuper = length(object@supernodes)-1,

              entries = vector_dc( "double", length(object@entries)), #lnz
              rowpointers = object@rowpointers,#xlnz
              snode = object@snmember,
              xsuper = object@supernodes,
              cachesize = object@memory[3],
              ierr = 0,

              ## INTENT = c("r", "r",
              ##            "r", "r", "r",
              ##            "r", "r", "r", "r", "r",
              ##            "rw", "rw", "r", "r", "rw", "rw"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)

  if(u$ierr>1) stop("Internal error in 'update.spam.chol.NgPeyton' code ", u$ierr,call.=FALSE)

  if(u$ierr == 1) {
    if (getOption("spam.cholupdatesingular") == "null")
      return(NULL)
    else if (getOption("spam.cholupdatesingular") == "error")
      stop("Singularity problem when updating a Cholesky Factor.")
    else if (getOption("spam.cholupdatesingular") == "warning")
      warning("Singularity problem when updating a Cholesky Factor.\n'object' not updated.")
    else
      stop("'cholupdatesingular' should be 'error', 'null' or 'warning'.")
  }  else {
    slot(object, "entries", check=FALSE) <- u$entries
  }
  invisible(object)
}




solve.spam <- function (a, b,  Rstruct=NULL, ...) {
  nrow <- a@dimension[1]
  ncol <- a@dimension[2]
  if (ncol != nrow)      stop("only square matrices can be inverted")

  if (missing(b)) {
    b <- diag(1, ncol)
  }  else {
    if(!is.matrix(b)) b <- as.matrix(b)
  }
  p <- dim(b)[2]
  if(nrow!=dim(b)[1])stop("'b' must be compatible with 'a'")

  # if we have a spam matrix, we calculate the Cholesky factor
  if (is(a,"spam"))
    if (is(Rstruct, "spam.chol.NgPeyton"))
        a <- update.spam.chol.NgPeyton(Rstruct, a, ...)
    else a <- chol.spam(a, ...)


  if (is(a,"spam.chol.NgPeyton")) {
      # The following is a fast way to perform:
      #     z <- backsolve(a,forwardsolve( t(a),b))
    nsuper <- as.integer( length(a@supernodes)-1)
      if( getOption("spam.force64") || .format.spam(a)$package != "spam" )
          SS <- .format64()
      else
          SS <- .format32

      z <- .C64("backsolves",
                SIGNATURE = c(rep(SS$signature,5),"double", rep(SS$signature, 4),
                              rep("double",3)),

                m = nrow,
                nsuper,
                p,
                a@colindices,
                a@colpointers,
                a@entries,  # "double"
                a@rowpointers,
                a@invpivot,
                a@pivot,
                a@supernodes,

                vector_dc("double",nrow),
                sol = vector_dc("double",nrow*p),
                as.vector(b,"double"),

                INTENT = c( rep( "r", 10), rep( "rw", 3)),
                NAOK = getOption("spam.NAOK"),
                PACKAGE = SS$package )$sol
  } else z <- backsolve(a, forwardsolve( t(a),b))
    # see the helpfile for a comment about the 't(a)' construct.

  if ( p!=1)    dim(z) <- c(nrow,p)
  return( z)
}

chol2inv.spam <- function (x, ...) {
    ## print("chol2inv.spam")
    nrow <- x@dimension[1]

    if (is(x,"spam.chol.NgPeyton")) {
        y <- vector("double",nrow*nrow)
        y[1L + 0L:(nrow - 1L) * (nrow + 1L)] <- 1.0

        if( getOption("spam.force64") || .format.spam(x)$package !="spam" )
            SS <- .format64()
        else
            SS <- .format32


        z <- .C64("backsolves",
          ##         subroutine backsolves(m,nsuper,nrhs,lindx,xlindx,lnz,
     ## &                   xlnz,invp,perm,xsuper,newrhs,sol,b)
                  SIGNATURE = c( rep( SS$signature, 5), "double", rep( SS$signature, 4),
                                rep( "double", 3)),

                  m = nrow, #r
                  length(x@supernodes)-1, #r
                  nrow, #r
                  x@colindices,
                  x@colpointers,
                  x@entries,   # "double"
                  x@rowpointers,
                  x@invpivot, #r
                  x@pivot, #r
                  x@supernodes, #r

                  vector_dc("double",nrow), #rw
                  sol = vector_dc("double",nrow*nrow), #w
                  y, #r

                  INTENT = c( rep( "r", 10),
                             "rw", "w", "r"),
                  NAOK = getOption("spam.NAOK"),
                  PACKAGE = SS$package)$sol

        dim(z) <- c(nrow,nrow)
    } else z <- backsolve.spam(x, forwardsolve.spam( t(x), diag(nrow)))
    return( z)
}

backsolve.spam <- function(r, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
# r: spam.chol.NgPeyton structure as returned by chol.spam or a spam object
# x: rhs a vector or a matrix in dense form
# dimensions:  ( m x n) ( n x p)

    if( getOption("spam.force64") || .format.spam(r)$package !="spam" )
        SS <- .format64()
    else
        SS <- .format32


  m <- r@dimension[1]
  if(is.vector(x)) {
    n <- length(x)
    p <- 1
  } else {
    x <- as.matrix(x)
    n <- dim(x)[1]
    p <- dim(x)[2]
  }

  # we separate between "spam.chol.NgPeyton" and "spam"
  if (is(r,"spam.chol.NgPeyton")) {
    if (n!=m) stop("Cholesky factor 'r' not compatible with 'x'")
    nsuper <- length(r@supernodes)-1
    if (!getOption("spam.dopivoting")) {

      z <- .C64("backsolvef",
                SIGNATURE=c(rep(SS$signature, 5), "double", SS$signature, SS$signature,
                            "double"),

                m,
                nsuper,
                p,
                r@colindices,

                r@colpointers,
                r@entries,   # "double"
                r@rowpointers,
                r@supernodes,

                sol = vector_dc("double",m*p),

                INTENT=c("r", "r", "r", "r", "r", "r", "r", "r", "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE=SS$package)$sol

    }else{
      z <- .C64("pivotbacksolve",
                SIGNATURE=c(rep(SS$signature, 5), "double", rep(SS$signature, 4),
                            "double", "double", "double"),

                m,
                nsuper,
                p,
                r@colindices,

                r@colpointers,
                r@entries, #  "double"
                r@rowpointers,
                r@invpivot,

                r@pivot,
                r@supernodes,
                vector("double",m),
                sol = vector_dc("double",m*p),
                x,            #  "double"

                INTENT=c("r", "r", "r", "r",
                    "r", "r", "r", "r",
                    "r", "r", "r", "w",
                    "r"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE=SS$package)$sol
    }
  } else {
    if (n!=m) stop("Triangular matrix 'r' not compatible with 'x'")
    # solve R sol = x
    z <- .C64("spamback",
              SIGNATURE=c(SS$signature, SS$signature, "double", "double",
                          "double", SS$signature, SS$signature),
              m=m,
              unused=p,
              sol = vector_dc("double",m*p),
              x=x,

              al=r@entries,
              jal=r@colindices,
              ial=r@rowpointers,

              INTENT=c("rw", "r", "w", "r", "r", "r", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE=SS$package)
    if (z$m<0) stop(gettextf("singular matrix in 'backsolve'. Last zero in diagonal [%d]",
            -z$m), domain = NA)
     else z <- z$sol

  }

  if (p>1) dim(z) <- c(m,p)
  return(z)
}

forwardsolve.spam <- function(l, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
#  l: spam.chol.NgPeyton structure as returned by chol.spam
#         or an ordinary lower triangular spam matrix
#  x: rhs a vector a matrix in dense form
#  dimensions:  ( m x n) ( n x p)
#  if (!any(is.null(c(upper.tri,k,transpose ))))
#    warning("'k', 'upper.tri' and 'transpose' argument do not have any effect here")
    if( getOption("spam.force64") || .format.spam(l)$package !="spam" )
        SS <- .format64()
    else
        SS <- .format32

  m <- l@dimension[1]
  if(is.vector(x)) {
    n <- length(x)
    p <- 1L
  } else {
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    p <- dim(x)[2]
  }

  # we separate between "spam.chol.NgPeyton" and "spam"
  if (is(l,"spam.chol.NgPeyton")) {
    if(n!=m) stop("Cholesky factor 'l' not compatible with 'x'")
    nsuper <- length(l@supernodes)-1
    if (!getOption("spam.dopivoting")) {
      z <- .C64("forwardsolvef",
                SIGNATURE=c(SS$signature, SS$signature, SS$signature, SS$signature,
                            SS$signature, "double", SS$signature, SS$signature,
                            "double"),

                m,
                nsuper,
                p,
                l@colindices,

                l@colpointers,
                l@entries,
                l@rowpointers,
                l@supernodes,

                sol = vector("double",m*p),

                INTENT=c("r", "r", "r", "r",
                    "r", "r", "r", "r",
                    "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE=SS$package)$sol
    }else{
      z <- .C64("pivotforwardsolve",
                SIGNATURE=c(rep(SS$signature,5), "double", rep(SS$signature, 4),
                            "double", "double", "double"),
                m,
                nsuper,
                p,
                l@colindices,

                l@colpointers,
                l@entries,     # "double"
                l@rowpointers,
                l@invpivot,

                l@pivot,
                l@supernodes,
                vector("double",m),
                sol = vector_dc("double",m*p),

                x,

                INTENT=c("r", "r", "r", "r",
                    "r", "r", "r", "r",
                    "r", "r", "r", "w",
                    "r"),
                PACKAGE=SS$package)$sol
    }
  } else {
    if (n!=m) stop("Triangular matrix 'l' not compatible with 'x'")
    # solve L sol = x
    z <- .C64("spamforward",
              SIGNATURE=c(SS$signature, SS$signature, "double", "double",
                          "double", SS$signature, SS$signature),

              m=m,
              p,
              sol = vector_dc("double",m*p),
              x=x,

              al=l@entries,
              jal=l@colindices,
              ial=l@rowpointers,

              INTENT=c("rw", "r", "w", "r", "r", "r", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE=SS$package)
    if (z$m<0) stop(gettextf("singular matrix in 'forwardsolve'. First zero in diagonal [%d]", -z$m), domain = NA)
    else z <- z$sol
  }
  if (p>1)
    dim(z) <- c(m,p)
  return(z)
}




########################################################################

