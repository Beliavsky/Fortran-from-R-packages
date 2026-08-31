# HEADER ####################################################
# This is file spam/R/toepliz.R.                            #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



########################################################################

"circulant.spam" <- function(x, n=NULL, eps = getOption("spam.eps"))
{
  if (!(is.vector(x)|is.list(x)) )
    stop("'x' is not a vector or a list")

  force64 <- getOption("spam.force64")

  if( is.list(x)) {
    if (!identical(length(x),2L))
      stop("Argument 'x' needs to be a list with two elements")
    if (is.null(n))
      stop("'n' needs to be given")
    ind <- x[[1]]
    x <- x[[2]]
    sel <- (ind <= n)&(abs(x)>eps)
    ind <- ind[sel]
    x <- x[sel]
   
  }else{
    n <- length(x)
    ind <- (1:n)[abs(x) > eps]
    x <- x[ind]
  }
  
  n <- as.integer(n)
  len <- as.integer(length( ind)[1]) # see ?length@value
  if(identical(len,0))
    return(.newSpam(
        # rowpointers = c(1, rep_len64(2, n)), 
        dimension = c(n, n),
        force64 = force64))
#      subroutine circulant(nrow,len, x,j, a,ja,ia)
  nz <- n*len
  ## z <- .Fortran("circulant",
  ##               as.integer(n),
  ##               as.integer(len),
  ##               as.double(x),
  ##               as.integer(ind),
  ##               entries= vector("double", nz),
  ##               colindices = vector("integer", nz),
  ##               rowpointers = vector("integer",  n + 1),
  ##               NAOK = getOption("spam.NAOK"),
  ##               PACKAGE = "spam")
  if(force64 || nz > 2147483647 || n+1 > 2147483647)
      SS <- .format64()
  else
      SS <- .format32
  z <- .C64("circulant",
             SIGNATURE = c(SS$signature, SS$signature, "double", SS$signature,
                 "double", SS$signature, SS$signature),
             
             n,
             len,
             x,
             ind,
            
             entries = vector_dc( "double", nz),
             colindices = vector_dc( SS$type, nz),
             rowpointers = vector_dc( SS$type,  n + 1),

             INTENT = c("r", "r", "r", "r",
                 "rw", "rw", "rw" ),
             NAOK = getOption("spam.NAOK"),
             PACKAGE = SS$package)

                
  ## newx <- new("spam")
  ## slot(newx, "entries", check = FALSE) <- z$entries
  ## slot(newx, "colindices", check = FALSE) <- z$colindices
  ## slot(newx, "rowpointers", check = FALSE) <- z$rowpointers
  ## slot(newx, "dimension", check = FALSE) <- c(n, n)
  ## return(newx)
  return(.newSpam(
      entries = z$entries,
      colindices = z$colindices,
      rowpointers = z$rowpointers,
      dimension = c(n,n),
      force64 = force64
      ))
}

toeplitz.spam <- function(x,y=NULL, eps = getOption("spam.eps"))
{
  force64 <- getOption("spam.force64")

  if (!is.vector(x)) 
    stop("'x' is not a vector")
  n <- length(x)

  if (!is.null(y)){
    if (!identical(length(y),n))
      stop("Length of 'y' and 'x' do not match")
    fullx <- c(rev(y[-1]),x)
  } else {
    fullx <- c(rev(x[-1]),x)
  }

  ind <- (1:(2*n-1))[abs(fullx) > eps]
  fullx <- fullx[ind]

 
  n <- as.integer(n)
  len <- as.integer(length( ind)[1]) # see ?length@value
  if(identical(len,0L)){
      ## print("degenerate")
      return(.newSpam(
          # rowpointers =  c(1, rep_len64(2, n)),
          dimension = c(n, n),
          force64 = force64))
    ## return(new("spam", rowpointers = c(1L, rep.int(2L, n)), 
    ##            dimension = as.integer(c(n, n))))
  }
#      subroutine toeplitz(nrow,len, x,j, a,ja,ia,kk)
  nz <- n*len
    ## z <- .Fortran("toeplitz",
    ##             as.integer(n),
    ##             as.integer(len),
    ##             as.double(fullx),
    ##             as.integer(ind),
    ##             entries= vector("double", nz),
    ##             colindices = vector("integer", nz),
    ##             rowpointers = vector("integer",  n + 1),
    ##             nnz=as.integer(1),
    ##             NAOK = getOption("spam.NAOK"), PACKAGE = "spam")
  if(force64 || n+1 > 2147483647 || nz > 2147483647 )
      SS <- .format64()
  else
      SS <- .format32
  
  z <- .C64("toeplitz",
            SIGNATURE = c(SS$signature, SS$signature, "double", SS$signature,
                "double", SS$signature, SS$signature, SS$signature),

            n,
            len,
            fullx,
            ind,
            
            entries = vector_dc("double", nz),
            colindices = vector_dc( SS$type, nz),
            rowpointers = vector_dc( SS$type,  n + 1),
            nnz = 1,

            INTENT = c("r", "r", "r", "r",
                "w", "w", "w", "rw"),
            NAOK = getOption("spam.NAOK"),
            PACKAGE = SS$package)

                
  ## newx <- new("spam")
  ## slot(newx, "entries", check = FALSE) <- z$entries[1:z$nnz]
  ## slot(newx, "colindices", check = FALSE) <- z$colindices[1:z$nnz]
  ## slot(newx, "rowpointers", check = FALSE) <- z$rowpointers
  ## slot(newx, "dimension", check = FALSE) <- c(n, n)
  ## return(newx)
  return(.newSpam(
      entries = z$entries[1:z$nnz],
      colindices = z$colindices[1:z$nnz],
      rowpointers = z$rowpointers,
      dimension = c(n, n),
      force64 = force64))
}
