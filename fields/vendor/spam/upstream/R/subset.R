# HEADER ####################################################
# This is file spam/R/subset.R.                             #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################




# SUBSETTING
##########################################################################################

# notice the drop catch...
#   I don"t know the best and official way, but it works as it is here...

setMethod("[", signature(x = "spam",
			 i = "missing", j = "missing", drop = "ANY"),
	  function (x, i, j,..., drop) { # cat("missmiss")
           x})

setMethod("[",signature(x="spam",i="vector",j="missing", drop = "logical"),
	  function (x, i, j,..., drop) {  #cat("   log call was", deparse(match.call()), "\n")
            if (nargs()==3) {
              subset_rows.spam(x, i,drop=drop)
            } else {
              subset_rows.spam(x, i,,drop=drop)
            }}
          )

setMethod("[",signature(x="spam",i="vector",j="missing", drop = "missing"),
	  function (x, i, j,..., drop) { #cat("    mis call was", deparse(match.call()), "\n")
            if (nargs()==2) {
              subset_rows.spam(x, i)
            } else {
              subset_rows.spam(x, i,)
              }})

setMethod("[",signature(x="spam",i="vector",j="vector", drop = "ANY"),
	  function (x, i, j,..., drop) { # cat("vecvec")
            subset.spam(x,rw=i,cl=j,drop=drop)} )

setMethod("[",signature(x="spam",i="missing",j="vector", drop = "ANY"),
	  function (x, i, j,...,drop) { # cat("missvec")
            subset.spam(x,rw=1:x@dimension[1],cl=j,drop=drop)} )

setMethod("[",signature(x="spam",i="matrix",j="missing", drop = "missing"),
	  function (x, i, j,..., drop) {subset.spam(x,rw=i) })

setMethod("[",signature(x="spam",i="matrix",j="missing", drop = "logical"),
	  function (x, i, j,..., drop) {subset.spam(x,rw=i,drop=drop) })

setMethod("[",signature(x="spam",i="matrix",j="matrix", drop = "ANY"),
	  function (x, i, j,..., drop) {subset.spam(x,rw=cbind(c(i),c(j)),drop=drop) })

setMethod("[",signature(x="spam",i="spam",j="missing", drop = "ANY"),
	  function (x, i, j,..., drop=getOption("spam.drop"))
{
                                        # drop is not implemented yet
    ## print("subset spam spam")
  dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]
  if ( i@dimension[1]>nrow | i@dimension[2]>ncol)
    stop("subscript out of bounds",call.=FALSE)
  ## z <- .Fortran("amask",
  ##               nrow=as.integer(nrow),
  ##               ncol=as.integer(ncol),
  ##               a=as.double(x@entries),
  ##               colindices=as.integer(x@colindices),
  ##               rowpointers=as.integer(x@rowpointers),
  ##               jmask=as.integer(i@colindices),
  ##               imask=as.integer(c(i@rowpointers,rep(i@rowpointers[length(i@rowpointers)],nrow+1-length(i@rowpointers)))),
  ##               c=as.double(x@entries),
  ##               jc=as.integer(x@colindices),
  ##               ic=as.integer(x@rowpointers),
  ##               iw=logical(ncol),
  ##               nzmax=as.integer(length(i@colindices)),
  ##               ierr=0L,
    ##               NAOK=getOption("spam.NAOK"),PACKAGE="spam") # some copying is required!!!!
    if(.format.spam(x)$package == "spam64")
        SS <- .format64()
    else
        SS <- .format32

    z <- .C64("amask",
             ##     subroutine amask (nrow,ncol,a,ja,ia,jmask,imask,
     ## *                  c,jc,ic,nzmax,ierr)
              SIGNATURE = c( SS$signature, SS$signature,
                            "double", SS$signature, SS$signature,
                            SS$signature, SS$signature,
                            "double", SS$signature, SS$signature,
                             SS$signature, SS$signature),

              nrow = nrow,
              ncol = ncol,

              a = x@entries,
              colindices = x@colindices,
              rowpointers = x@rowpointers,

              jmask = i@colindices,
              imask = c(i@rowpointers, rep_len64(i@rowpointers[length(i@rowpointers)], nrow+1-length(i@rowpointers))),

              c = x@entries,
              jc = x@colindices,
              ic = x@rowpointers,

              nzmax = length( i@colindices),
              ierr = 0,

              INTENT = c( "rw", "rw",
                         "r", "r", "r",
                         "r", "r",
                         "w", "w", "w",
                         "r", "w" ),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package) # some copying is required!!!!

  nz <- z$ic[nrow+1]-1
    if (nz==0) {
        ## print("nz = 0")
        return( numeric(0))
    }
    if (drop) {
        ## print("drop")
    ic <- unique( z$ic[1:(z$nrow+1)])
    dimx <- as.integer(c(length(ic)-1,max(z$jc[1:nz])))
    } else {
        ## print("notdrop")
    ic <-z$ic[1:(z$nrow+1)]
  }
  ## return(new("spam",entries=z$c[1:nz],colindices=z$jc[1:nz],rowpointers=ic,
    ##              dimension=dimx))
    return(.newSpam(
        entries = z$c[1:nz],
        colindices = z$jc[1:nz],
        rowpointers = ic,
        dimension = dimx))
}      )

setMethod("[", signature(x = "spam", i = "ANY", j = "ANY", drop = "ANY"),
	  function(x,i,j, drop)
          stop("Invalid or not-yet-implemented 'spam' subsetting"))

# the proper S3 subsetting causes problems...
#  "[.spam" <- function (x, rw, cl,drop=getOption("spam.drop")) {subset.spam(x,rw=rw,cl=cl,drop) }


"subset_rows.spam" <-
function (x, i, ..., drop=getOption("spam.drop"))
  # approach: we extract rows (nargs=2) or elements (nargs=3)
  #  i is a vector of integers or logical!
  # nargs idea from Matrix!
{
  nA <- nargs()+missing(drop)
#  cat("subset_rows.spam call was", deparse(match.call())," ",nargs(), " " , nA, "\n")
 dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]
  mini <- min(i, na.rm=TRUE)
  maxi <- max(i, na.rm=TRUE)
  if (mini<0 & maxi>0)
    stop("Negative and positive subscripts mixed")

  SS <- .format.spam(x)

  if(nA==3) { # extract elements
    if (is.logical(i)) {
      inefficiencywarning( "Logical subsetting may be inefficient, is this really what you want?",
                          prod(dimx))
      return(.C64("spamcsrdns",
                  SIGNATURE = c(SS$signature,
                                "double", SS$signature, SS$signature,
                                "double"),
                  nrow,

                  entries = x@entries,
                  colindices = x@colindices,
                  rowpointers = x@rowpointers,

                  res = vector_dc("double", prod(dimx)),
 #                 INTENT = c("r","r","r","r","w"), # not checked   !!!
                  NAOK=getOption("spam.NAOK"),PACKAGE = SS$package)$res[i])
    }

    if (mini<0) {
      inefficiencywarning( "Negative subsetting may be inefficient, is this really what you want?",
                          prod(dimx))
      return(.C64("spamcsrdns",
                  SIGNATURE = c(SS$signature,
                                "double", SS$signature, SS$signature,
                                "double"),
                  nrow,

                  entries=x@entries,
                  colindices=x@colindices,
                  rowpointers=x@rowpointers,

                  res=vector_dc("double",prod(dimx)),
#                  INTENT = c("r","r","r","r","w"),  # not checked   !!!
                  NAOK=getOption("spam.NAOK"),PACKAGE = SS$package)$res[i])
    }
    # eliminate zeros,
    # force too large to NA, keep NAs
    i <- i[i>0]
    ind <- !(is.na(i)|(i> (nrow*ncol)))
    ii <- i[ind]-1
    i <- ii %% nrow+1
    j <- ii %/% nrow+1
    nir <- length(i)
    z <- vector("double",length(ind))
    z[!ind] <- NA
    z[ind] <- .C64("getallelem",
                   SIGNATURE = c(SS$signature, SS$signature, SS$signature,
                                 "double", SS$signature, SS$signature,
                                 SS$signature, "double"),
                   nir,
                   i,
                   j,

                   x@entries,
                   x@colindices,
                   x@rowpointers,

                   iadd = vector_dc(SS$type,nir),
                   allelem = vector_dc("double",nir),

                   ## INTENT = c("r","r","r",
                   ##            "r","r","r",
                   ##            "w","w"),   # not checked   !!!
                   NAOK=getOption("spam.NAOK"),PACKAGE=SS$package)$allelem
# getallelem(nir,ir,jr,a,ja,ia,alliadd,allelem)
    return(z)
  }
  if(nA==4) {

    if (is.logical(i)) {  # logical
      if( length(i) > nrow)  stop("(subscript) logical subscript too long",call.=FALSE)

      i <- seq_len( nrow)[i]
    }  else {

      i <- i[i!=0]     # eliminate zero lines

    if (maxi>x@dimension[1])
      stop("subscript out of bounds",call.=FALSE)

      # negative values:
      if ( maxi <= 0 )       i <- seq_len( nrow)[i]
    }


    ni <- length(i)
    if (ni==0) return(numeric(0))   # zero elements...

    if (any(is.na(i))) {
      stop("'NA's in subsetting vector have been eliminated.")
#      i <- i[!is.na(i)]
    }


    nz <- (sum(x@rowpointers[i+1]-x@rowpointers[i]))
    if (nz==0) {#trap zero matrix
      if (drop==TRUE && (ni==1 || ncol==1)) return( vector("double",max(ni,ncol)))
      else
        return(new("spam",rowpointers=c(1L,rep.int(2L,ni )),
                   dimension = c(ni,ncol)))
    }  else {
 #          subroutine getlines(a,ja,ia, ni, i, bnz, b,jb,ib)
     z <- .C64("getlines",
               SIGNATURE = c("double", SS$signature, SS$signature,
                             SS$signature,  SS$signature, SS$signature,
                             "double", SS$signature, SS$signature),
               x@entries, x@colindices, x@rowpointers,

               ni,
               i,
               newnz=nz,

               entries=vector_dc("double", nz),
               colindices=vector_dc(SS$type, nz),
               rowpointers=vector_dc(SS$typ, ni+1),

               NAOK=getOption("spam.NAOK"),
               PACKAGE=SS$package)
  #    print(c(nz,z$newni,is.integer(nz), is.integer(z$newni),z$newni!=ni))
     if(z$newnz!=nz) stop(gettextf("Subsetting error, please report %d, %d",z$newnz,nz))
    }
#    print(c(drop,ni,ncol,(drop==TRUE && (ni==1 || ncol==1) )))
    if (drop==TRUE && (ni==1 || ncol==1))
      # this is essentially a c() call
      return(.C64("spamcsrdns",
                  SIGNATURE = c(SS$signature,
                                "double", SS$signature, SS$signature,
                                "double"),
                  nrow=ni,

                  entries=z$entries,
                  colindices=z$colindices,
                  rowpointers=z$rowpointers,

                  res=vector_dc("double",prod(ni,ncol)),
                  NAOK=getOption("spam.NAOK"),
                  PACKAGE = SS$package)$res)
    else {
        return(.newSpam(
            entries=z$entries,
            colindices=z$colindices,
            rowpointers=z$rowpointers,
            dimension=c(ni,ncol),
            force64=getOption("spam.force64")))
    }

  }
  stop("incorrect number of dimensions")
}

## subset_rows.spam <- function(x, i, drop){

##   dimx <- x@dimension
##   nrow <- dimx[1]
##   ncol <- dimx[2]

##     ni <- length(i)
##     if (ni==0) return(numeric(0))   # zero elements...

##     if (any(is.na(i))) {
##       stop(""NA"s in subsetting vector not allowed.")
##     }


##     nz <- sum(x@rowpointers[i+1]-x@rowpointers[i])
##     if (nz==0) { #trap zero matrix
##       if (drop==TRUE && (ni==1 || ncol==1))
##         return( vector("double",max(ni,ncol)))
##       else
##         return(new("spam",rowpointers=c(1,rep(2L,ni )),
##                    dimension = c(ni,ncol)))
##     }

##     SS <- .format.spam(x)

##      z <- .C64("getlines",
##                SIGNATURE=c("double", SS$signature, SS$signature, SS$signature,
##                    SS$signature, SS$signature, "double", SS$signature,
##                    SS$signature),

##                x@entries,
##                x@colindices,
##                x@rowpointers,
##                ni,

##                i, #int64
##                newnz=nz,
##                entries=vector_dc("double",nz),
##                colindices=vector_dc(SS$type,nz),

##                rowpointers=vector_dc(SS$type,ni+1),

##                INTENT=c("r", "r", "r", "r",
##                    "r", "rw", "w", "w",
##                    "w"),
##                NAOK = getOption("spam.NAOK"),
##                PACKAGE=SS$package)
##      if(z$newnz!=nz) stop(gettextf("Subsetting error, please report %d, %d",z$newnz,nz))

##     # --- CHANGED ---

##   newx <- .newSpam(entries=z$entries,
##                    colindices=z$colindices,
##                    rowpointers=z$rowpointers,
##                    dimension=c(ni,ncol))

##     if (drop==TRUE && (ni==1 || ncol==1)) {

##         x <- newx
##         SS <- .format.spam(x)

##         dimx <- x@dimension

##         result <- .C64("spamcsrdns",
##                        SIGNATURE=c(SS$signature, "double", SS$signature, SS$signature, "double"),
##                        nrow=dimx[1],
##                        entries=x@entries,
##                        colindices=x@colindices,
##                        rowpointers=x@rowpointers,
##                        res=vector("double",prod(dimx)),  # TODO: Expects zeros, because it
##                                         # only overwrites the non-zero elements! (use RW)

##                        INTENT=c("r", "r", "r", "r", "rw"),
##                        NAOK = getOption("spam.NAOK"),
##                        PACKAGE = SS$package
##                        )
##         return(result$res)
##     }else {
##       return(newx)
##     }

##   stop("incorrect number of dimensions")
## }

"subset.spam" <- function (x,rw,cl,...,drop=getOption("spam.drop"))
{
  # we separate into cases where:
  # (A) rw matrix:
  #     1: logical: transformation to spam and extract structure
  #     2: two column matrix: extract (i,j) as given by the lines.
  #     3: all else extract   x[ c( rw)]
  # (B) rw and cl one element: ((i,j)
  # (C) rw and cl vectors:  (i1:i2,j1:j2)               [i1<=i2, j1<=j2]
  #                         (c(i1,...,ii),c(j1,...,jj)) [arbitrary block]
#  if (missing(drop)) drop <- getOption("spam.drop")
#  print(drop)
  dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]

  SS <- .format.spam(x)

  if (is.matrix(rw)) {
    if (is.logical(rw)) {
      return( x[as.spam.matrix(rw)] )
    }
    if (dim(rw)[2]==2) {
      ir <- rw[,1]
      jr <- rw[,2]
    } else  {
      ir <- c(rw-1) %% nrow + 1
      jr <- c(rw-1) %/% nrow + 1
    }
    if ( (min(ir)<1)|(max(ir)>x@dimension[1])|(min(jr)<1)|(max(jr)>x@dimension[2]))
      stop("subscript out of bounds",call.=FALSE)
    nir <- length(ir)
    return(.C64("getallelem",
                SIGNATURE = c(SS$signature, SS$signature, SS$signature,
                              "double", SS$signature, SS$signature,
                              SS$signature,"double"),
                nir,
                ir,
                jr,

                x@entries, x@colindices, x@rowpointers,

                vector_dc(SS$type, nir),
                allelem=vector_dc("double",nir),
                NAOK=getOption("spam.NAOK"), PACKAGE=SS$package)$allelem)

  }
  # negative values:
  if ( max(rw)<0 )       rw <- seq_len( nrow)[rw]
  if ( max(cl)<0 )       cl <- seq_len( ncol)[cl]

  # logical
  if (is.logical(rw))    rw <- seq_len( nrow)[rw]
  if (is.logical(cl))    cl <- seq_len( ncol)[cl]

  if (length(cl)==0) stop("You should subset at least one element for the columns",call.=FALSE)
  if (length(rw)==0) stop("You should subset at least one element for the rows",call.=FALSE)

  if ( (min(rw)<1)|(max(rw)>x@dimension[1])|(min(cl)<1)|(max(cl)>x@dimension[2]))
    stop("subscript out of bounds",call.=FALSE)

  if (length(rw)==1 & length(cl)==1){
                                        # function to extract only one element
    return(.C64("getelem",
                SIGNATURE = c(SS$signature, SS$signature,
                              "double", SS$signature, SS$signature,
                              SS$signature, "double"),
                rw,
                cl,

                x@entries, x@colindices, x@rowpointers,

                iadd=vector(SS$type, 1),
                elem=vector("double", 1),
                PACKAGE=SS$package)$elem)
  }
  if (is.vector(rw) && is.vector(cl)) {
    nrw <- length(rw)   # length returns an integer, so is a product therof
    ncl <- length(cl)
    diffrw <- diff(rw)
    diffcl <- diff(cl)
    nz <- ( min( (1+sum(diff(sort(rw))==0))*(1+sum(diff(sort(cl))==0))*
                length(x@entries), prod(nrw,ncl)))  # very pessimistic
    if (all(diffrw==1) & all(diffcl==1)) {
      z <- .C64("submat",
                SIGNATURE = c(SS$signature,
                              SS$signature, SS$signature, SS$signature, SS$signature,
                              "double", SS$signature, SS$signature,
                              SS$signature, SS$signature,
                              "double", SS$signature, SS$signature),
                job=1, # need values as well

                i1=rw[1],
                i2=rw[nrw],
                j1=cl[1],
                j2=cl[ncl],

                x@entries, x@colindices, x@rowpointers,

                nr=0,
                nc=0,

                entries=vector_dc("double",nz),
                colindices=vector_dc(SS$type,nz),
                rowpointers=vector_dc(SS$type,nrw+1),

                NAOK=getOption("spam.NAOK"),PACKAGE = SS$package)
      nz <- z$rowpointers[z$nr+1]-1
    } else {
      z <- .C64("getblock",
                SIGNATURE = c("double", SS$signature, SS$signature,
                              SS$signature, SS$signature, SS$signature, SS$signature, SS$signature,
                              "double", SS$signature, SS$signature),
                x@entries, x@colindices, x@rowpointers,

                nr=nrw,
                rw,
                nc=ncl,
                cl,
                nz=nz,

                entries=vector("double",nz),
                colindices=vector(SS$type,nz),rowpointers=vector(SS$type,nrw+1),
                NAOK=getOption("spam.NAOK"),
                PACKAGE = SS$package)
      nz <- z$nz
    }
    if (nz==0) {#trap zero matrix
      if (drop==TRUE && (z$nr==1 || z$nc==1)) return( vector("double",max(z$nr,z$nc)))
      else
        return(new("spam",rowpointers=c(1L,rep.int(2L,z$nr )),
                   dimension = c(z$nr,z$nc)))
    }

    if (drop==TRUE && (z$nr==1 || z$nc==1))
      # this is essentially a c() call
      return(.C64("spamcsrdns",
                  SIGNATURE = c(SS$signature, "double",
                                SS$signature, SS$signature, "double"),
                  nrow=z$nr,
                  entries=z$entries[1:nz],

                  colindices=z$colindices[1:nz],
                  rowpointers=z$rowpointers[1:(z$nr+1)],

                  res=vector_dc("double", prod(z$nr,z$nc)),
                  NAOK=getOption("spam.NAOK"), PACKAGE = SS$package)$res)
    else {
        return(.newSpam(
            entries=z$entries[1:nz],
            colindices=z$colindices[1:nz],
            rowpointers=z$rowpointers[1:(z$nr+1)],
            dimension=c(z$nr,z$nc),
            force64=getOption("spam.force64")))
    }

  }
  stop("invalid or not-yet-implemented 'spam' subsetting")
}



#subset.rows.spam <- function(x, i, ..., drop=getOption("spam.drop")) {
subset.rows.spam <- function(...) {
    .Defunct('spam:::subset_rows.spam', package = 'spam', msg = "'subset.rows.spam' is defunct. Use 'spam:::subset_rows.spam' instead.\n")
}
