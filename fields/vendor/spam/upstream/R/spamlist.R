# HEADER ####################################################
# This is file spam/R/spamlist.R.                           #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



as.spam.list <- function(x, eps = getOption("spam.eps")) {
    spam.list( x,  eps=eps)
}

spam.list <-  function(x, nrow, ncol, eps = getOption("spam.eps")) {
    force64 <- getOption("spam.force64")
    
    if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
    if (!is.list(x)|(length(x)<2)|(length(x)>3))
        stop("Argument 'x' needs to be a list with two or three elements")
                                        # two cases: list of length
                                        # -  two (matrix with two columns called ind* and the elements)
                                        # -  three (each one column called i*, j*.   

    if (identical(length(x),2L)) {
        indnr <- pmatch("ind",names(x)) 
        if (is.na(indnr)) stop("Argument 'x' needs an element called 'indices'")
        elenr <- ifelse( identical( indnr,1L), 2L, 1L)
        
        nz <- length( x[[elenr]])

        dimx <- dim(x[[indnr]])
        if (is.null(dimx)||(dimx[2] != 2))  stop("Indices should have two columns")
        if (dimx[1] != nz) stop("Number of indices does not match with number of elements")
        
        ir <- as.integer(x[[indnr]][,1])
        jc <- as.integer(x[[indnr]][,2])

        if(force64 || length(x[[elenr]]) > 2147483646)
            SS <- .format64()
        else
            SS <- .format32
        
    } else {
        inr <- pmatch("i",names(x)) 
        jnr <- pmatch("j",names(x))
        
        if (is.na(inr)||is.na(jnr)) stop("Argument 'x' needs elements called 'i' and 'j'")
        elenr <- c(1:3)[-c(inr,jnr)]
        nz <- length( x[[elenr]])
        
        ir <- as.integer(x[[inr]])
        jc <- as.integer(x[[jnr]])

        if ((length(ir) != nz)||(length(jc) != nz))
            stop("Number of indices does not match with number of elements")

        if(force64 || length(x[[elenr]]) > 2147483646)
            SS <- .format64()
        else
            SS <- .format32
    }
    
    if (nz == 0)
        return(.newSpam(
            # rowpointers = c(1,rep_len64(2, nrow)),
            dimension = c(nrow,ncol),
            force64 = force64))
    if (any( ir <= 0) || any( jc <= 0))
        stop("Indices need to be positive")
    if (any(!is.finite(x[[elenr]]))) {
        warning("'NA/NaN/Inf' coerced to zero")
        x[[elenr]][!is.finite(x[[elenr]])] <- 0
    }
    nrow <- as.integer(ifelse(missing(nrow),max(ir),nrow))
    ncol <- as.integer(ifelse(missing(ncol),max(jc),ncol))
    ## z <- .Fortran(ifelse(toupper(getOption("spam.listmethod")=="PE"),"triplet3csr","triplet2csr"),
    ##               nrow=as.integer(nrow), ncol=as.integer(ncol),
    ##               nz=as.integer(nz),
    ##               as.double(x[[elenr]]),as.integer(ir),as.integer(jc),
    ##               entries=vector("double",nz),
    ##               colindices=vector("integer",nz),
    ##               rowpointers=vector("integer",nrow+1), as.double(eps),
    ##               NAOK=TRUE, PACKAGE = "spam"
    ##               )
    
    z <- .C64(ifelse(toupper(getOption("spam.listmethod")=="PE"),"triplet3csr","triplet2csr"),
              ## subroutine triplet3csr(nrow,ncol,nnz,a,ir,jc,ao,jao,iao,eps)
              ## subroutine triplet2csr(nrow,ncol,nnz,a,ir,jc,ao,jao,iao,eps)
              SIGNATURE = c( SS$signature, SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  "double", SS$signature, SS$signature,
                  "double"),
              
              nrow = nrow,
              ncol = ncol,
              nz = nz,
              
              x[[elenr]],
              ir,
              jc,
              
              entries = vector_dc( "double", nz),
              colindices = vector_dc( SS$type, nz),
              rowpointers = vector_dc( SS$type, nrow+1),
              
              eps,

              INTENT = c("r", "r", "rw",
                  "rw", "rw", "rw",
                  "rw", "rw", "rw",
                  "r"),
              NAOK=TRUE,
              PACKAGE = SS$package )

    

                                        #  print(z)
    if (z$nz == 0){
    ## if (identical(z$nz, 0)){
        ## print("special case")
        return(.newSpam(
            # rowpointers = c(1, rep_len64(2,nrow)),
            dimension = c(nrow, ncol),
            force64 = force64))
         ## return(new("spam",rowpointers=c(1L,rep.int(2L,nrow)), dimension=c(nrow,ncol)))
    }
   
    ## newx <- new("spam")
    ## slot(newx,"entries",check=FALSE) <- z$entries[1:z$nz]
    ## slot(newx,"colindices",check=FALSE) <- z$colindices[1:z$nz]
    ## slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
    ## slot(newx,"dimension",check=FALSE) <- c(nrow,ncol)
    ## return(newx)
    return(.newSpam(
        entries = z$entries[1:z$nz],
        colindices = z$colindices[1:z$nz],
        rowpointers = z$rowpointers,
        dimension = c(nrow,ncol),
        force64 = force64))
}

setMethod("as.spam", "list", as.spam.list) #  { function(x,eps) spam.list(x,eps=eps)})
setMethod("spam", "list", spam.list)
