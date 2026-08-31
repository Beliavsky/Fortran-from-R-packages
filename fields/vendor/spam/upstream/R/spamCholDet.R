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


# The main function here is very similar to chol.spam().
# We work with different files such that changes can be easily ported via `meld`.
# Both call z <- .spam.chol.basis(). Only, chol() has the verbose argument.

determinant.spam <- function(x, logarithm=TRUE, ...){

    
  argList <- list(...)
  if (!is.null(argList$memory)||!is.null(argList$memory)||!is.null(argList$pivot))
        warning("Arguments 'memory', 'eps' and 'pivot' will become defunct in the future. \nIf you really need these, pass via a Choleski factorization first.")
    
  z <- .chol.spam.basis(x, pivot="MMD", memory=list(), eps=getOption("spam.eps"), verbose=FALSE) 
    
  tmp <- 2* sum( log( z$lnz[ z$xlnz[ -(z$nrow+1)]]))

  logdet <- list()     
  if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)
  attr(logdet$modulus,"logarithm") <- logarithm

  logdet$sign <- 1
  attr(logdet,"class") <- "det"

  return(logdet)
}

determinant.spam.chol.NgPeyton <- function(x, logarithm = TRUE,...)
{
  tmp <- sum( log(x@entries[ x@rowpointers[-(x@dimension[1]+1)]]))
 
  logdet <- list()
  if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)
  attr(logdet$modulus,"logarithm") <- logarithm

  logdet$sign <- 1
  attr(logdet,"class") <- "det"

  return(logdet)
}

setMethod("determinant","spam",               determinant.spam)
setMethod("determinant","spam.chol.NgPeyton", determinant.spam.chol.NgPeyton)

## The ``Right Thing'' to do :
## base::det() calls [base::]determinant();
## our det() should call our determinant() :
det <- base::det
environment(det) <- environment()## == asNamespace("Matrix")
######################################################################
########################################################################

