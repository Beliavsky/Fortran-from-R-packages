
stdEr.maxLik <- function(x, eigentol=1e-12, ...) {
   ## if(!inherits(x, "maxLik"))
   ##    stop("'stdEr.maxLik' called on a non-'maxLik' object")
   ## Here we should actually coerce the object to a 'maxLik' object, dropping all the subclasses...
   ## Instead, we force the program to use maxLik-related methods
   if(!is.null(vc <- vcov(x, eigentol=eigentol))) {
      s <- sqrt(diag(vc))
      names(s) <- names(coef(x))
      return(s)
   }
   # if vcov is not working, return NULL
   return(NULL)
}

