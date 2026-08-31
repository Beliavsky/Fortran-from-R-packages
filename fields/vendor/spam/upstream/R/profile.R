# HEADER ####################################################
# This is file spam/R/profile.R.                            #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################


.format32 <- list(
    name      = "32-bit",
    type      = "integer", 
    signature = "integer",  
    package   = "spam")

.format64 <- list(
    name      = "64-bit",
    type      = "numeric",
    signature = "int64",
    package   = "spam64")


.format64 <- function(){
    if (!isNamespaceLoaded("spam64")) {
        stop("Large (64-bit) sparse matrices detected. Please load the required package 'spam64' and see the help page '?large_matrix'.")
    }
    list(
        name      = "64-bit",
        type      = "numeric",
        signature = "int64",
        package   = "spam64")
}

.format.spam <- function(x, ... , validate = getOption("spam.validate") ) {
    objects <- c(list(x), list(...))

    if (validate) for(o in objects)  stopifnot(validate_spam(o))
    
    for(o in objects){
        ## If both pointer vectors are of the same type,
        ## use this type to determine the format
        if(identical(typeof(o@colindices), typeof(o@rowpointers))) {
            if(identical(typeof(o@colindices), "double")){
                return(.format64())
            }
            next
        }
        
        ## As fallback use the length of the entries vector and the dimension
        if(nrow(o) > 2147483647 || ncol(o) > 2147483647 || 
           length(o@entries) > 2147483647){
            return(.format64())
        }
    }   
    return(.format32)
}


spam.Version <- function() {
  release <- utils::packageDescription("spam",field="Version")
  date <- utils::packageDescription("spam",field="Date")
  list(status="",
              major=sub("-","",substr(release,1,4)),
              minor=substr(sub("-","",substr(release,5,7)),1,1),
              year=substr(date,1,4),
              month=substr(sub("200.-","",date),1,2),
              day=sub("200.-..-","",date),
              version.string= paste("Spam version ",
                utils::packageDescription("spam",field="Version")," (",
                utils::packageDescription("spam",field="Date"),")",sep="")
              )
}
                
spam.version <- spam.Version()
class(spam.version) <- "simple.list"
  

"inefficiencywarning" <- function(msg,size) {
  maxsize <- if (is.logical(getOption("spam.inefficiencywarning"))) {
    ifelse(getOption("spam.inefficiencywarning"),1,Inf) } else { 
    getOption("spam.inefficiencywarning")
  }
  if (size>maxsize) warning(msg, call. = FALSE)
}
    
".onAttach" <- function (lib, pkg) {
   packageStartupMessage( spam.version$version.string," is loaded.",
       "\nType 'help( Spam)' or 'demo( spam)' for a short introduction ",
       "\nand overview of this package.",
       "\nHelp for individual functions is also obtained by ",
       "adding the\nsuffix '.spam' to the function name, e.g. 'help( chol.spam)'.")
}

.onLoad <- function(libname, pkgname) {
  
  default_options <- list(  spam.eps=.Machine$double.eps,   # smaller than this is considered as zero
     spam.force64=FALSE,
     spam.validate=FALSE,            # validate the spam object before calling a native routine for
                                     # increased stability.
     
     spam.drop=FALSE,                # drop passed to subset functions
     
     spam.printsize=100,             # the max size which we print as regular matrices
     spam.imagesize=10000,           # the max size which we display as regular matrices
     spam.cex=1200,                  # scaling factor for scatter displays

     spam.structurebased=TRUE,      # calculating on nonzero entries only...
     
     spam.inefficiencywarning=1e6,  # tell when something inefficient is done
     
     spam.trivalues=FALSE,           # with upper./lower/.tri return values (TRUE) or only structure?
     spam.listmethod="PE",           # method to be used when using spam.list

     spam.NAOK = FALSE,
     spam.safemodevalidity=TRUE,     # verify while S4 construction
     spam.dopivoting=TRUE,           # what type of back/forwardsolve?
     spam.cholsymmetrycheck=TRUE,     # Should symmetry be tested in the cholesky factorization
     spam.cholpivotcheck=TRUE,        # Should the pivot be tested?
     spam.cholupdatesingular="warning",     # ("error", "warning","NULL")
     spam.cholincreasefactor=c(1.25,1.25),
     spam.nearestdistincreasefactor=1.3,
     spam.nearestdistnnz=c(500^2,500)  )

	# Load the default settings:
	options(default_options)
}

powerboost <- function(flag = "on") {
    if (tolower(flag) %in% c("true","on","an","ein")) {
        options(spam.NAOK = TRUE,
                spam.safemodevalidity = FALSE,
                spam.cholsymmetrycheck = FALSE,
                spam.cholpivotcheck = FALSE,
                spam.eps = 1e-8)
    } else {
         options(spam.NAOK = FALSE,
                spam.safemodevalidity = TRUE,
                spam.cholsymmetrycheck = TRUE,
                spam.cholpivotcheck = TRUE,
                spam.eps = .Machine$double.eps)
    }
    invisible(NULL)
}
