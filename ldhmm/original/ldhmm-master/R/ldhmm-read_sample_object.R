#' Read sample ldhmm object
#' 
#' This utility is used to read sample ldhmm object so that the user doesn't need to go through
#' lengthy optimization process to obtain a trained HMM for advanced features.
#'
#' @param symbol Character for the symbol of the time series. Default is \code{spx-daily-m10}
#' @param extdata_dir optionally specify user's own extdata folder
#'
#' @return The ldhmm object
#'
#' @keywords data
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#'
#' @examples
#' hs <- ldhmm.read_sample_object() # SPX daily 10-state HMM
#'
### <======================================================================>
"ldhmm.read_sample_object" <- function(symbol = "spx-daily-m10", extdata_dir=NULL)
{
    file <- paste(symbol, "p3.RData", sep="")
    f <- .ldhmm.locate_file(file, extdata_dir=extdata_dir)
    if (length(f) > 0 & file.exists(f)) {
        hs <- NULL
        load(file=f) # this loads the object hs
        return(hs)
    }
    stop(paste("Failed to locate RData file for", file))
}
### <---------------------------------------------------------------------->
".ldhmm.locate_file" <- function(filename, extdata_dir=NULL)
{
    if(! (is.null(extdata_dir) || is.na(extdata_dir))) {
        f1 <- file.path(extdata_dir, filename)
        if(length(f1) > 0 & file.exists(f1)) return(f1)
    }
    # find the sample data location inside package
    f1 <- system.file("extdata", filename, package = "ldhmm")
    if(length(f1) > 0 & file.exists(f1)) return(f1)
    
    # during development, this is where it is!
    f2 <- system.file("inst", "extdata", filename, package = "ldhmm")
    if(length(f2) > 0 & file.exists(f2)) return(f2)
    return("")
}
### <---------------------------------------------------------------------->
