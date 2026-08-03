#' Read csv file of sample data
#' 
#' This is a helper utility to read sample csv file into data frame.
#' The main use for external users is to read the option data
#' since it has a different format than other price timeseries data.
#'
#' @param symbol Character for the symbol of the time series. Default: dji
#' @param extdata_dir optionally specify user's own extdata folder
#'
#' @return The data.frame object
#'
#' @keywords data
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @importFrom utils read.csv
#'
#'
#' @examples
#' dji <- ldhmm.read_csv_by_symbol("dji")
#' spx <- ldhmm.read_csv_by_symbol("spx")
#'
### <======================================================================>
"ldhmm.read_csv_by_symbol" <- function(symbol = "dji", extdata_dir=NULL)
{
    c <- .data_config(symbol, extdata_dir=extdata_dir)
    if(nrow(c) != 1){
        stop(paste("Unknown symbol", symbol, "!\n"))
    }
    file <- paste(c$symbol, "_archive_", c$cols, ".csv", sep="")
    .read_file(file, extdata_dir=extdata_dir)
}
### <---------------------------------------------------------------------->
".locate_file" <- function(filename, extdata_dir=NULL)
{
    if(! (is.null(extdata_dir) || is.na(extdata_dir))) {
        f1 <- file.path(extdata_dir, filename)
        if(length(f1) > 0 & file.exists(f1)) return(f1)
    }
    # find the sample data location inside package
    f1 <- system.file("extdata", filename, package = "ldhmm")
    if(length(f1) > 0 & file.exists(f1)) return(f1)
    
    # during development, this is where it is!
    # f2 <- system.file("inst", "extdata", filename, package = "ldhmm")
    # if(length(f2) > 0 & file.exists(f2)) return(f2)
    return("")
}
### <---------------------------------------------------------------------->
".read_file" <- function(filename, extdata_dir=NULL)
{
        
    # regular file
    f <- .locate_file(filename, extdata_dir=extdata_dir)
    if (length(f) > 0 & file.exists(f)) {
        return(read.csv(f, header=TRUE, stringsAsFactors=FALSE))
    }
    
    # try zip file (arranged to save storage)
    f <- .locate_file(paste(filename, "zip", sep="."), extdata_dir=extdata_dir)
    if (length(f) > 0 & file.exists(f)) {
        df <- read.csv(unz(f, filename), header=TRUE, stringsAsFactors=FALSE)
        return(df)
    }
    stop(paste("Failed to locate file for", filename))
}
### <---------------------------------------------------------------------->

