#' Utility to download time series from FRED
#' 
#' This utility downloads time series from FRED. It serves as a data source for daily data,
#' e.g. SP500 for S&P 500, and VIXCLS for CBOE VIX index. This can be concatenated to the 
#' static data to provide daily updates.
#'
#' @param symbol character, the name of the time series
#' @param col_out character, the name of the output closing price column. Default: "Close"
#' @param do.logr logical, if \code{TRUE} (default), produce xts object of logr; otherwise, just the \code{col_out} column.
#'                Be aware that, because logr uses diff, the first day close will be deleted.
#'
#' @return The xts object for the time series
#'
#' @keywords data
#'
#' @export
#'
#' @importFrom utils download.file
#'
#' @examples
#' \dontrun{
#' ldhmm.fred_data("VIXCLS")
#' }
### <======================================================================>
ldhmm.fred_data <- function (symbol, col_out="Close", do.logr=TRUE)
{

    tmp <- tempfile()
    FRED.URL <- "https://fred.stlouisfed.org/series"
    URL <- paste(FRED.URL, "/", symbol, "/downloaddata/", symbol, ".csv", sep="")
    utils::download.file(URL, destfile=tmp, quiet=TRUE)
    fr <- read.csv(tmp, na.strings=".")
    unlink(tmp)

    ts <- xts(as.matrix(fr[,-1]),
              as.Date(fr[,1],origin='1970-01-01'),
              src='FRED', symbol=symbol,
              updated=Sys.time())
    
    dim(ts) <- c(NROW(ts),1)
    colnames(ts) <- col_out
    ts <- na.exclude(ts)
    if (!do.logr) return (ts)
    
    # derive log returns
    ts$logr <- suppressWarnings(diff(log(ts)))
    ts <- ts[!is.na(ts$logr)] # remove NA
    return (ts)
}
### <---------------------------------------------------------------------->
    
