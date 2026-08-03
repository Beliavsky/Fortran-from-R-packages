#' Read sample data
#' 
#' Read sample data by specifying the symbol. The two utilities, \code{ecd.data} and \code{ecd.data.arr},
#' serves for slightly different purpose. 
#' \code{ecd.data} works off the xts object that has two rows: 
#' the prices and log-returns indexed by the dates. 
#' \code{ecd.data.arr} and \code{ecd.data.ts} separate the data into list of three vectors: x is the log-return, p is the prices, and d is the dates.
#' And allows for more sophisticated call for range of dates, and different ways of slice and lag.
#' \code{ecd.data.arr} takes symbol as input, while \code{ecd.data.ts} takes an xts object.
#' 
#' @param symbol character, the symbol of the time series. Default: dji
#' @param ts xts, the time series
#' @param start.date,end.date Date or character of ISO format (YYYY-MM-DD), to specify the date range, 
#'                            default is from 1950-01-01 to 2015-12-31. 
#'                            Set start.date and end.date to NULL or "" if you wish to get the entire time series.
#' @param on character, specify the calendar interval, days, weeks, months. Default is \code{days}.
#' @param lag integer, specify the lags of return calculation, default is 1.
#' @param drop integer, specify number of largest outliners to drop, default is 0.
#' @param repeated logical, specify whether to use repeated sampling or unique sampling, default is \code{TRUE}.
#'                 Using "repeated" sampling can reduce noise due to insufficient sample size. This is particularly useful for larger lags.
#' @param do.kurtosis logical, if specified, calculate mean, sd, var, skewness, and kurtosis, default is \code{FALSE}.
#' @param cache logical, use R's options memory to cache xts data, default is \code{TRUE}.
#' 
#' @return \code{ecd.data} returns an xts object for the time series, with two columns - "Close" and "logr".
#'         \code{ecd.data.arr} and \code{ecd.data.ts} return a list of three vectors: x is the log-return, p is the prices, and d is the dates.
#'
#' @keywords timeseries xts sample-data
#'
#' @importFrom utils head
#' @importFrom xts endpoints
#' @importFrom stats na.exclude
#'
#' @export ecd.data
#' @export ecd.data.arr
#' @export ecd.data.ts
#'
#' @examples
#' dji <- ecd.data()
#' wti <- ecd.data("wti")
#' spx <- ecd.data.arr("spx", lag=5)
### <======================================================================>
"ecd.data" <- function(symbol = "dji")
{
    if(grepl("option", symbol)) {
        stop(paste("Option is not supported for symbol", symbol))
    }
    
    # basic settings
    
    dt <- "Date"
    col_out <- "Close"
    
    c <- .ecd.data_config(symbol)
    if(nrow(c) != 1){
        stop(paste("Unknown symbol", symbol, "for sample data!"))
    }
    
    df <- ecd.read_csv_by_symbol(symbol)
    
    ts <- ecd.df2ts(df, date_format = c$date_format,
                    dt = dt, col_in = c$col_in, col_out = col_out)
    xtsAttributes(ts) <- list(symbol = symbol)
    ts
}
### <---------------------------------------------------------------------->
#' @rdname ecd.data
ecd.data.arr <- function(symbol="dji", 
                         start.date="1950-01-01", end.date="2015-12-31", 
                         on="days", lag=1, drop=0, repeated=TRUE,
                         cache=TRUE, do.kurtosis=FALSE) {

    # get data
    force <- (!cache) # force to read from file
    if (is.null(getOption("ecd.data.symbol"))) force <- TRUE
    else if (getOption("ecd.data.symbol") != symbol) force <- TRUE
    
    ts <- NULL
    if (!force) ts <- getOption("ecd.data.ts")
    else {
        ts <- ecd.data(symbol)
        options("ecd.data.symbol"=symbol)
        options("ecd.data.ts"=ts)
    }

    ecd.data.ts(ts, start.date=start.date, end.date=end.date, 
                on=on, lag=lag, drop=drop, repeated=repeated,
                do.kurtosis=do.kurtosis)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.data
ecd.data.ts <- function(ts, start.date="1950-01-01", end.date="2015-12-31", 
                        on="days", lag=1, drop=0, repeated=TRUE,
                        do.kurtosis=FALSE) {
    
    # date range
    if (is.null(start.date)) start.date <- "" 
    if (is.null(end.date)) end.date <- ""
    if (start.date=="") start.date <- min(index(ts))
    if (end.date=="") end.date <- max(index(ts))
    
    if (is(start.date, "Date")) start.date <- as.character(start.date)
    if (!is(start.date, "character")) stop("start.date must be in ISO-date string or Date")
    if (is(end.date, "Date")) end.date <- as.character(end.date)
    if (!is(end.date, "character")) stop("end.date must be in ISO-date string or Date")
    
    ts1 <- ts[paste(start.date, "/", end.date, sep="")]
    ep <- xts::endpoints(ts1, on=on, k=1)
    ts2 <- ts1[ep[1:length(ep)]]
    ts2$logr <- log(ts2$Close) - stats::lag(log(ts2$Close))
    ts2 <- stats::na.exclude(ts2)
    
    d <- zoo::index(ts2)
    x <- as.numeric(ts2$logr)
    p <- as.numeric(ts2$Close)
    
    if (lag > 1) {
        NS <- length(x)
        I.from <- 1:(NS-lag)
        I.to <- (1+lag):NS
        x <- log(p[I.to]/p[I.from])
        d <- d[I.to]
        p <- p[I.to]
        
        if (! repeated) {
            U <- seq(1, length(I.to), by=lag)
            x <- x[U]
            d <- d[U]
            p <- p[U]
        }
    }
    if (drop > 0) {
        `%notin%` <- function (x, table) match(x, table, nomatch = 0L) == 0L
        z <- abs(x)
        z1 <- rev(z[order(z)])
        J <- which(abs(x) %notin% utils::head(z1, drop))
        x <- x[J]
        d <- d[J]
        p <- p[J]
    }
    rs <- list(d=d, x=x, p=p)
    
    if (do.kurtosis) {
        rs$mean <- mean(x, na.rm=TRUE)
        rs$sd <- sd(x, na.rm=TRUE)
        rs$var <- var(x, na.rm=TRUE)
        rs$skewness <- moments::skewness(x, na.rm=TRUE)
        rs$kurtosis <- moments::kurtosis(x, na.rm=TRUE)
    }
    return(rs)
}
### <---------------------------------------------------------------------->
