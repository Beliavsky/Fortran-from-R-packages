### Tools for working with data sets ###########################################

##' @title Download data via quantmod's getSymbols()
##' @param x vector of, for example, ticker symbols (if src = "yahoo") or
##'        "EUR/USD" if (src = "oanda")
##' @param from start date as character string (e.g. 2015-01-01); if NULL,
##'        the earliest available date is picked
##' @param to end date; today unless otherwise specified
##' @param src source of the data
##' @param FUN function to apply to the downloaded data:
##'        - if data is NA (could not be retrieved): none
##'        - if provided: the given function
##'        - if not provided: Ad() if src = "yahoo"; Cl() if src = "google";
##'          none otherwise (e.g. if src = "oanda")
##' @param verbose logical indicating whether progress monitoring is done
##' @param warn logical indicating whether a warning is given showing the error
##'        message when fetching x fails
##' @param ... additional arguments passed to getSymbols()
##' @return (n, d)-matrix of data (an xts object)
##' @author Marius Hofert
##' @note - One could do...
##'         nenv <- new.env()
##'         getSymbols(c("^GSPC", "AAPL"), env = nenv, from = "2015-11-28", to = "2015-12-07")
##'         ... but we would not get progress monitoring then
##'       - In case of failure to get the data, note that we don't need to
##'         throw a warning, as getSymbols() produces them (unfortunately even
##'         if 'warnings = FALSE')
##'       - getSymbols calls underlying methods getSymbols.yahoo etc.
##'         Characters 'from' and 'to' are accepted by getSymbols.google,
##'         getSymbols.oanda and getSymbols.yahoo, but not by getSymbols.FRED
##'         => For FRED, we pick out the range separately (drawback: still the
##'            whole data is downloaded first -- so only use this hack for FRED
##'            (= Federal Reserve Economic Data) data)
get_data <- function(x, from = NULL, to = NULL,
                     src = c("yahoo", # "quandl", # no data available anymore for free
                             "oanda", "FRED", "google"),
                     FUN = NULL, verbose = TRUE, warn = TRUE, ...)
{
    ## Checking
    if(is.factor(x)) x <- as.character(x)
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.matrix(x)) x <- as.vector(x)
    stopifnot((d <- length(x)) >= 1)
    if(is.null(from)) from <- "1900-01-01" # to get all data available
    if(is.null(to)) to <- as.character(Sys.Date()) # to get all data available
    stopifnot(is.character(from) || inherits(from, "Date"),
              is.character(to) || inherits(to, "Date"),
              is.character(src), is.logical(verbose))
    src <- match.arg(src)

    ## Distinguish univariate/multivariate data
    if(d == 1) {

        ## Get the data
        start <- as.Date(from)
        end <- as.Date(to)
        time.diff <- difftime(end, start, units = "days")[[1]] # in days
        num.blocks.5y <- time.diff/(5*365) # number of blocks of 5y to fit in [start, end]
        if(src == "oanda" && num.blocks.5y > 1) { # get the data in blocks (5y max)
            periods <- seq(start, end, length.out = ceiling(num.blocks.5y)+1) # start,..., end (>= 2 dates of class Date)
            num.periods <- length(periods)-1 # >= 1 periods to get data from
            dat <- NULL
            for(k in seq_len(num.periods)) {
                from <- as.character(periods[k])
                to <- as.character(if(k == num.periods) periods[k+1] else periods[k+1]-1)
                if(verbose) cat("Getting", x, "from", from, "to", to, "\n") # progress output
                dat. <- get_data(x, from = from, to = to, src = src, FUN = FUN, verbose = verbose, ...) # recursion
                if(!(length(dat.) == 1 && is.na(dat.)))
                    dat <- rbind(dat, dat.) # exclude NAs (if no data available for that time period)
            }
            if(is.null(dat)) dat <- NA # as in case of an error -- also no data available
        } else { # if src != "oanda" or "oanda" but only one block, we can get it directly
            if(!is.character(from)) from <- as.character(from)
            if(!is.character(to)) to <- as.character(to)
            dat <-
                ## Archived (not working anymore) due to Quandl being archived
                ## if(src == "quandl") { # quandl
                ##     tryCatch(Quandl(x, type = "xts", start_date = from, end_date = to, ...),
                ##              error = function(e) e)
                ## } else { # !quandl
                    tryCatch(getSymbols(x, from = from, to = to, src = src, auto.assign = FALSE, ...),
                             error = function(e) e)
                ## }
            if(is(dat, "simpleError")) {
                if(warn) warning("Error when fetching ",x,": ",conditionMessage(dat)," (will use NA instead)")
                dat <- NA
            }
        }

        ## Applying FUN, naming
        if(identical(dat, NA)) {
            res <- NA
            names(res) <- x # use ticker symbol as name
        } else { # getting the data worked fine
            ## Apply FUN
            res <- if(is.null(FUN)) { # if not given, apply a useful default
                if(src == "yahoo") {
                    Ad(dat) # does a grep on colnames(x)
                } else if(src == "google") {
                    Cl(dat) # no Ad() available, use Cl()
                } else dat # nothing (important for src = "oanda")
            } else { # FUN given
                stopifnot(is.function(FUN))
                dat. <- FUN(dat)
                if(ncol(dat.) != 1)
                    stop("'FUN' has to return a single time series")
                dat.
            }
            if(ncol(res) == length(x)) {
                colnames(res) <- x # use ticker symbol as name
            } # otherwise, don't rename columns automatically
        }

    } else { # d > 1
        res <- vector("list", length = d)
        for(j in 1:d) {
            if(verbose) cat("Getting", x[j], "\n")
            res[[j]] <- get_data(x[j], from = from, to = to, src = src, FUN = FUN,
                                 verbose = verbose, ...) # recursion
        }
        res <- do.call(merge, res)
    }

    ## Return
    if(src == "FRED") { # "FRED" doesn't accept the above 'from' and 'to' => pick out range manually
        rows <- index(res)
        res[from <= rows & rows <= to, ]
    } else res
}
