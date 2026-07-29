### Working with returns #######################################################

##' @title Computing Returns or the Inverse Transformation
##' @param x matrix or vector of values to be turned into returns
##'        (if inverse = FALSE) or returns to be turned into the original data
##'        (if inverse = TRUE)
##' @param method method string; available are
##'        "logarithmic": log-returns (X_t = log(S_t/S_{t-1}))
##'        "simple":      simple returns (X_t = (S_t-S_{t-1})/S_{t-1})
##'        "diff":        differences (X_t = S_t-S_{t-1})
##' @param inverse logical indicating whether the inverse transformation
##'        (data from given returns) shall be computed (if TRUE, this
##'        requires 'start' to be specified)
##' @param start if inverse = TRUE, the last available value of the time
##'        series
##' @param start.date start date to be used if inverse = TRUE (currently only
##'        for 'xts' objects)
##' @return matrix (or vector) containing the returns or their 'inverses'
##' @author Marius Hofert
##' @note - For negative log-returns, use -returns(x) or
##'         returns(-x, inverse = TRUE, start = ...)
##'       - For percentage returns, use 100 * returns(x) or
##'         returns(x/100, inverse = TRUE, start = ...)
returns <- function(x, method = c("logarithmic", "simple", "diff"), inverse = FALSE,
                    start, start.date)
{
    ## Basics
    if(!is.matrix(x)) x <- cbind(x)
    is.zoo <- inherits(x, "zoo")
    method <- if(!missing(method)) {
                  match.arg(method, several.ok = TRUE) # would return all methods if 'method' wasn't provided
              } else {
                  method[1] # not provided => we set it to the default
              }
    len.method <- length(method)

    ## Main
    if(len.method == 1) {
        d <- ncol(x) # only grab here as this might be a recursive call
        switch(method,
               "logarithmic" = { # Logarithmic returns
                   if(inverse) {
                       ## Check whether 'start' (= S_0) has been provided
                       stopifnot(!missing(start), length(start) == d)
                       ## S_t = S_{t-1} * exp(X_t) = ... = S_{0} * exp(X_1 + X_2 + .. + X_t)
                       x.csum <- apply(x, 2, cumsum) # note: 'xts' lost here
                       x.csum <- rbind(rep(0, d), x.csum) # include 0 for S_0 * exp(0) = S_0 below
                       res <- rep(start, each = nrow(x.csum)) * exp(x.csum) # S_0, S_1, ..., S_t
                       res <- drop(res) # drops 1-column matrices to vectors
                       ## Return
                       sdate <- if(missing(start.date)) NA else start.date
                       if(is.zoo) as.xts(res, order.by = c(as.Date(sdate), index(x))) else res
                   } else {
                       ## X_t = log(S_t/S_{t-1})
                       res <- apply(x, 2, function(x.) diff(log(x.)))
                       res <- drop(res) # drops 1-column matrices to vectors
                       if(is.zoo) as.xts(res, order.by = index(x)[-1]) else res
                   }
               },
               "simple" = { # Simple returns
                   if(inverse) {
                       ## Check whether 'start' (= S_0) has been provided
                       stopifnot(!missing(start), length(start) == d)
                       ## S_t = S_{t-1} * (1 + X_t) = S_{t-2} * (1 + X_{t-1}) * (1 + X_t)
                       ##     = ... = S_0 * prod_{s = 1}^t (1 + X_s)
                       ##     (= exp(log(S_0) + sum_{s = 1}^t log1p(X_s)))
                       x.cprod <- apply(1 + x, 2, cumprod) # note: 'xts' lost here
                       x.cprod <- rbind(rep(1, d), x.cprod) # include 1 for S_0 * (1 + 0) = S_0 below
                       res <- rep(start, each = nrow(x.cprod)) * x.cprod # S_0, S_1, ..., S_t
                       res <- drop(res) # drops 1-column matrices to vectors
                       ## Return
                       sdate <- if(missing(start.date)) NA else start.date
                       if(is.zoo) as.xts(res, order.by = c(as.Date(sdate), index(x))) else res
                   } else {
                       ## X_t = (S_t-S_{t-1})/S_{t-1}
                       res <- apply(x, 2, function(x.) diff(x.)/head(x., n = -1))
                       res <- drop(res) # drops 1-column matrices to vectors
                       if(is.zoo) as.xts(res, order.by = index(x)[-1]) else res
                   }
               },
               "diff" = { # Differences
                   if(inverse) {
                       ## Check whether 'start' (= S_0) has been provided
                       stopifnot(!missing(start), length(start) == d)
                       ## S_t = S_{t-1} + X_t = S_{t-2} + X_{t-1} + X_t = ...
                       ##     = S_0 + sum_{s = 1}^t X_s
                       x.csum <- apply(x, 2, cumsum) # note: 'xts' lost here
                       x.csum <- rbind(rep(0, d), x.csum) # include 0 for S_0 + 0 = S_0 below
                       res <- rep(start, each = nrow(x.csum)) + x.csum # S_0, S_1, ..., S_t
                       res <- drop(res) # drops 1-column matrices to vectors
                       ## Return
                       sdate <- if(missing(start.date)) NA else start.date
                       if(is.zoo) as.xts(res, order.by = c(as.Date(sdate), index(x))) else res
                   } else {
                       ## X_t = S_t-S_{t-1}
                       res <- apply(x, 2, function(x.) diff(x.))
                       res <- drop(res) # drops 1-column matrices to vectors
                       if(is.zoo) as.xts(res, order.by = index(x)[-1]) else res
                   }
               },
               stop("Wrong 'method'"))
    } else { # length(method) > 1 => iteration over all margins separately
        if(len.method != ncol(x))
            stop("length(method) must be 1 or ncol(x)")
        res.lst <- lapply(seq_len(len.method), function(j) {
            if(inverse) {
                returns(x[,j], method = method[j], inverse = inverse,
                        start = as.numeric(start)[j], start.date = start.date)
            } else {
                returns(x[,j], method = method[j], inverse = inverse)
            }})
        res <- do.call(cbind, res.lst) # merge; note: x was rectangular before, so cbind() will work
        colnames(res) <- colnames(x) # ... since they were lost
        res
    }
}
returns_qrmtools <- returns # alias