### Plot of fitted GEV shape as a function of the block size ###################

##' @title Fitted GEV Shape as a Function of the Block Size
##' @param x numeric vector of data
##' @param blocksize numeric vector of block sizes for which to fit a GEV to
##'        the block maxima
##' @param estimate.cov logical indicating whether confidence intervals are
##'        computed
##' @param conf.level confidence level of the confidence intervals
##' @param CI.col color of the confidence interval region; can be NA
##' @param lines.args arguments passed to the underlying lines()
##' @param xlim see ?plot
##' @param ylim see ?plot
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param xlab2 label of the secondary x-axis
##' @param plot logical indicating whether a plot is done
##' @param ... additional arguments passed to the underlying plot()
##' @return invisibly returns the block sizes, the list of corresponding block
##'         maxima and the list of fitted GEV distributions.
##' @author Marius Hofert
GEV_shape_plot <- function(x, blocksize = tail(pretty(seq_len(length(x)/20), n = 64), -1),
                           estimate.cov = TRUE, conf.level = 0.95,
                           CI.col = adjustcolor(1, alpha.f = 0.2),
                           lines.args = list(), xlim = NULL, ylim = NULL,
                           xlab = "Block size", ylab = NULL,
                           xlab2 = "Number of blocks", plot = TRUE, ...)
{
    ## Checks
    stopifnot(length(blocksize) >= 2, is.logical(estimate.cov), 0 <= conf.level, conf.level <= 1,
              is.logical(plot))

    ## Compute block maxima for block size n (last one not filled)
    x <- as.numeric(x)
    block_maxima <- function(n) sapply(split(x, ceiling(seq_along(x)/n)), max)
    blockmax <- lapply(blocksize, function(n) block_maxima(n)) # block maxima data for each block size
    fits <- lapply(blockmax, fit_GEV_MLE, estimate.cov = estimate.cov) # fit a GEV for each block size
    ## Extract the fitted shape parameters and compute CIs
    xi <- sapply(fits, function(f) f$par[["shape"]])
    if(estimate.cov) {
        xi.SE <- sapply(fits, function(f) f$SE[["shape"]])
        q <- qnorm(1-(1-conf.level)/2)
        xi.CI.low <- xi - xi.SE * q
        xi.CI.up  <- xi + xi.SE * q
    }

    ## Plot
    if(plot) {
        if(is.null(xlim)) xlim <- range(blocksize)
        if(is.null(ylim))
            ylim <- if(is.na(CI.col) || !estimate.cov) range(xi) else range(xi, xi.CI.low, xi.CI.up)
        if(is.null(ylab)) ylab <- paste0("Estimated GEV shape parameter",
                                         if(estimate.cov) paste0(" with ", 100*conf.level,
                                                                 "% confidence intervals") else "")
        plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...) # plot region
        if(estimate.cov)
            polygon(x = c(blocksize, rev(blocksize)), y = c(xi.CI.low, rev(xi.CI.up)),
                    col = CI.col, border = NA) # CI
        do.call(lines, args = c(list(x = blocksize, y = xi), lines.args)) # actual plot
        pb <- axTicks(1) # used to be pretty(blocksize)
        axis(3, at = pb, labels = sapply(pb, function(n) ceiling(length(x)/n)))
        mtext(xlab2, side = 3, line = 3)
    }

    ## Return
    invisible(list(block.sizes = blocksize, block.maxima = blockmax, GEV.fits = fits))
}
