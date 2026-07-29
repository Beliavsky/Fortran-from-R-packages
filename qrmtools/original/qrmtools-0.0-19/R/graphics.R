### General graphical tools ####################################################

##' @title Image Indicating NAs in a Data Set
##' @param x matrix (ideally an xts object)
##' @param col colors for NA and non-NA, respectively
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param text see mtext()
##' @param side see mtext()
##' @param line see mtext()
##' @param adj see mtext()
##' @param ... additional arguments passed to image()
##' @return invisible()
##' @author Marius Hofert
NA_plot <- function(x, col = c("black", "white"), xlab = "Time", ylab = "Component",
                    text = "Black: NA; White: Available data", side = 4, line = 1, adj = 0,
                    ...)
{
    stopifnot(is.matrix(x))
    x. <- if(inherits(x, "xts")) {
        index(x) # use the time points
    } else {
        rn <- rownames(x)
        if(is.null(rn)) seq_len(nrow(x)) else rn # if available, use row names, otherwise numbers
    }
    image(x = x., y = seq_len(ncol(x)), z = is.na(x),
          col = rev(col), xlab = xlab, ylab = ylab, ...)
    if(inherits(text, "call") || nchar(text) > 0)
        mtext(text, side = side, line = line, adj = adj)
    invisible()
}

##' @title Plot of a Matrix
##' @param x matrix
##' @param ran range over which to plot the colors (can be used to force (-1, 1),
##'        for example)
##' @param ylim y-axis limits in reverse order (for the rows to appear 'top down')
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param scales see levelplot(); if NULL, labels and ticks are omitted
##' @param at see levelplot()
##' @param colorkey see levelplot()
##' @param col default colors for the color key
##' @param col.regions see levelplot()
##' @param ... additional arguments passed to levelplot()
##' @return the plot, a Trellis object
##' @author Marius Hofert
##' @note - Another option would be:
##'         corrplot::corrplot(err, method = "color", col = grey(seq(0.4, 1, length.out = 200)),
##'                            tl.col = "black", is.corr = FALSE)
##'       - Check via examples on ?matrix_plot
matrix_plot <- function(x, ran = range(x, na.rm = TRUE), ylim = rev(c(0.5, nrow(x) + 0.5)),
                        xlab = "Column", ylab = "Row",
                        scales = list(alternating = c(1,1), tck = c(1,0),
                                      x = list(at = pretty(1:ncol(x)), rot = 90),
                                      y = list(at = pretty(1:nrow(x)))),
                        at = NULL, colorkey = NULL, col = c("royalblue3", "white", "maroon3"),
                        col.regions = NULL, ...)
{
    stopifnot(is.matrix(x), (nr <- nrow(x)) >= 1,
              is.null(scales) || is.list(scales), is.numeric(at) || is.null(at),
              is.list(colorkey) || is.null(colorkey))
    ## Determine colors for the color key
    if(all(ran >= 0)) { # all non-NA entries >= 0
        if(is.null(at)) at <- seq(0, ran[2], length.out = 200)
        if(is.null(col.regions))
           col.regions <- colorRampPalette(c(col[2], col[3]), space = "Lab")(200)
    } else if(all(ran <= 0)) { # all non-NA entries <= 0
        lcol <- length(col)
        stopifnot(2 <= lcol, lcol <= 3)
        if(is.null(at)) at <- seq(ran[1], 0, length.out = 200)
        if(is.null(col.regions))
           col.regions <- colorRampPalette(c(col[1], col[2]), space = "Lab")(200)
    } else { # entries < 0 && entries > 0
        lcol <- length(col)
        stopifnot(2 <= lcol, lcol <= 3)
        if(lcol == 2) col <- c(0, col)
        ## Scale so that 0 gets the 'middle' color col[2]
        frac1 <- -ran[1]/diff(ran)
        frac2 <- ran[2]/diff(ran) # => frac1 + frac2 = 1 (and both in [0,1])
        if(is.null(at)) at <- seq(ran[1], ran[2], length.out = 200)
        if(is.null(col.regions))
            col.regions <- c(colorRampPalette(col[1:2], space = "Lab")(floor(200*frac1)),
                             colorRampPalette(col[2:3], space = "Lab")(ceiling(200*frac2)))
    }
    if(min(x, na.rm = TRUE) < at[1] || max(x, na.rm = TRUE) > at[length(at)])
        stop("'x' values outside the range spanned by 'at'. Choose 'at' appropriately.")
    ## Plot
    ## Note: t() does not affect the above as we only determine ranges/max/min there
    levelplot(t(x), ylim = ylim, # correct (x-axis = column; y-axis = row (and decreasing))
              xlab = xlab, ylab = ylab, col.regions = col.regions,
              scales = if(is.null(scales)) list(alternating = c(0,0), tck = c(0,0)) else scales,
              at = at, colorkey = if(is.null(colorkey)) list(at = at) else colorkey, ...)
}

##' @title Density Plot of the Values from a Lower Triangular Matrix
##' @param x matrix
##' @param xlab x-axis label; see plot()
##' @param main title; see plot()
##' @param text see mtext()
##' @param side see mtext()
##' @param line see mtext()
##' @param adj see mtext()
##' @param ... additional arguments passed to plot()
##' @return invisible()
##' @author Marius Hofert
matrix_density_plot <- function(x, xlab = "Entries in the lower triangular matrix",
                                main = "", text = NULL, side = 4, line = 1, adj = 0, ...)
{
    if(!is.matrix(x)) x <- rbind(x, deparse.level = 0L)
    d <- ncol(x)
    x.vec <- x[lower.tri(x)] # grab out values from the lower triangular matrix
    dens.x <- density(x.vec) # density
    plot(dens.x, main = main, xlab = xlab, ...) # plot
    rug(x.vec, col = "black") # rugs
    if(is.null(text))
        text <- substitute("Dimension: "*d.*"; Sample size: "*n.*"; Bandwidth: "*b.,
                           list(d. = d, n. = dens.x$n, b. = format(dens.x$bw, digits = 3, scientific = TRUE)))
    if(inherits(text, "call") || nchar(text) > 0)
        mtext(text, side = side, line = line, adj = adj)
    invisible()
}

##' @title P-P Plot
##' @param x data (a vector of convertible to such)
##' @param FUN hypothesized *distribution* function
##' @param pch plot symbol
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
##' @note as.vector() required since sort(x) == x for an 'xts' object!
pp_plot <-  function(x, FUN = pnorm, pch = 20,
                     xlab = "Theoretical probabilities", ylab = "Sample probabilities",
                     ...)
{
    p <- ppoints(length(x)) # theoretical probabilities of sorted data
    y <- FUN(sort(as.vector(x))) # hypothesized quantiles of sorted data
    plot(p, y, pch = pch, xlab = xlab, ylab = ylab, ...)
    abline(a = 0, b = 1)
    invisible()
}

##' @title Q-Q Plot
##' @param x data (a vector or convertible to such)
##' @param FUN hypothesized *quantile* function (vectorized)
##' @param method method used to construct the Q-Q line:
##'        "theoretical": theoretically true line which helps
##'                       deciding whether x comes from exactly
##'                       the distribution specified by FUN
##'        "empirical": qqline() which interpolates
##'                     (F^-(0.25), F_n^-(0.25)) and (F^-(0.75)), F_n^-(0.75))
##'                     and which helps deciding whether x comes from
##'                     a location-scale transformed distribution
##'                     specified by FUN.
##'        Example:
##'        z <- rnorm(1000)
##'        qq_plot(z, qnorm) # fine
##'        mu <- 3
##'        sig <- 2
##'        z. <- mu+sig*z
##'        qq_plot(z., qnorm) # not fine
##'        qq_plot((z.-mean(z.))/sd(z.), qnorm) # fine again
##'        qq_plot(z., qnorm, method = "empirical") # fine again
##' @param pch plot symbol
##' @param do.qqline logical indicating whether a Q-Q line is plotted
##' @param qqline.args list containing additional arguments passed to the
##'        underlying abline() functions
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
##' @note - as.vector() required since sort(x) == x for an 'xts' object!
##'       - This is a convenience function which does the same as
##'         qqplot(FUN(ppoints(length(x))), x)
##'         qqline(x, distribution = function(p) FUN(p))
##'         ... but has nicer default labels, does the Q-Q line by default
##'         and uses the *theoretical* Q-Q line by default, not an
##'         empirically estimated one via
##'         x <- FUN(probs)
##'         y <- quantile(x, probs = c(0.25, 0.75))
##'         slope <- diff(y) / diff(x)
##'         int <- y[1] - slope * x[1]
##'         abline(a = int, b = slope)
qq_plot <-  function(x, FUN = qnorm, method = c("theoretical", "empirical"),
                     pch = 20, do.qqline = TRUE, qqline.args = NULL,
                     xlab = "Theoretical quantiles", ylab = "Sample quantiles",
                     ...)
{
    x. <- x[!is.na(x)] # grab out non-NA data
    q <- FUN(ppoints(length(x.))) # theoretical quantiles of sorted data
    y <- sort(as.vector(x.)) # compute order statistics (sample quantiles)
    plot(q, y, pch = pch, xlab = xlab, ylab = ylab, ...)
    if(do.qqline) {
        method <- match.arg(method)
        switch(method,
        "theoretical" = { # intercept 0, slope 1
            if(is.null(qqline.args)) qqline.args <- list(a = 0, b = 1) # defaults
            do.call(abline, qqline.args)
        },
        "empirical" = { # estimates intercept and slope
            if(is.null(qqline.args)) qqline.args <- list()
            arg.lst <- c(list(x., distribution = FUN), qqline.args)
            do.call(qqline, arg.lst)
        },
        stop("Wrong 'method'"))
    }
    invisible()
}

##' @title Plot of a Step Function
##' @param x numeric vector of x-values
##' @param y numeric vector of corresponding y-values
##' @param y0 y-value of the graph extending to the left of the first x-value
##' @param x0 smallest x-value
##' @param x1 largest x-value
##' @param method character string indicating the method to be used
##' @param log character indicating whether to plot in log-scale
##' @param verticals logical indicating whether vertical lines are to be plotted
##' @param do.points logical indicating whether points are to be plotted
##' @param add logical indicating whether the current plot is added to the last one
##' @param col color
##' @param main title
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param plot.args list with additional arguments passed to the underlying plot()
##' @param segments.args list with additional arguments passed to the underlying segments()
##' @param points.args list with additional arguments passed to the underlying points()
##' @return nothing (plot by side-effect)
##' @author Marius Hofert
step_plot <- function(x, y, y0 = NA, x0 = NA, x1 = NA, method = c("edf", "eqf"), log = "",
                      verticals = NA, do.points = NA, add = FALSE,
                      col = par("col"), main = "", xlab = "x", ylab = "Function value at x",
                      plot.args = NULL, segments.args = NULL, points.args = NULL)
{
    ## Check x and y and sort (as plot.stepfun())
    nx <- length(x)
    ny <- length(y)
    if(nx != ny)
        stop("'x' and 'y' need to be numeric vectors of the same length")
    n <- nx
    ord <- order(as.numeric(x))
    x. <- x[ord]
    y. <- y[ord] # order as x to not mess up function
    if(!is.na(x0) && x0 > x.[1])
        stop("'x0' must be less than or equal to the smallest x")
    if(!is.na(x1) && x1 < x.[nx])
        stop("'x1' must be greater than or equal to the largest x")

    ## Verticals and points
    stopifnot(is.logical(verticals), is.logical(do.points))
    if(is.na(verticals)) verticals <- n >= 100
    if(is.na(do.points)) do.points <- n <  100
    method <- match.arg(method)
    if(method == "edf") {
        if((verticals && is.na(y0)) || (!is.na(x0) && is.na(y0)))
            stop("'y0' needs to be provided") # ... as we otherwise cannot plot the first vertical line
    }
    if(!is.na(y0) && y0 > min(y.))
        stop("'y0' must be less than or equal to the smallest y")

    ## Determine x values
    x.. <- c(x0, x., x1) # extend (possibly by NA)
    xlim <- range(x.., na.rm = TRUE)
    n.. <- length(x..) # update

    ## Determine y values
    if(is.na(y0)) y0 <- y.[1] # default smallest y = smallest y provided
    y0. <- min(y0, y.[1]) # the smaller of the two
    y.. <- switch(method,
                  "edf" = {
                      c(y0., y., y.[n]) # extend (possibly by NA); note that 'y1' is always y.[n]
                  },
                  "eqf" = {
                      c(y0., y.[-1], y.[n], y.[n]) # need to shift all horizontal and vertical segments upwards one 'step'
                  },
                  stop("Wrong 'method'"))
    ylim <- range(y.., na.rm = TRUE)

    ## Plot
    if(!add) {
        nms <- if(is.null(plot.args)) "" else names(plot.args)
        if(grepl("xlim", nms) || grepl("ylim", nms))
            stop("'xlim' and 'ylim' are determined by step_plot()")
        plt.args <- c(list(x = NA, log = log, xlim = xlim, ylim = ylim,
                           col = col, main = main, xlab = xlab, ylab = ylab),
                      plot.args)
        do.call(plot, args = plt.args) # plot region
    }
    sgmts.args <- c(list(x0 = x..[-n..], y0 = y..[-n..], x1 = x..[-1], y1 = y..[-n..], col = col),
                    segments.args)
    do.call(segments, args = sgmts.args) # horizontal segments
    if(verticals) {
        sgmts.args. <- c(list(x0 = x..[-1], y0 = y..[-n..], x1 = x..[-1], y1 = y..[-1], col = col),
                         segments.args)
        do.call(segments, args = sgmts.args.) # vertical segments
    }
    if(do.points) {
        nms <- if(is.null(points.args)) "" else names(points.args)
        if(!grepl("pch", nms)) points.args <- c(list(pch = 20), points.args)
        pnts.args <- c(list(x = x., y = y., col = col), points.args)
        do.call(points, args = pnts.args) # only put points where original data was, not extended one
    }
}

##' @title Plot of (an) Empirical Distribution Function(s)
##' @param x vector or list of vectors; if a list, each element corresponds
##'        to the x-values of one empirical distribution function
##' @param y0 see step_plot()
##' @param x0 see step_plot()
##' @param x1 see step_plot()
##' @param log see step_plot()
##' @param verticals see step_plot()
##' @param do.points see step_plot()
##' @param col color or vector of such
##' @param main see step_plot()
##' @param xlab see step_plot()
##' @param ylab see step_plot()
##' @param ... additional arguments passed to step_plot()
##' @return nothing (plot by side-effect)
##' @author Marius Hofert
edf_plot <- function(x, y0 = 0, x0 = NA, x1 = NA, log = "",
                     verticals = NA, do.points = NA, col = par("col"),
                     main = "", xlab = "x", ylab = "Distribution function at x", ...)
{
    if(!is.list(x)) x <- list(x)
    len <- length(x)
    stopifnot(len >= 1)
    xran <- range(unlist(x))
    if(is.na(x0)) x0 <- xran[1]
    if(is.na(x1)) x1 <- xran[2]
    if(length(col) != len) col <- rep(col, length.out = len)

    ## Loop over all plots
    for(l in 1:len) {
        x. <- sort(as.numeric(x[[l]]))
        y. <- ecdf(x.)(x.)
        step_plot(x., y., y0 = y0, x0 = x0, x1 = x1, log = log,
                  verticals = verticals, do.points = do.points,
                  add = l >= 2, col = col[l], main = main, xlab = xlab, ylab = ylab, ...)
    }
}

##' @title Plot of (an) Empirical Quantile Function(s)
##' @param x vector or list of vectors; if a list, each element corresponds
##'        to the x-values of one empirical quantile function
##' @param y0 see step_plot()
##' @param x0 see step_plot()
##' @param x1 see step_plot()
##' @param log see step_plot()
##' @param verticals see step_plot()
##' @param do.points see step_plot()
##' @param col color or vector of such
##' @param main see step_plot()
##' @param xlab see step_plot()
##' @param ylab see step_plot()
##' @param ... additional arguments passed to step_plot()
##' @return nothing (plot by side-effect)
##' @author Marius Hofert
eqf_plot <- function(x, y0 = NA, x0 = 0, x1 = 1, log = "",
                     verticals = NA, do.points = NA, col = par("col"),
                     main = "", xlab = "x", ylab = "Quantile function at x", ...)
{
    if(!is.list(x)) x <- list(x)
    len <- length(x)
    stopifnot(len >= 1)
    if(length(col) != len) col <- rep(col, length.out = len)

    ## Loop over all plots
    for(l in 1:len) {
        y. <- sort(as.numeric(x[[l]]))
        x. <- ecdf(y.)(y.)
        step_plot(x., y., y0 = y0, # default: NA (=> quantile function is the smallest y on (0, 1/n])
                  x0 = x0, x1 = x1, method = "eqf", log = log,
                  verticals = verticals, do.points = do.points,
                  add = l >= 2, col = col[l], main = main, xlab = xlab, ylab = ylab, ...)
    }
}
