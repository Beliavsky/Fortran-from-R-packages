### Hill estimator #############################################################

##' @title Hill Estimator
##' @param x vector of numeric data; those x corresponding to k[1]:k[2]
##'        must be > 0.
##' @param k vector of length 2, determining the smallest and largest number
##'        of order statistics of x to compute the Hill estimator for (the
##'        smallest needs to be >= 2). If of length 1, k is expanded
##'        by length(x).
##' @param conf.level confidence level for the asymptotic CIs
##' @return (k[1]:k[2], 5)-matrix with columns:
##'         - k = the k's used
##'         - k.prob = empirical probabilities corresponding to the k's
##'         - tail.index = estimated tail indices for the different k's
##'         - CI.low = lower CI for the tail indices
##'         - CI.up  = upper CI for the tail indices
##' @author Marius Hofert
Hill_estimator <- function(x, k = c(10, length(x)), conf.level = 0.95)
{
    ## Basics
    if(!is.vector(x)) x <- as.vector(x) # convert time series to vector (otherwise sort() fails)
    n <- length(x)
    if(length(k) == 1) k <- c(k, n)
    stopifnot(is.numeric(x), length(k) == 2, 2 <= k, diff(k) > 0,
              0 < conf.level, conf.level < 1)
    if(k[2] > n) # more meaningful error message
        stop("k[2] = ",k[2]," must be <= length(x) = ",n)

    ## Build ingredients
    k.min <- k[1]
    k.max <- k[2]
    k.ran <- k.min:k.max
    k.prob <- 1 - (k.ran-1) / n # k.ran == 2 corresponds to 1 - 1/n, k.ran = n to 1/n
    x. <- sort(x, decreasing = TRUE) # sorted in decreasing order
    if(x.[k.max] <= 0) # check smallest x.
        stop("sort(x) > 0 violated for some indices in range(k)")
    one.to.k.max <- seq_len(k.max)
    lx <- log(x.[one.to.k.max]) # fine (all arguments > 0)

    ## Hill estimator
    ## This form is frequently found in the literature (although
    ## EKM (1997, p. 190, Eq. (4.12)) use the one with k ~> k-1)
    lx.cmean.k.max <- cumsum(lx) / one.to.k.max # compute all cumulative means from 1 to k.max
    tail.index <- 1 / (lx.cmean.k.max[k.ran] - lx[k.ran])

    ## CI (see evir::hill or QRM::hillPlot)
    SE <- tail.index / sqrt(k.ran)
    q <- qnorm(1-(1-conf.level)/2)
    CI.low <- tail.index - q * SE
    CI.up  <- tail.index + q * SE

    ## Return
    cbind(k = k.ran, k.prob = k.prob, tail.index = tail.index,
          CI.low = CI.low, CI.up = CI.up)
}


### Hill plot ##################################################################

##' @title Hill Plot (Estimated Tail Index as a Function of the Number
##'        of Order Statistics)
##' @param x see ?Hill_estimator
##' @param k see ?Hill_estimator
##' @param conf.level see ?Hill_estimator
##' @param Hill.estimator object as returned by Hill_estimator()
##' @param log see ?plot
##' @param xlim see ?plot
##' @param ylim see ?plot
##' @param xlab see ?plot
##' @param ylab see ?plot
##' @param CI.col color of the confidence interval region; can be NA
##' @param lines.args arguments passed to the underlying lines()
##' @param xaxis2 logical indicating whether the 3rd axis is
##'        plotted
##' @param xlab2 label of the secondary x-axis
##' @param ... additional arguments passed to the underlying plot()
##' @return Hill plot (by side-effect)
##' @author Marius Hofert
Hill_plot <- function(x, k = c(10, length(x)), conf.level = 0.95, Hill.estimator = NULL,
                      log = "x", xlim = NULL, ylim = NULL,
                      xlab = "Order statistics", ylab = "Tail index",
                      CI.col = adjustcolor(1, alpha.f = 0.2), lines.args = list(),
                      xaxis2 = TRUE, xlab2 = "Empirical probability",
                      ...)
{
    ## Ingredients
    if(is.null(Hill.estimator))
        Hill.estimator <- Hill_estimator(x, k = k, conf.level = conf.level)
    k <- Hill.estimator[,"k"]
    k.prob <- Hill.estimator[,"k.prob"]
    tail.index <- Hill.estimator[,"tail.index"]
    CI.low <- Hill.estimator[,"CI.low"]
    CI.up  <- Hill.estimator[,"CI.up"]
    if(is.null(xlim)) xlim <- range(k)
    if(is.null(ylim))
        ylim <- if(is.na(CI.col)) range(tail.index) else range(tail.index, CI.low, CI.up)

    ## Hill plot with CI
    plot(NA, log = log, xlim = range(k), ylim = ylim, xlab = xlab, ylab = ylab, ...) # setup plot region
    polygon(x = c(k, rev(k)), y = c(CI.low, rev(CI.up)), col = CI.col, border = NA) # CI
    do.call(lines, args = c(list(x = k, y = tail.index), lines.args)) # Hill plot

    ## Third axis
    if(xaxis2) {
        ## Get tick locations of first axis
        at <- axTicks(1) # ticks in same places as x-axis

        ## Determine the corresponding probabilities
        ## Note: R may round order statistics down ('at' for the smallest:
        ## 0 instead of 1; 'at' for the largest: >> largest k)
        labs. <- numeric(length(at)) # same length
        labs.[at <= 0] <- 1 # probability level 100%
        labs.[at > max(k)] <- 0 # probability level 0%
        ii <- 0 < at & at <= max(k) # indices for which we have probabilities corresponding to k
        labs.[ii] <- k.prob[at[ii]]
        labs <- format(signif(labs., 4)) # round (fine)

        ## Plot third axis
        axis(3, at = at, labels = labs)
        mtext(xlab2, side = 3, line = 3)
    }
}
