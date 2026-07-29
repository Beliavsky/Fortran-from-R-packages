### Mean excess function and plot ##############################################

##' @title (Sample) Mean Excess Function
##' @param x numeric vector of data.
##' @param omit number of unique last observations to be omitted
##'        from the sorted data (as mean excess estimator becomes unreliable for
##'        these observations as thresholds)
##' @return Two-column matrix giving the sorted data without the omit-largest
##'         unique values (first column) and the corresponding values of the
##'         sample mean excess function (second column).
##' @author Marius Hofert
##' @note - Main idea:
##'         1) e_n(u) = (sum_{i=1}^n (X_{(i)}-u) I(X_{(i)} > u)) / sum_{i=1}^n I(X_{(i)} > u), u < X_{(n)}
##'         2) e_n(X_{(k)}) = sum_{i=k+1}^n (X_{(i)}-X_{(k)}) / (n-k), k in {1,..,n-1}
##'                         = ( (1/(n-k)) sum_{i=k+1}^n X_{(i)} ) - X_{(k)}
##'            ... but this holds only if there are *no* ties.
##'         3) Make the data unique, compute 2), then assign each duplicated value
##'            the corresponding value of the mean excess function.
mean_excess_np <- function(x, omit = 3)
{
    ## Check
    stopifnot(omit >= 1) # last one has to be omitted (otherwise e_n(X_{(i)}) not defined)

    ## Sort data and make it unique
    x <- sort(as.numeric(x)) # sorted data
    runs <- rle(x) # runs (make sense here because of x is sorted)
    xu <- runs$values # unique x values (sorted)
    nu <- length(xu) # how many unique x values
    ## Omit the last so-many unique values
    nuo <- nu - omit # how many we consider
    if(nuo < 1)
        stop("Not enough unique x values for computing the sample mean excess function.")
    ## => nuo >= 1 => nu >= 2.

    ## Compute e_n(X_{k}) for k = 1,..,n based on unique X values and
    ## n = nu as described above
    ## Note: This works (and thus shows the method) but is slow:
    ##       enu <- vapply(1:nuo, function(k) mean(xu[(k+1):nu]-xu[k]), NA_real_)
    ## => Idea: compute partial sums first and then do 'look-up'
    csx <- rev(cumsum(xu[nu:2])) # sum(xu[2:nu]), sum(xu[3:nu]), ..., xu[nu-1] + xu[nu], xu[nu]
    enu <- vapply(1:nuo, function(k) (csx[k]/(nu-k)) - xu[k], NA_real_) # note: nuo = nu - omit <= nu - 1 => nuo + 1 <= nu
    ## => Checked (coincides with the above computation of 'enu')

    ## Expand the 'enu' values (to plot all points according to their correct number of appearances)
    xu.times <- runs$lengths[1:nuo] # numbers of how often the first nuo-many unique values appear
    en <- rep(enu, times = xu.times) # expand (corresponding mean excesses e_n(X_{(i)}))

    ## Return
    cbind(x = x[1:sum(xu.times)], mean.excess = en) # X_{(i)}'s up to the last unique 'omit'-many data points
}

##' @title Mean Excess Function of a GPD
##' @param x evaluation points (thresholds)
##' @param shape GPD shape parameter xi
##' @param scale GPD scale parameter beta
##' @return vector of length as x, giving the mean excess function of the
##'         specified GPD (or NA if scale + shape * x <= 0).
##' @author Marius Hofert
##' @note Call with 'x - <chosen threshold>' to compare against mean excess plot
##'       as a function of the (general) thresholds.
mean_excess_GPD <- function(x, shape, scale)
{
    stopifnot(shape < 1)
    num <- scale + shape * x
    res <- rep(NA, length(x))
    res[num > 0] <- num / (1 - shape)
    res
}

##' @title (Sample) Mean Excess Plot
##' @param x numeric vector of data.
##' @param omit number of unique last observations to be omitted
##'        from the sorted data (as mean excess plot becomes unreliable for
##'        these observations as thresholds)
##' @param xlab x-axis label; see plot()
##' @param ylab y-axis label; see plot()
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
##' @note - Note that QRM::MEplot() only considers *unique* values. In particular,
##'         each value in the plot then appears equally often (which is actually
##'         wrong when using alpha transparency). This is also visible for the
##'         Danish data: str(MEplot(danish, omit = 3)) => 1645 points and not 2164
##'         as we obtain.
##'       - evd::mrlplot(QRM::danish) essentially uses type = "l" with pointwise
##'         asymptotic CIs, but does not seem too useful. qqtest's idea for CIs
##'         would be good but then probably too slow.
##'       - We could (and did once) pass 'shape', 'scale', 'q', 'length.out',
##'         'lines.args', 'log' (see tail_plot()) and add a fitted GPD mean
##'         excess function this way.
##'         => Not worth it, also doesn't look too good (Danish fire) and at
##'            that point one typically doesn't have a fitted GPD yet anyways.
##'       - The mean excess function of the fitted GPD would be:
##'         mean_excess_GPD(q-<chosen threshold>, shape = shape, scale = scale)
mean_excess_plot <- function(x, omit = 3,
                             xlab = "Threshold", ylab = "Mean excess over threshold", ...)
{
    plot(mean_excess_np(x, omit = omit), xlab = xlab, ylab = ylab, ...)
    invisible()
}


### Plot of fitted GPD shape as a function of the threshold ####################

##' @title Fitted GPD Shape as a Function of the Threshold
##' @param x numeric vector of data
##' @param thresholds numeric vector of thresholds for which to fit a GPD to
##'        the excesses
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
##' @return invisibly returns the thresholds, the list of corresponding threshold
##'         excesses and the list of fitted GPDs.
##' @author Marius Hofert
GPD_shape_plot <- function(x, thresholds = seq(quantile(x, 0.5), quantile(x, 0.99), length.out = 65),
                           estimate.cov = TRUE, conf.level = 0.95,
                           CI.col = adjustcolor(1, alpha.f = 0.2),
                           lines.args = list(), xlim = NULL, ylim = NULL,
                           xlab = "Threshold", ylab = NULL,
                           xlab2 = "Excesses", plot = TRUE, ...)
{
    ## Checks
    stopifnot(length(thresholds) >= 2, is.logical(estimate.cov), 0 <= conf.level, conf.level <= 1,
              is.logical(plot))

    ## Compute threshold excesses
    x <- as.numeric(x)
    excesses <- function(u) x[x > u] - u

    ## Peaks-over-threshold for each considered threshold
    exc <- lapply(thresholds, function(u) excesses(u))
    fits <- lapply(exc, fit_GPD_MLE, estimate.cov = estimate.cov) # fit a GPD for each threshold
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
        if(is.null(xlim)) xlim <- range(thresholds)
        if(is.null(ylim))
            ylim <- if(is.na(CI.col) || !estimate.cov) range(xi) else range(xi, xi.CI.low, xi.CI.up)
        if(is.null(ylab)) ylab <- paste0("Estimated GPD shape parameter",
                                         if(estimate.cov) paste0(" with ", 100*conf.level,
                                                                 "% confidence intervals") else "")
        plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...) # plot region
        if(estimate.cov)
            polygon(x = c(thresholds, rev(thresholds)), y = c(xi.CI.low, rev(xi.CI.up)),
                    col = CI.col, border = NA) # CI
        do.call(lines, args = c(list(x = thresholds, y = xi), lines.args)) # actual plot
        pu <- axTicks(1) # used to be pretty(thresholds)
        axis(3, at = pu, labels = sapply(pu, function(u) sum(x > u))) # *corresponding* excesses
        mtext(xlab2, side = 3, line = 3)
    }

    ## Return
    invisible(list(thresholds = thresholds, excesses = exc, GPD.fits = fits))
}


### Plot tail estimator: empirical vs Smith ####################################

##' @title GPD-Based Tail Estimator in the POT Method
##' @param q numeric vector of evaluation points (quantiles)
##' @param threshold threshold u
##' @param p.exceed probability of exceeding the threshold u;
##'        for Smith estimator, this would be mean(x > threshold) for
##'        'x' being the data.
##' @param shape GPD shape parameter xi
##' @param scale GPD scale parameter beta
##' @return Estimator evaluated at q
##' @author Marius Hofert
tail_estimator_GPD <- function(q, threshold, p.exceed, shape, scale)
{
    stopifnot(q >= threshold, 0 <= p.exceed, p.exceed <= 1, scale > 0)
    ## Idea:
    ## 1) bar{F}(x) = bar{F}(u) bar{F}_u(x-u) for x >= u (exceedance).
    ##    Or: bar{F}(u + y) = bar{F}(u) bar{F}_u(y) for y (= x-u) >= 0 (excess)
    ## 2) bar{F}(u) ~= bar{F}_n(u) = Nu/n (number of exceedances / number of all samples)
    ## 3) bar{F}_u(y) ~= bar{F}_GPD(y) for evaluation points y = q - u
    ## => bar{F}(x) = bar{F}(u + y = q) = bar{F}_n(u) bar{F}_{u}(y = q - u) (by 1))
    ##              = Nu/n * pGPD(q - u), shape, scale, lower.tail = FALSE) (by 2), 3))
    p.exceed * pGPD(q - threshold, shape = shape, scale = scale, lower.tail = FALSE)
}

##' @title Plot of the Empirical Tail Estimator (Possibly Overlaid with Smith Tail Estimator)
##' @param x numeric vector of data points
##' @param threshold threshold u (above which tail exceedance df is computed)
##' @param shape NULL or GPD shape parameter xi (typically obtained from fit_GPD_MLE())
##' @param scale NULL or GPD scale parameter beta (typically obtained from fit_GPD_MLE())
##' @param q NULL or a sequence of quantiles >= threshold where to evaluate
##'        the Smith estimator. Can also be of length 1 in which case the plot
##'        is extended to the right of the largest exceedance until q if q is larger
##'        than the largest exceedance.
##' @param length.out length of q
##' @param lines.args list providing additional arguments for the underlying
##'        lines()
##' @param log logical indicating whether logarithm scale is used
##' @param xlim x-axis limits
##' @param ylim y-axis limits
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param ... additional arguments passed to the underlying plot()
##' @return invisible()
##' @author Marius Hofert
##' @note - Using points is good (no lines) as this is more a Q-Q plot or
##'         mean excess plot (not comparable with edf_plot())
##'       - QRM::plotTail() evaluated inverse of Smith estimator at separate
##'         ppoints() (ppoints.gpd()), but that doesn't always match well with
##'         the data: QRM::plotTail(QRM::fit.GPD(fire, threshold = 50), ppoints.gpd = ppoints(4))
tail_plot <- function(x, threshold, shape = NULL, scale = NULL,
                      q = NULL, length.out = 129, lines.args = list(),
                      log = "xy", xlim = NULL, ylim = NULL,
                      xlab = "x", ylab = "Tail probability at x", ...)
{
    ## Compute empirical tail estimator
    ## Note: Unless evaluated at the excesses themselves, there is no point in using
    ##       the EVT-based 'decomposition'.
    ## Idea:
    ## 1) bar{F}(x) = bar{F}(u) bar{F}_u(x-u) for x >= u (exceedance).
    ##    Or: bar{F}(u + y) = bar{F}(u) bar{F}_u(y) for y (= x-u) >= 0 (excess)
    ## 2) bar{F}(u) ~= bar{F}_n(u) = Nu/n (number of exceedances / number of all samples)
    ## 3) bar{F}_u(y) ~= bar{F}_{u,n}(y) => bar{F}_u(y_{(i)}) ~= 1-F_{u,n}(y_{(i)})
    ##    = 1-(1:Nu)/Nu ~= rev(ppoints(Nu)), where y_{(1)} <= ... <= y_{(Nu)} are
    ##    the sorted excesses.
    ## => bar{F}_n(x_{(i)}) = bar{F}_n(u + y_{(i)}) = bar{F}_n(u) bar{F}_{u,n}(y_{(i)}) (by 1))
    ##                      = Nu/n * rev(ppoints(Nu)) (by 2), 3))
    x <- as.numeric(x)
    exceed <- sort(x[x > threshold]) # sorted exceedances; x-values later
    Nu <- length(exceed) # number of exceedances
    Fn.bar.u <- Nu/length(x) # tail estimate at threshold u (= bar{F}_n(u))
    Fn.bar.exceed <- Fn.bar.u * rev(ppoints(Nu)) # bar{F}_n(x_{(i)}); see above; y-values later

    ## Compute GPD tail estimator
    doGPD <- !is.null(shape) && !is.null(scale)
    if(doGPD) {
        stopifnot(scale > 0)
        ## Determine q (of exceedances where to evaluate Smith estimator; range should match 'exceed')
        qlen <- length(q)
        if(is.null(q) || qlen == 1) {
            ## Could extend q[1] 'to the left' if qlen == 2, but that doesn't make
            ## much sense since GPD(x-u) is 0 below the threshold.
            up <- if(qlen == 1) max(exceed[Nu], q) else exceed[Nu]
            q <- seq2(exceed[1], up, length.out = length.out, log = grepl("x", log))
        }
        ## Compute Smith estimator (classical GPD-based tail estimator in the POT method)
        y.Smith <- tail_estimator_GPD(q, threshold = threshold, p.exceed = Fn.bar.u,
                                      shape = shape, scale = scale)
        if(is.null(xlim)) xlim <- range(exceed, q) # only equal to range(exceed) if q is not user-provided
        if(is.null(ylim)) ylim <- range(Fn.bar.exceed, y.Smith)
    } else {
        if(is.null(xlim)) xlim <- range(exceed)
        if(is.null(ylim)) ylim <- range(Fn.bar.exceed)
    }
    ## Plot
    plot(exceed, Fn.bar.exceed, xlim = xlim, ylim = ylim, log = log, xlab = xlab, ylab = ylab, ...)
    if(doGPD) {
        do.call(lines, args = c(list(x = q, y = y.Smith), lines.args))
        invisible(list(np    = cbind(exceed = exceed, Fn.bar.exceed = Fn.bar.exceed),
                       Smith = cbind(q = q, estimator = y.Smith)))
    } else invisible(cbind(exceed = exceed, Fn.bar.exceed = Fn.bar.exceed))
}
