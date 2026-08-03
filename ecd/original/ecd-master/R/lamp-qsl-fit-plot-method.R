#' Plot the fit to asset returns using quartic stable lambda distribution
#'
#' Plot the fit to asset returns using quartic stable lambda distribution
#'
#' @param key character, the key(s) to retrieve configuration for plot.
#'            If key is provided as single-row data frame, then it will be used directly as config.
#' @param debug logical, if true, print debug information. Default is FALSE.
#' @param plot.type character, type of plot: pdf, log-pdf, qqplot. Default is \code{c("pdf","log-pdf")}.
#' @param qqplot.n numeric, specify number of QSLD simulations for qqplot utility, default is 1000000.
#' @param extdata_dir optionally specify user's own extdata folder, default is NULL.
#' @param filename character, optionally specify user's own config file name, default is NULL.
#'
#' @return returns a list of each key and its data and config blocks as nested list
#'
#' @keywords QSLD
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom moments kurtosis
#' @importFrom moments skewness
#' @importFrom stats qqplot
#'
#' @export lamp.qsl_fit_plot
#'
### <======================================================================>
lamp.qsl_fit_plot <- function(key, debug=FALSE,
                                   plot.type=c("pdf", "log-pdf"), qqplot.n=1000000,
                                   extdata_dir=NULL, filename=NULL) {
    c <- NULL
    if (is(key, "character")) {
        if(length(key) > 1) {
            dt <- list()
            for (k in key) {
                dt1 <- lamp.qsl_fit_plot(k, debug=debug, plot.type=plot.type,
                                              extdata_dir=extdata_dir, filename=filename)
                dt[k] <- list(dt1[[k]])
            }
            return(invisible(dt))
        }
        c <- lamp.qsl_fit_config(key, extdata_dir=extdata_dir, filename=filename)
    }
    else c <- key # else we assume period contains config

    if (! is(c, "data.frame")) stop("config is not a data frame")
    if (debug) print(c)
    if (NROW(c) != 1) stop("config is not a single-row data frame")

    symbol <- c$symbol
    t <- c$t
    nu0 <- c$nu0
    theta <- c$theta
    convo <- c$convo
    beta.a <- c$beta.a
    start.date <- if (!is.na(c$start.date)) c$start.date else "1950-01-01"
    end.date <- if (!is.na(c$end.date)) c$end.date else "2015-12-31"

    idx <- ecd.data.arr(symbol, on=c$ts.on, lag=c$lag, drop=c$drop,
                        start.date=start.date, end.date=end.date)
    kgq <- kqsl(t=t, nu0=nu0/100, theta=theta/100, convo=convo, beta.a=beta.a)
    kll <- kstdlap(t=t, convo=convo, beta=beta.a*sqrt(t))

    x <- idx$x
    mean_x <- mean(x)
    sd_x <- sd(x)
    #
    mu <- mean_x - kgq[1]

    # plot
    h <- hist(x, breaks=c$breaks, plot=FALSE)

    px0 <- seq(min(h$mids), max(h$mids), length.out=1000)
    px <- px0 - mu
    p <- dqsl(px, t=t, nu0=nu0/100, theta=theta/100, convo=convo, beta.a=beta.a)
    pdf0 <- max(p)

    h$density[is.infinite(log(h$density))] <- NaN
    yh <- log(h$density)
    ymax <- max(yh, na.rm=TRUE)
    ymin <- min(log(h$density), na.rm=TRUE)


    .bounded_k34 <- function(x, pdf) {
        dx <- diff(x)[1]
        m1 <- sum(x*pdf*dx)
        m2 <- sum((x-m1)^2*pdf*dx)
        m3 <- sum((x-m1)^3*pdf*dx)
        m4 <- sum((x-m1)^4*pdf*dx)
        c(m3/m2^(3/2), m4/m2^2) # skewness, kurtosis
    }
    k34b <- .bounded_k34(px, p)

    .plot_sd_lines <- function(n=c(1,2), col="green") {
        abline(v=mean_x, lty=2, col=col)
        for (i in 1:6) {
            if (i %in% n) {
                abline(v=mean_x+i*sd_x,   lty=2, col=col)
                abline(v=mean_x-i*sd_x,   lty=2, col=col)
            }
        }
    }
    # PDF
    if ("pdf" %in% plot.type) {
        plot(h$mids, h$density, cex=0.4, col="blue",
             xlim=c(-4,4)*sd_x+mean_x,
             ylim=c(0, exp(ymax)),
             xlab="log-returns", ylab="PDF",
             main=sprintf("Index: %s Period: %s", symbol, c$period))

        lines(px0, p, col="red")
        .plot_sd_lines()

        x0 <- 1.4*sd_x
        y0 <- exp(ymax)
        sz <- 0.8
        ss <- "$\\widehat{P}_{\\varLambda0}$"
        if (debug) ss <- "PDF(0)"
        text(x0, y0*0.90, pos=4, cex=sz, label=sprintf("data stats"))
        text(x0, y0*0.80, pos=4, cex=sz, label=sprintf("mean  %.4f", mean_x))
        text(x0, y0*0.70, pos=4, cex=sz, label=sprintf("stdev %.4f", sd_x))
        text(x0, y0*0.60, pos=4, cex=sz, label=sprintf("skewness %.2f", moments::skewness(x)))
        text(x0, y0*0.50, pos=4, cex=sz, label=sprintf("kurtosis %.2f", moments::kurtosis(x)))
        text(x0, y0*0.40, pos=4, cex=sz, label=sprintf("%s %.3f", ss, exp(ymax)*sd_x))
        text(x0, y0*0.30, pos=4, cex=sz, label=sprintf("years %d", floor(length(x)/250)))

        x0 <- -3.8*sd_x
        y0 <- exp(ymax)
        text(x0, y0*0.90, pos=4, cex=sz, label=sprintf("fit stats"))
        text(x0, y0*0.80, pos=4, cex=sz, label=sprintf("$\\mu=$  %.4f", mu))
        text(x0, y0*0.70, pos=4, cex=sz, label=sprintf("stdev %.4f", sqrt(kgq[2])))
        if (kgq[6] > moments::kurtosis(x)) {
            text(x0, y0*0.60, pos=4, cex=sz, label=sprintf("skewness %.2f / %.2f", kgq[5], k34b[1])) # kll[5]
            text(x0, y0*0.50, pos=4, cex=sz, label=sprintf("kurtosis %.2f / %.2f", kgq[6], k34b[2]))
        } else {
            text(x0, y0*0.60, pos=4, cex=sz, label=sprintf("skewness %.2f", kgq[5])) # kll[5]
            text(x0, y0*0.50, pos=4, cex=sz, label=sprintf("kurtosis %.2f", kgq[6]))
        }
        text(x0, y0*0.40, pos=4, cex=sz, label=sprintf("%s %.3f", ss, pdf0*kgq[2]^0.5))
        text(x0, y0*0.30, pos=4, cex=sz, label=sprintf("$m=$ %.1f", convo))
    }

    # log(PDF)
    if ("log-pdf" %in% plot.type) {
        plot(h$mids, log(h$density), cex=0.4, col="blue",
             ylim=c(ymin-1, ymax),
             xlab="log-returns", ylab="log(PDF)",
             main=sprintf("Index: %s Period: %s", symbol, c$period))
        lines(px0, log(p), col="red")
        .plot_sd_lines()

        y0 <- ymax-0.5
        text(mean_x+2*sd_x, y0, pos=4, cex=sz, label=sprintf("+2sd"))
        text(mean_x-2*sd_x, y0, pos=2, cex=sz, label=sprintf("-2sd"))
    }
    # qqplot
    if ("qqplot" %in% plot.type) {
        # simulate
        Q <- rqsl(qqplot.n, t=t, nu0=nu0/100, theta=theta/100, convo=convo, beta.a=beta.a, mu=mu)
        qset <- h$mids
        qset2 <- c(h$mids,Q)
        stats::qqplot(Q, idx$x, cex=0.4, col="blue",
               xlim=c(min(qset), max(qset)),
               ylim=c(min(qset), max(qset)),
               xlab="Theoretical Quantiles from QSLD (log-rtns)",
               ylab="Empirical Quantiles (log-rtns)",
               main=sprintf("Index: %s Period: %s", symbol, c$period))
        xx <- seq(min(qset2), max(qset2), length.out=200)
        lines(xx, xx, col="red")

        .plot_sd_lines(c(2,4), col="orange")

        sz <- 0.8
        y0 <- min(qset)*0.8
        text(mean_x+4*sd_x, y0, pos=4, cex=sz, label=sprintf("+4sd"))
        text(mean_x-4*sd_x, y0, pos=2, cex=sz, label=sprintf("-4sd"))
    }

    dt <- list()
    idx$conf <- c
    dt[key] <- list(idx)
    return(invisible(dt))
}
### <---------------------------------------------------------------------->
