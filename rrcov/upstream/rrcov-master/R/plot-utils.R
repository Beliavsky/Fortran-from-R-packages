######
##  VT::14.01.2020
##
##
##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/R/rrcov", load_code=roxygen2:::load_installed)
##
#' Calculates the points for drawing a confidence ellipsoid
#'
#' @description A simple function to calculate the points of
#'   a confidence ellipsoid, by default \code{dist=qchisq(0.975, 2)}
#' @param loc location vector
#' @param cov a \code{pXp} covariance matrix
#' @param crit the confidence level, default is \code{crit=0.975}
#' @return  A matrix with two columns containing the calculated points.
#'
#' @examples
#'
#' data(hbk)
#' cc <- cov.wt(hbk)
#' e1 <- getEllipse(loc=cc$center[1:2], cov=cc$cov[1:2,1:2])
#' e2 <- getEllipse(loc=cc$center[1:2], cov=cc$cov[1:2,1:2], crit=0.99)
#' plot(X2~X1, data=hbk,
#'     xlim=c(min(X1, e1[,1], e2[,1]), max(X1,e1[,1], e2[,1])),
#'     ylim=c(min(X2, e1[,2], e2[,2]), max(X2,e1[,2], e2[,2])))
#' lines(e1, type="l", lty=1, col="red")
#' lines(e2, type="l", lty=2, col="blue")
#' legend("topleft", legend=c(0.975, 0.99), lty=1:2, col=c("red", "blue"))
#'
#' @export
#' @author Valentin Todorov, \email{valentin.todorov@@chello.at}
#'
getEllipse <- function(loc=c(0, 0), cov=matrix(c(1,0,0,1), ncol=2), crit=0.975) {
    ## A simple function to calculate the points of a confidence ellipsoid,
    ##  by default \code{dist=qchisq(0.975, 2)}
    ## input: data set location and covariance estimate, cutoff

    if (length(d <- dim(cov)) != 2 || (p <- d[1]) != d[2])
        stop("'cov' must be p x p  cov-matrix defining an ellipsoid")

    dist <- sqrt(qchisq(crit, p))
    A <- solve(cov)
    eA <- eigen(A)
    ev <- eA$values

    lambda1 <- max(ev)
    lambda2 <- min(ev)
    eigvect <- eA$vectors[, order(ev)[2]]

    z <- seq(0, 2 * pi, by = 0.01)
    z1 <- dist/sqrt(lambda1) * cos(z)
    z2 <- dist/sqrt(lambda2) * sin(z)
    alfa <- atan(eigvect[2]/eigvect[1])
    r <- matrix(c(cos(alfa),  - sin(alfa), sin(alfa), cos(alfa)), ncol = 2)
    z <- t(t(cbind(z1, z2) %*% r) + loc)    #   xmin <- min(x, z[, 1])

    z
}

## Internal class defining a default legend parameters to be used
##  in plots with both robust and classical lines in the same panel
##  (e.g. scree plot, tolerance ellipses)
setClass(".Legend", representation( leg = "logical",        # whether to draw a legend
                                    pos = "character",      # position of the legend
                                    txt = "vector",         # texts
                                    col = "vector",         # colors
                                    lty = "vector",         # line styles
                                    pch = "vector"),        # characters
                    prototype = list(leg = TRUE,
                                     pos = "topright",
                                     txt = c("robust", "classical"),
                                     col = c("red", "blue"),
                                     lty = c("solid", "solid"),
                                     pch = c(1,24)))


myscreeplot <- function(rcov, ccov) {
    .myscreeplot(rcov, ccov)
}

.myscreeplot <- function(rcov, ccov) {
    er <- NULL
    ec <- NULL
    ee <- NULL
    if(!missing(ccov))
        ee <- ec <- getEvals(ccov)
    if(!missing(rcov))
        ee <- er <- getEvals(rcov)
    if(is.null(ee))
        stop("Both parameters are NULL")

    leg <- new(".Legend")
    eall <- if(!is.null(ec) && !is.null(er))
                c(er, ec)
            else if(!is.null(ec))
                ec
            else er

    ylim <- c(min(eall), max(eall))

    plot(ee, ylim=ylim, ylab="Eigenvalues", xlab="Index", type="n")
    if(!is.null(ec) && !is.null(er))
        legend(leg@pos, leg@txt, pch = leg@pch, lty = leg@lty, col = leg@col)

    if(!is.null(er))
        lines(er, type="o", pch = leg@pch[1], lty = leg@lty[1], col=leg@col[1])
    if(!is.null(ec))
        lines(ec, type="o", pch = leg@pch[2], lty = leg@lty[2], col=leg@col[2])

    title(main="Scree plot")
}

.myddplot <- function(md, rd, cutoff, id.n,
    main="Distance-Distance Plot",
    xlab="Mahalanobis distance",
    ylab="Robust distance",
    labs=1:length(md),
    ...)
{
    ##  Distance-Distance Plot:
    ##  Plot the vector y=rd (robust distances) against
    ##  x=md (mahalanobis distances). Identify by a label the id.n
    ##  observations with largest rd. If id.n is not supplied, calculate
    ##  it as the number of observations larger than cutoff. Use cutoff
    ##  to draw a horisontal and a vertical line. Draw also a dotted line
    ##  with a slope 1.
    n <- length(md)
    if(missing(id.n))
        id.n <- length(which(rd>cutoff))

    plot(md, rd, xlab=xlab, ylab=ylab, type="p", ...)
    .label(md, rd, id.n, labs=labs)
    abline(0, 1, lty=2)
    abline(v=cutoff)
    abline(h=cutoff)

    title(main=main)
}

.mydistplot <- function(x, cutoff, classic = FALSE, id.n,
    ylim=NULL,
    main="Distance Plot",
    xlab="Index",
    ylab=paste(if(classic) "Mahalanobis" else "Robust", "distance"),
    labs=1:length(x),
    ...)
{
    ## VT::10.11.2007 - change "Squared Robust distance" to "Robust distance"
    ## VT::10.11.2007 - Add parameter ylim to make possible robust and
    ##  classical plot to use the same y-scale

    ##  Index Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the observation indexes. Identify by a label the id.n
    ##  observations with largest value of x. If id.n is not supplied,
    ##  calculate it as the number of observations larger than cutoff.
    ##  Use cutoff to draw a horisontal line.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes

    n <- length(x)
    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    plot(x, ylab=ylab, xlab=xlab, type="p", ylim=ylim, ...)
    .label(1:n, x, id.n, labs=labs)
    abline(h=cutoff)

    title(main=main)
}

.qqplot <- function(x, p, cutoff, classic=FALSE, id.n,
    main=eval(substitute(expression(paste(chi^2, " QQ-Plot")))),
    xlab=eval(substitute(expression(paste("Sqrt of the quantiles of the ", chi[p]^2, " distribution")), list(p=p))),
    ylab=paste(if(classic) "Mahalanobis" else "Robust", "distance"),
    labs=1:length(x),
    ...)
{
    ##  Chisquare QQ-Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the square root of the quantiles of the chi-squared distribution
    ##  with p degrees of freedom.
    ##  Identify by a label the id.n observations with largest value of x.
    ##  If id.n is not supplied, calculate it as the number of observations
    ##  larger than cutoff.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes


    ##  parameters and preconditions

    n <- length(x)

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

    x <- sort(x, index.return=TRUE)
    ix <- x$ix
    x <- x$x

    ## xlab <- "Square root of the quantiles of the chi-squared distribution"
    plot(qq, x, xlab=xlab, ylab=ylab, type="p", ...)
    if(id.n > 0) {
        ind <- (n-id.n+1):n
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(qq[ind] + xrange/50, x[ind], labs[ix[ind]])
    }
    abline(0, 1, lty=2)
    title(main=main)
}

.label <- function(x, y, id.n=3, labs=1:length(x)) {
    if(id.n > 0) {
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        n <- length(y)
        ind <- sort(y, index.return=TRUE)$ix
        ind <- ind[(n-id.n+1):n]
        text(x[ind] + xrange/40, y[ind], labs[ind], cex=0.9)
    }
}

.tolellipse <- function(rcov, ccov, cutoff = NULL, id.n = NULL, tol = 1e-07,
    main = "Tolerance ellipse (97.5%)",
    xlab = "",
    ylab = "",
    labs, legend.position="bottomright",
    ...)
{

    leg <- new(".Legend")

    if(missing(rcov) || is.null(rcov)){
        if(missing(ccov) || is.null(ccov))
            stop("Location and scatter matrix must be provided!")

        ## there is no robust location/scatter object
        rcov <- ccov
        leg@leg <- FALSE
    }
    if(is.null(data <- getData(rcov)))
        stop("No data provided!")
    n <- dim(data)[1]
    p <- dim(data)[2]
    if(p != 2)
        stop("Dimension must be 2!")

    r.cov <- getCov(rcov)
    r.loc <- getCenter(rcov)
    if(length(r.loc) == 0 ||  length(r.cov) == 0)
        stop("Invalid 'rcov' object: attributes center and cov missing!")
    z1 <- getEllipse(loc = r.loc, cov = r.cov)
    rd <- sqrt(getDistance(rcov))
    x1 <- c(min(data[, 1], z1[, 1]), max(data[,1],z1[,1]))
    y1 <- c(min(data[, 2], z1[, 2]), max(data[,2],z1[,2]))
    classic <- FALSE

    if(!missing(ccov) && !is.null(ccov)){
        c.cov <- getCov(ccov)
        c.loc <- getCenter(ccov)
        if(length(c.loc) == 0 ||  length(c.cov) == 0)
            stop("Invalid 'ccov' object: attributes center and cov missing!")
        classic <- TRUE
        z2 <- getEllipse(loc=c.loc, cov=c.cov)
        md <- sqrt(getDistance(ccov))
        x1 <- c(min(data[, 1], z1[, 1], z2[, 1]), max(data[,1],z1[,1], z2[,1]))
        y1 <- c(min(data[, 2], z1[, 2], z2[, 2]), max(data[,2],z1[,2], z2[,2]))
    }

    ## Note: the *calling* function may pass a 'missing' value
    if(missing(cutoff) || is.null(cutoff))
        cutoff <- sqrt(qchisq(0.975, df = 2))
    if(missing(id.n) || is.null(id.n))
        id.n <- sum(rd > cutoff)
    if(missing(labs) || is.null(labs))
        labs <- 1:length(rd)

    ind <- sort(rd, index.return=TRUE)$ix
    ind <- ind[(n-id.n+1):n]

##  1. Robust tolerance
##  define the plot, plot a box, plot the "good" points,
##  plot the outliers either as points or as numbers depending on outflag,
##  plot the ellipse, write a title of the plot
    plot(data, xlim = x1, ylim = y1, xlab = xlab, ylab = ylab, ...)
    box()

    ## VT:::03.08.2008
    if(id.n > 0){
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(data[ind, 1] + xrange/50, data[ind, 2], labs[ind])
    }

    if(leg@leg)
        points(z1, type = "l", lty=leg@lty[1], col=leg@col[1])
    title(main)

##  2. Classical tolerance ellipse and legend
    if(classic){
        points(z2, type = "l", lty=leg@lty[2], col=leg@col[2])
        if(leg@leg)
            legend(legend.position, leg@txt, pch=leg@pch, lty=leg@lty, col=leg@col)
    }

    invisible()
}

## Plot a scatter plot of the data 'x' and superimpose a
##  (cutoff-)tollernce ellipse.
##  The differences to tolEllipsePlot() in robustbase are:
##  - customizable (titles, limits, labels, etc)
##  - can take either a Cov object or a list (aka cov.wt or covMcd)
##  -
.myellipse <- function (x, xcov,
    cutoff = NULL,
    id.n = NULL,
    classic = FALSE,
    tol = 1e-07,
    xlab = "",
    ylab = "",
    main="Tolerance ellipse (97.5%)",
    sub="",
    leg.txt=c("robust", "classical"), leg.col=c("red", "blue"),
    leg.lty=c("solid", "dashed"),
    xlim, ylim,
    off.x, off.y,
    ...)
{
    if (is.data.frame(x))
        x <- data.matrix(x)
    if (!is.matrix(x) || !is.numeric(x))
        stop("x is not a numeric dataframe or matrix.")
    n <- dim(x)[1]
    p <- dim(x)[2]
    if (p != 2)
        stop("Dimension {= ncol(x)} must be 2!")

    if(missing(xcov))
        xcov <- CovMcd(x)
    if(is(xcov, "Cov"))
    {
        x.loc <- getCenter(xcov)
        x.cov <- getCov(xcov)
    } else {
        if (!is.numeric(xcov$center) || !is.numeric(xcov$cov))
            stop("argument 'xcov' must have numeric components 'center' and 'cov'")

        x.loc <- xcov$center
##        x.cov <- n/(n - 1) * xcov$cov
        x.cov <- xcov$cov
    }

    xM <- colMeans(x)
    z1 <- getEllipse(loc=xM, cov=n/(n - 1) * cov.wt(x)$cov)
    z2 <- getEllipse(loc=x.loc, cov=x.cov)
    ## VT::09.06.200
    ## if no need of the classic ellipse - set it to the robust one
    ## otherwise the xlim and ylim will be set to fit both ellipses

##    if(classic == FALSE)
##        z1 <- z2

    x1 <- c(min(x[, 1], z1[, 1], z2[, 1]), max(x[, 1], z1[, 1],
        z2[, 1]))
    y1 <- c(min(x[, 2], z1[, 2], z2[, 2]), max(x[, 2], z1[, 2],
        z2[, 2]))
    md <- sqrt(mahalanobis(x, xM, cov(x), tol = tol))
    rd <- sqrt(mahalanobis(x, x.loc, x.cov, tol = tol))
    if (missing(cutoff) || is.null(cutoff))
        cutoff <- sqrt(qchisq(0.975, df = 2))
    if (missing(id.n) || is.null(id.n)) {
        id.n <- sum(rd > cutoff)
        if(n <= 10)             # if less than 10 observations, show all
            id.n <- n
    }

    ind <- sort(rd, index.return = TRUE)$ix
    ind <- ind[(n - id.n + 1):n]

    if(missing(xlim))
        xlim <- x1
    if(missing(ylim))
        ylim <- y1

    plot(x, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main,
        sub=sub, ...)

    box()
    
    ## Add the labels
    if(id.n > 0) {
        if(missing(off.x))
        {
            xrange <- par("usr")
            xrange <- xrange[2] - xrange[1]
            off.x <- xrange/50
        }
        if(missing(off.y))
        {
            xrange <- par("usr")
            xrange <- xrange[4] - xrange[3]
            off.y <- xrange/50
        }
        labels <- if(!is.null(rownames(x))) rownames(x)[ind] else ind
        text(x[ind, 1] + off.x, x[ind, 2] + off.y, labels)
    }

    ## Draw the ellipse(s)
    points(z2, type="l", lty=leg.lty[1], col=leg.col[1])
    if(classic) {
        points(z1, type="l", lty=leg.lty[2], col=leg.col[2])
        legend("bottomright", leg.txt, lty=leg.lty, col=leg.col)
    }
    
    invisible()
}

##  Draw pairwise scatter plots for the data set 'x'
##  - upper triangle - scatter plot with classical and robust 0.975-ellipses
##  - histograms on the diagonal
##  - lower triangle - robust (MCD) and classical correlations
##
##  - x     - data
##  - main  - caption of the plot
##
.rrpairs <- function(obj, main="", sub="", xlab="", ylab="", ...){
    hcol       <- "#377eb8"   # colour for histogram - dark blue
    dcol       <- "red"       # color of the density line
    ecol.class <- "blue"      # colour for classical ellipse
    ecol.rob   <- "red"       # colour for robust ellipse

    ## quick and dirty simulation of panel.number() which seems not to work
    ##  in the panel functions of pairs()
    ##  Retursns list of the corresponding i and j
    ##
    ##  Example call: which.ij(hbk[,1], hbk[,3], getData(CovMcd(hbk)))
    ##
    which.ij <-function(x, y, data)
    {
        i <- j <- 0
        for(k in 1:ncol(data))
        {
            ifi <- all.equal(x, data[,k], check.attributes=FALSE)
            ifj <- all.equal(y, data[,k], check.attributes=FALSE)
            if(i == 0 && !is.character(ifi) && ifi)
                i <- k
            if(j == 0 && !is.character(ifj) && ifj)
                j <- k
            if(i != 0 & j != 0)
                break
        }
        list(i=i, j=j)
    }

    panel.hist <- function(x, ...)
    {
        usr <- par("usr"); on.exit(par(usr=usr))
        par(usr = c(usr[1:2], 0, 1.5) )

        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col=hcol, ...)
    }

    panel.hist.density <- function(x,...)
    {
        usr <- par("usr"); on.exit(par(usr=usr))
        par(usr = c(usr[1:2], 0, 1.5) )

        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col=hcol)

        tryd <- try( d <- density(x, na.rm=TRUE, bw="nrd", adjust=1.2), silent=TRUE)

        ## VT::11.08.2022: fix error "Found if() conditions comparing class() to string"
        ##  if(class(tryd) != "try-error")
        if(!is(tryd, "try-error"))
        {
            d$y <- d$y/max(d$y)
            lines(d, col=dcol)
        }
    }

    panel.cor <- function(x, y, digits=2, ...)
    {
        ix <- which.ij(x, y, getData(obj))

        usr <- par("usr"); on.exit(par(usr=usr))
        par(usr = c(0, 1, 0, 1))

        r <- cor(x, y)
        rr <- getCorr(obj)[ix$i,ix$j]

        prefix <- ""
        rprefix <- ""
        ff <- 0.35

        txt  <- format(c(r, 0.123456789), digits=digits)[1]
        txt  <- paste(prefix, txt, sep="")
        rtxt <- format(c(rr, 0.123456789), digits=digits)[1]
        rtxt <- paste(rprefix, rtxt, sep="")
        txt  <- if(isClassic(obj)) txt else paste0("(", txt, ")")   # parentheses around the classical correlation
                                                                    # only if the object is a robust one
        cex  <- ff/strwidth(txt)
        text(0.5, 0.3, txt, cex = cex, col=ecol.class)
        if(!isClassic(obj))
            text(0.5, 0.5, rtxt, cex = cex, col=ecol.rob)
    }

    panel.ellipse <- function(x, y, ...)
    {
        usr <- par("usr"); on.exit(par(usr=usr))

        X <- cbind(x,y)
        C.ls <- cov(X)
        m.ls <- colMeans(X)

        ix <- which.ij(x, y, getData(obj))
        ## cat("\npanel.ellipse: ", ix$i, ix$j,"\n")

        C.rr <- getCov(obj)[c(ix$i,ix$j), c(ix$i,ix$j)]
        m.rr <- getCenter(obj)[c(ix$i,ix$j)]

        e.class <- getEllipse(loc=m.ls, cov=C.ls, crit=0.975)
        e.rob <- getEllipse(loc=m.rr, cov=C.rr, crit=0.975)

        xmin <- min(c(min(x), min(e.class[,1]), min(e.rob[,1])))
        xmax <- max(c(max(x), max(e.class[,1]), max(e.rob[,1])))
        ymin <- min(c(min(y), min(e.class[,2]), min(e.rob[,2])))
        ymax <- max(c(max(y), max(e.class[,2]), max(e.rob[,2])))

        ff <- 0.1
        xmin <- xmin - ff*(xmax-xmin); #print(ff*(xmax-xmin))
        xmax <- xmax + ff*(xmax-xmin); #print(ff*(xmax-xmin))
        ymin <- ymin - ff*(ymax-ymin); #print(ff*(ymax-ymin))
        ymax <- ymax + ff*(ymax-ymin); #print(ff*(ymax-ymin))
        par(usr = c(xmin, xmax, ymin, ymax))

        points(x,y, ...)
        lines(e.class, col=ecol.class, lty="dashed")
        if(!isClassic(obj))
            lines(e.rob, col=ecol.rob)
    }

    ## get the data
    x <- getData(obj)

    ## VT::27.04.2011 - fix the names of the variables on the
    ##  diagonal - labels= in the call to pairs().
    ##  Plot also density line
    ##
    pairs(x, main = main, sub=sub,
        lower.panel=panel.cor,
        diag.panel=panel.hist.density,
        upper.panel=panel.ellipse,
        labels=names(getCenter(obj)),
        ...)
}

## Distance plot:
##  Plot the robust and classical distances against the the index
##  obj - A CovRobust object,
##  getData(obj)    - data frame or matrix
##
.xydistplot <- function(obj, cutoff,
    main="Distance Plot",
    xlab="Index",
    ylab="Mahalanobis distance",
    col="darkred",
    labs,
    ...)
{
    myPanel <- function(x, y, subscripts, cutoff, id.n, ...) {
        panel.xyplot(x, y, ...)
        panel.abline(h=cutoff,lty="dashed")

        n <- length(y)
        if(missing(id.n))
            id.n <- length(which(y > cutoff))

        if(id.n > 0){
            ind <- sort(y, index.return=TRUE)$ix
            ind <- ind[(n-id.n+1):n]

            xrange<-range(y)
            adj <- (xrange[2]-xrange[1])/20
            ltext(x[ind] + adj, y[ind] + adj, ind, cex=0.85)
        }
    }

    X <- getData(obj)
    n <- nrow(X)
    p <- ncol(X)
    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))
    if(missing(labs) || is.null(labs))
        labs <- 1:length(getDistance(obj))

    ## VT::12.11.2018 - In case of MRCD we can use also classic regularized estimate
    dd1 <- sqrt(getDistance(obj))                   # robust distances
    if(inherits(obj, "CovMrcd"))
    {
        vv <- CovMrcd(X, alpha=1)                   # classical REGULARIZED center and covariance
        dd2 <- sqrt(getDistance(vv))                # classical distances
    } else
    {
        vv  <- cov.wt(X)
        dd2 <- sqrt(mahalanobis(X,vv$center,vv$cov))    # classical distances
    }
    dd  <- c(dd1, dd2)                               # a vector with both robust and classical distances

    ind <- c(1:n, 1:n)                              # 1, 2, ..., n, 1, 2, ...n      -
    gr  <- as.factor(c(rep(1,n), rep(2,n)))         # 1, 1, ..., 1, 2, 2, ...2      -   n x 1, n x 2
    levels(gr)[1] <- "Robust"
    levels(gr)[2] <- "Classical"

    xyplot(dd~ind|gr,
                cutoff=cutoff,
                panel = myPanel,
                xlab=xlab,
                ylab=ylab,
                main=main,
                col=col,
                ...)
}

## QQ-Chi plot:
##  Plot QQ plot of the robust and classical distances against the
##  quantiles of the Chi2 distr
##  X - data frame or matrix
##
.xyqqchi2 <- function(obj, cutoff,
    main=eval(substitute(expression(paste(chi^2, " QQ-Plot")))),
    xlab=eval(substitute(expression(paste("Sqrt of the quantiles of the ", chi[p]^2, " distribution")), list(p=p))),
    ylab="Mahalanobis distance",
    col="darkred",
    labs,
    ...)
{
    myPanel <- function(x, y, subscripts, cutoff, id.n, ...)
    {
        y <- sort(y, index.return=TRUE)
        iy <- y$ix
        y <- y$x
        panel.xyplot(x, y, ...)
        panel.abline(0,1,lty="dashed")

        n <- length(y)
        if(missing(id.n))
            id.n <- length(which(y > cutoff))
        if(id.n > 0){
            ind <- (n-id.n+1):n

            xrange<-range(y)
            adj <- (xrange[2]-xrange[1])/50
            ltext(x[ind] + adj, y[ind] + adj, iy[ind], cex=0.85)
        }
    }

    X <- getData(obj)
    n <- nrow(X)
    p <- ncol(X)
    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))
    if(missing(labs) || is.null(labs))
        labs <- 1:length(getDistance(obj))

    ## VT::12.11.2018 - In case of MRCD we can use also classic regularized estimate
    dd1 <- sqrt(getDistance(obj))                   # robust distances
    if(inherits(obj, "CovMrcd"))
    {
        vv  <- CovMrcd(X, alpha=1)                   # classical REGULARIZED center and covariance
        dd2 <- sqrt(getDistance(vv))                 # classical Mahalanobis distances
    } else
    {
        vv  <- cov.wt(X)                             # classical center and covariance
        dd2 <- sqrt(mahalanobis(X,vv$center,vv$cov)) # classical Mahalanobis distances
    }
    dd  <- c(dd1, dd2)                               # a vector with both robust and classical distances

    qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))
    ind<-c(qq, qq)
    gr<-as.factor(c(rep(1,n), rep(2,n)))             # 1, 1, ...., 1, 2, 2, ..., 2   - n x 1, n x 2
    levels(gr)[1]<-"Robust"
    levels(gr)[2]<-"Classical"

    xyplot(dd~ind|gr,
        cutoff=cutoff,
        panel = myPanel,
        xlab=xlab,
        ylab=ylab,
        main=main,
        col=col,
        ...)
}
