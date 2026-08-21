#' Compute crude probability of death
#'
#' Estimates the crude probability of death due to disease and due to
#' population reasons
#'
#' NOTE: The follow-up time must be specified in days. The \code{ratetable}
#' being used may have different variable names and formats than the user's
#' data set, this is dealt with by the \code{rmap} argument. For example, if
#' age is in years in the data set but in days in the \code{ratetable} object,
#' age=age*365.241 should be used. The calendar year can be in any date format
#' (Date and POSIXt are allowed), the date formats in the
#' \code{ratetable} and in the data may differ.
#'
#' Note that numerical integration is required to calculate the variance
#' estimator. The integration precision is set with argument \code{precision},
#' which defaults to daily intervals, a default that should give enough
#' precision for any practical purpose.
#'
#' The area under the curve is calculated on the interval [0,\code{tau}].
#'
#' Function \code{summary} may be used to get the output at specific points in
#' time.
#'
#' @aliases cmp.rel print.cmp.rel
#' @param formula a formula object, with the response as a \code{Surv} object
#' on the left of a \code{~} operator, and, if desired, terms separated by the
#' \code{+} operator on the right. If no strata are used, \code{~1} should be
#' specified.
#'
#' NOTE: The follow-up time must be in days.
#' @param data a data.frame in which to interpret the variables named in the
#' \code{formula}.
#' @param ratetable a table of event rates, organized as a \code{ratetable}
#' object, such as \code{slopop}.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param tau the maximum follow-up time of interest, all times larger than
#' \code{tau} shall be censored. Equals maximum observed time by default
#' @param conf.int the level for a two-sided confidence interval on the
#' survival curve(s). Default is 0.95.
#' @param precision the level of precision used in the numerical integration of
#' variance. Default is 1, which means that daily intervals are taken, the
#' value may be decreased to get a higher precision or increased to achieve a
#' faster calculation. The calculation intervals always include at least all
#' times of event and censoring as border points.
#' @param add.times specific times at which the value of estimator and its
#' variance should be evaluated. Default is all the event and censoring times.
#' @param rmap an optional list to be used if the variables are not organized
#' and named in the same way as in the \code{ratetable} object. See details
#' below.
#' @return An object of class \code{cmp.rel}. Objects of this class have
#' methods for the functions \code{print} and \code{plot}. The \code{summary}
#' function can be used for printing output at required time points. An object
#' of class \code{cmp.rel} is composed of several lists, each pertaining the
#' cumulative hazard function for one risk and one strata. Each of the lists
#' contains the following objects: \item{time}{the time-points at which the
#' curves are estimated} \item{est}{the estimate} \item{var}{the variance of
#' the estimate} \item{lower}{the lower limit of the confidence interval}
#' \item{upper}{the upper limit of the confidence interval} \item{area}{the
#' area under the curve calculated on the interval [0,\code{tau}]}
#' \item{index}{indicator of event and censoring times among all the times in
#' the output. The times added via paramater \code{add.times} are also
#' included} \item{add.times}{the times added via parameter \code{add.times}}
#' @seealso \code{rs.surv}, \code{summary.cmp.rel}
#' @references Package: Pohar Perme, M., Pavlic, K. (2018) "Nonparametric
#' Relative Survival Analysis with the R Package relsurv". Journal of
#' Statistical Software. 87(8), 1-27, doi: "10.18637/jss.v087.i08"
#' @keywords survival
#' @examples
#'
#'
#' data(slopop)
#' data(rdata)
#' #calculate the crude probability of death
#' #note that the variable year must be given in a date format and that
#' #age must be multiplied by 365.241 in order to be expressed in days.
#' fit <- cmp.rel(Surv(time,cens)~sex,rmap=list(age=age*365.241),
#' 		ratetable=slopop,data=rdata,tau=3652.41)
#' fit
#' plot(fit,col=c(1,1,2,2),xscale=365.241,xlab="Time (years)")
#' #if no strata are desired:
#' fit <- cmp.rel(Surv(time,cens)~1,rmap=list(age=age*365.241),
#' 		ratetable=slopop,data=rdata,tau=3652.41)
#'
#'
#'
cmp.rel <-  function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop,
     na.action,tau,conf.int=0.95,precision=1,add.times,rmap)

    #formula: for example Surv(time,cens)~1 #not implemented for subgroups - DO IT!
    #data: the observed data set
    #ratetable: the population mortality tables
    #conf.type: confidence interval calculation (plain, log or log-log)
    #conf.int: confidence interval
    #tau: max. cas do katerega racuna
{

    call <- match.call()
    if (!missing(rmap)) {
          rmap <- substitute(rmap)
    }
    rform <- rformulate(formula, data, ratetable, na.action,rmap)			#get the data ready
    data <- rform$data								#the data set
    se.fac <- sqrt(qchisq(conf.int, 1))						#factor needed for confidence interval

    if(missing(tau)) tau<-max(rform$Y)

    p <- rform$m								#number of covariates
    if (p > 0) #if covariates
    data$Xs <- strata(rform$X[, ,drop=FALSE ])					#make strata according to covariates
    else data$Xs <- rep(1, nrow(data))						#if no covariates, just put 1

    tab.strata <- table(data$Xs)						#unique strata values
    ntab.strata <- length(tab.strata)						#number of strata

    dtemp <- list(NULL)
    out <- as.list(rep(dtemp,ntab.strata*2))

    for (kt in 1:ntab.strata) {							#for each stratum
        inx <- which(data$Xs == names(tab.strata)[kt])				#individuals within this stratum

	extra <- as.numeric(seq(1,max(rform$Y[inx]),by=precision))
	if(!missing(add.times)) extra <- c(extra,as.numeric(add.times))
	tis <- sort(unique(pmin(tau,union(rform$Y[inx],extra))) )	#1-day long intervals used - to take into the account the continuity of the pop. part

	#if(!all.times)tis <- sort(unique(pmin(rform$Y[inx],tau)))					#unique times
	#else{
	#	tis <- sort(union(rform$Y[inx], as.numeric(1:floor(max(rform$Y[inx]))))) 	#1-day long intervals used - to take into the account the continuity of the pop. part
	#	tis <- unique(pmin(tis,tau))
   	#}
	k <- length(tis)

	out[[2*kt-1]]$time <- out[[2*kt]]$time <- c(0,tis)

	temp <- exp_prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=TRUE,cmp=T)		#calculate the values for each interval of time

	areae <- sum(temp$areae)/365.241		# sum(diff(c(0,tis))*temp$cumince)/365.241
	areap <- sum(temp$areap)/365.241		#sum(diff(c(0,tis))*temp$cumincp)/365.241

  	options(warn=-1)
  	out[[2*kt-1]]$est <- c(0,temp$cumince)
  	out[[2*kt-1]]$var <- c(0,temp$ve)
  	out[[2*kt-1]]$lower <- temp$cumince-se.fac*sqrt(temp$ve)
  	out[[2*kt-1]]$upper <- temp$cumince+se.fac*sqrt(temp$ve)
  	out[[2*kt-1]]$area <- areae

  	out[[2*kt]]$est <- c(0,temp$cumincp)
	out[[2*kt]]$var <- c(0,temp$vp)
	out[[2*kt]]$lower <- temp$cumincp-se.fac*sqrt(temp$vp)
	out[[2*kt]]$upper <- temp$cumincp+se.fac*sqrt(temp$vp)
	out[[2*kt]]$area <- areap
	options(warn=0)

	ne <- sum(temp$ve<0)
	if(ne>0) warning(paste(names(tab.strata)[kt],": The estimated variance of crude mortality is negative in ", ne, " out of ", length(temp$ve)," intervals"), call. = FALSE)

	if(!missing(add.times)){
	out[[2*kt-1]]$index <- out[[2*kt]]$index <- unique(c(1,which(tis %in% c(rform$Y[inx],add.times,tau))))
	out[[2*kt-1]]$add.times <- out[[2*kt]]$add.times <- add.times
	}
	else out[[2*kt-1]]$index <- out[[2*kt]]$index <- unique(c(1,which(tis %in% c(rform$Y[inx],tau))))

   }
   if(p>0)names(out) <- paste(rep(c("causeSpec","population"),ntab.strata),rep(names(tab.strata),each=2))
   else names(out) <- c("causeSpec","population")
   out$tau <- tau
   class(out) <- "cmp.rel"
   out
}




#' Plot the crude probability of death
#'
#' Plot method for cmp.rel. Plots the cumulative probability of death due to
#' disease and due to population reasons
#'
#' By default, the graph is plotted as a step function for the cause specific
#' mortality and as a piecewise linear function for the population mortality.
#' It is evaluated at all event and censoring times even though it constantly
#' changes also between these time points.
#'
#' If the argument \code{all.times} is set to \code{TRUE}, the plot is
#' evaluated at all times that were used for numerical integration in the
#' \code{cmp.rel} function (there, the default is set to daily intervals). If
#' only specific time points are to be added, this should be done via argument
#' \code{add.times} in \code{cmp.rel}.
#'
#' @param x a list, with each component representing one curve in the plot,
#' output of the function \code{cmp.rel}.
#' @param main the main title for the plot.
#' @param curvlab Curve labels for the plot. Default is \code{names(x)}, or if
#' that is missing, \code{1:nc}, where \code{nc} is the number of curves in
#' \code{x}.
#' @param ylim yaxis limits for plot.
#' @param xlim xaxis limits for plot (default is 0 to the largest time in any
#' of the curves).
#' @param wh if a vector of length 2, then the upper right coordinates of the
#' legend; otherwise the legend is placed in the upper right corner of the
#' plot.
#' @param xlab X axis label.
#' @param ylab y axis label.
#' @param lty vector of line types. Default \code{1:nc} (\code{nc} is the
#' number of curves in \code{x}). For color displays, \code{lty=1},
#' \code{color=1:nc}, might be more appropriate. If \code{length(lty)<nc}, then
#' \code{lty[1]} is used for all.
#' @param xscale Scale of the X axis. Default is in days (1).
#' @param col vector of colors. If \code{length(col)<nc}, then the
#' \code{col[1]} is used for all.
#' @param lwd vector of line widths. If \code{length(lwd)<nc}, then
#' \code{lwd[1]} is used for all.
#' @param curves Vector if integers, specifies which curves should be plotted.
#' May take values \code{1:nc}, where \code{nc} is the number of curves in
#' \code{x}. By default, all of the curves are plotted.
#' @param conf.int Vector if integers, specifies which confidence intervals
#' should be plotted. May take values \code{1:nc}, where \code{nc} is the
#' number of curves in \code{x}. By default, no confidence intervals are
#' plotted.
#' @param all.times By default, the disease specific mortality estimate is
#' plotted as a step function between event or censoring times. If set to
#' \code{TRUE}, the graph is evaluated at all estimated times.
#' @param ... additional arguments passed to the initial call of the plot
#' function.
#' @return No value is returned.
#' @seealso \code{rs.surv}
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' fit <- cmp.rel(Surv(time,cens)~sex,rmap=list(age=age*365.241),
#' 		ratetable=slopop,data=rdata,tau=3652.41)
#' plot(fit,col=c(1,1,2,2),xscale=365.241,conf.int=c(1,3))
#'
plot.cmp.rel <- function (x, main = " ", curvlab, ylim = c(0, 1), xlim, wh = 2,
    xlab = "Time (days)", ylab = "Probability", lty = 1:length(x), xscale=1,
    col = 1, lwd = par("lwd"), curves, conf.int, all.times=FALSE,...)
{
#wh= upper left coordinates of the legend, if of length 1, the legend is placed top left.
    tau <- x$tau
    x$tau <- NULL
    nc <- length(x)			#number of curves
    if (length(lty) < nc)		#if not enough different line types
        lty <- rep(lty[1], nc)
    else lty <- lty[1:nc]
    if (length(lwd) < nc)		#if not enough different line widths
        lwd <- rep(lwd[1], nc)
    else lwd <- lwd[1:nc]
    if (length(col) < nc)		#if not enough different colors
        col <- rep(col[1], nc)
    else col <- col[1:nc]
    if (missing(curvlab)) {		#if no curve labels desired
        if (mode(names(x)) == "NULL") { #and no curve labels prespecified
            curvlab <- as.character(1:nc)
        }
        else curvlab <- names(x)[1:nc]	#use prespecified if they exist
    }
    if (missing(xlim)) {		#if no limits desired
        xmax <- 0
        for (i in 1:nc) {
            xmax <- max(c(xmax, x[[i]][[1]]/xscale))	#take max time over all strata
        }
        xlim <- c(0, xmax)
    }
    if(all.times){
    for(it in 1:nc){
       x[[it]]$index <- 1:length(x[[it]][[1]])
    }
    }
    if(missing(curves))curves <- 1:nc
    if(missing(conf.int))conf.int <- NULL
    curves <- unique(curves)
    conf.int <- unique(conf.int)
    if(any((curves %in% 1:nc)==FALSE)) stop(paste("The curves argument should be specified as a vector of integers from 1 to", nc,sep=" "))
    if(any((conf.int %in% 1:nc)==FALSE)) stop(paste("The conf.int argument should be specified as a vector of integers from 1 to", nc,sep=" "))
    if(any((conf.int %in% curves)==FALSE)) stop("Confidence interval may only be plotted if the curve is plotted, see argument curves")

    col_nums <- floor(seq(from=95,to=50,length.out=length(conf.int)+2))
    col.conf.temp <-  sapply(col_nums,function(x)paste("gray",as.character(x),sep=""))
    col.conf.int <- rep("white",nc)
    col.conf.int[conf.int] <-  col.conf.temp[-c(1,length(col.conf.temp))]

    plot((x[[1]][[1]]/xscale)[x[[1]]$index], (x[[1]][[2]])[x[[1]]$index], type = "n", ylim = ylim, xlim = xlim,
        main = main, xlab = xlab, ylab = ylab, bty = "l", ...)		#plot estimates [[1]]=time, [[2]]=est
    if (length(wh) != 2) {
        wh <- c(xlim[1], ylim[2])
    }
    u <- list(...)
    if (length(u) > 0) {
        i <- pmatch(names(u), names(formals(legend)), 0)
        do.call("legend", c(list(x = wh[1], y = wh[2], legend = curvlab[curves],
            col = col[curves], lty = lty[curves], lwd = lwd[curves], bty = "n", bg = -999999),
            u[i > 0]))
    }
    else {
        do.call("legend", list(x = wh[1], y = wh[2], legend = curvlab[curves],
            col = col[curves], lty = lty[curves], lwd = lwd[curves], bty = "n", bg = -999999))
    }
    for(i in conf.int){
    	 if(i%%2==0)with(x[[i]],polygon(c(time[index][!is.na(lower[index])],rev(time[index][!is.na(upper[index])]))/xscale,c(lower[index][!is.na(lower[index])],rev(upper[index][!is.na(upper[index])])),col = col.conf.int[i] , border = FALSE))
    	 else with(x[[i]],my.poly(time[index][!is.na(lower[index])]/xscale,time[index][!is.na(upper[index])]/xscale,lower[index][!is.na(lower[index])],upper[index][!is.na(upper[index])],col = col.conf.int[i] , border = FALSE))
    }
    for (i in curves) {
        tip <- "s"
        if(i%%2==0)tip <- "l"
        lines((x[[i]][[1]]/xscale)[x[[i]]$index], (x[[i]][[2]])[x[[i]]$index], lty = lty[i], col = col[i],
            lwd = lwd[i], type=tip, ...)
    }
}

my.poly <- function(x1,x2,y1,y2,...){
	x1 <- rep(x1,each=2)[-1]
	y1 <- rep(y1,each=2)[-(2*length(y1))]
	x2 <- rep(x2,each=2)[-1]
	y2 <- rep(y2,each=2)[-(2*length(y2))]
	polygon(c(x1,rev(x2)),c(y1,rev(y2)),...)
}


print.cmp.rel <- function (x, ntp = 4, maxtime,scale=365.241, ...)
{
    tau <- x$tau
    x$tau <- NULL
    nc <- length(x)

    if (missing(maxtime)) {
        maxtime <- 0
        for (i in 1:nc) maxtime <- max(maxtime, x[[i]]$time)
    }
    tp <- pretty(c(0, maxtime/scale), ntp + 1)
    tp <- tp[-c(1, length(tp))]

    if(length(x[[1]]$add.times)>0 & length(x[[1]]$add.times)<5){
	tp <- sort(unique(c(tp,round(x[[1]]$add.times/scale,1))))
    }
    cat("Estimates, variances and area under the curves:\n")
    x$tau <- tau
    print(summary(x, tp,scale,area=TRUE), ...)
    invisible()
}



#' Summary of the crude probability of death
#'
#' Returns a list containing the estimated values at required times.
#'
#' The variance is calculated using numerical integration. If the required time
#' is not a time at which the value was estimated, the value at the last time
#' before it is reported. The density of the time points is set by the
#' \code{precision} argument in the \code{cmp.rel} function.
#'
#' @param object output of the function \code{cmp.rel}.
#' @param times the times at which the output is required.
#' @param scale The time scale in which the times are specified. The default
#' value is \code{1}, i.e. days.
#' @param area Should area under the curves at time \code{tau} be printed out?
#' Default is \code{FALSE}.
#' @param ... Additional arguments, currently not implemented
#' @return A list of values is returned.
#' @seealso \code{cmp.rel}
#' @keywords survival
#' @examples
#'
#' data(slopop)
#' data(rdata)
#' #calculate the crude probability of death and summarize it
#' fit <- cmp.rel(Surv(time,cens)~sex,rmap=list(age=age*365),
#'       ratetable=slopop,data=rdata,tau=3652.41)
#' summary(fit,c(1,3),scale=365.241)
#'
summary.cmp.rel <- function (object, times,scale=365.241,area=FALSE,...)
{
    tau <- object$tau
    object$tau <- NULL
    ng <- length(object)
    times <- sort(unique(times))*scale
    nt <- length(times)
    storage.mode(times) <- "double"
    storage.mode(nt) <- "integer"
    ind <- matrix(0, ncol = nt, nrow = ng)
    oute <- matrix(NA, ncol = nt, nrow = ng)
    outv <- oute
    outa <- matrix(NA,ncol=1,nrow=ng)
    storage.mode(ind) <- "integer"
    slct <- rep(TRUE, ng)
    for (i in 1:ng) {
        if (is.null((object[[i]])$est)) {
            slct[i] <- FALSE
        }
        else {
            z <- rep(NA,nt)
            for(kt in 1:nt)z[kt] <- rev(which(object[[i]][[1]]<=times[kt]))[1]
            ind[i, ] <- z
            oute[i, ind[i, ] > 0] <- object[[i]][[2]][z]
            outa[i,] <- object[[i]][[6]]
            if (length(object[[i]]) > 2)
                outv[i, ind[i, ] > 0] <- object[[i]][[3]][z]
        }
    }
    dimnames(oute) <- list(names(object)[1:ng], as.character(times/scale))
    dimnames(outv) <- dimnames(oute)
    rownames(outa) <- rownames(oute)
    colnames(outa) <- paste("Area at tau =",tau/scale)
    if(area)list(est = oute[slct, , drop = FALSE], var = outv[slct, ,
        drop = FALSE], area=outa[slct,,drop=FALSE])
    else list(est = oute[slct, , drop = FALSE], var = outv[slct, ,
        drop = FALSE])
}
