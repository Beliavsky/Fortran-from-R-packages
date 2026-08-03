#' Plotting HMM expected volatility for SPX overlaid with adjusted VIX
#' 
#' This utility plots the HMM expected volatility of SPX overlaid with the VIX index adjusted by a ratio.
#' The expected volatility is shown to have a long-term ratio of 0.79 relative to the VIX index. This plot will
#' show how HMM deviates from VIX in a shorter time window.
#' Optionally the insert shows the relation between the return and volatility indicated by each state. 
#' This plot is also called "volatility yield curve".
#'
#' @param object an ldhmm object with a stationary solution. If this is set to \code{NULL},
#'               an internal 10-state HMM object will be used.
#' @param days.pa a positive integer specifying trading days per year, default is 252.
#' @param start.date Date or character of ISO format (YYYY-MM-DD), 
#'                   specifying the start date of the plot, default is \code{NULL}, which is converted to 1.5 years ago.
#' @param end.date Date or character of ISO format (YYYY-MM-DD), 
#'                 specifying the end date of the plot, default is \code{NULL}, which means the latest date.
#' 
#' @param px.scale numeric, specifying the scaling factor when plotting price trend, default is 15.
#'                 The closing price is converted to cumulative return by the price of the first date.
#'                 Then plot from the mid-point of volatility axis with this scale.
#' @param px.origin numeric, specifying the starting value of the index price line,
#'                  the default is \code{NULL}, which will start the index price line from the middle of y-axis.
#'
#' @param vix.adj.ratio numeric, if specified, VIX index is adjusted and plotted, default is \code{NULL}.
#'                      Default is to use the long-term ratio between VIX and 10-state HMM, which is about 0.79.
#' 
#' @param insert.plot logical, if true, also plot the volatility-return as insert in upper-right corner, default is \code{TRUE}.
#' @param insert.viewport optional viewport for the insert, default is \code{NULL}, 
#'                        which is internally set to \code{grid::viewport(.8, .75, .3, .3)}.
#'                        
#' @keywords VIX
#'
#' @author Stephen H. Lihn
#' 
#' @export
#' 
#' @importFrom ggplot2 qplot
#' @importFrom grid viewport
#' @importFrom scales alpha
#' 
#' @examples
#' \dontrun{
#'     ldhmm.plot_spx_vix_obs(h)
#' }
### <======================================================================>
ldhmm.plot_spx_vix_obs <- function(object, days.pa=252, 
                                      start.date=NULL, end.date=NULL,
                                      px.origin=NULL,
                                      px.scale=NULL,
                                      vix.adj.ratio=NULL,
                                      insert.plot=TRUE, 
                                      insert.viewport=NULL)
{
    if (is.null(object)) object <- ldhmm.read_sample_object()
    
    if (is.null(start.date)) start.date <- Sys.Date()-365*1.5

    if (!is.null(start.date)) {
        if (is(start.date, "character")) start.date <- as.Date(start.date)
        if (! is(start.date, "Date")) stop("start.date must be a Date object")
    }
    if (!is.null(end.date)) {
        if (is(end.date, "character")) end.date <- as.Date(end.date)
        if (! is(end.date, "Date")) stop("end.date must be a Date object")
    }
    
    spx <- ldhmm.ts_log_rtn(on="days", end.date=end.date, fred.data=TRUE)
    vix <- ldhmm.ts_log_rtn("vix", on="days", end.date=end.date, fred.data=TRUE)
    hss <- ldhmm.decoding(object, spx$x)
    v <- ldhmm.decode_stats_history(hss, annualize=TRUE)[,"V"]
    
    # calculate the ratio between HMM vol and VIX
    get_vol_ratio <- function(d) {
        i <- which(spx$d==d)
        j <- which(vix$d==d)
        if (length(i)>0) v[i]/vix$p[j] else NaN
    }
    if (is.null(vix.adj.ratio)) {
        rr <- sapply(vix$d, get_vol_ratio)
        vix.adj.ratio <- mean(rr, na.rm=TRUE)
        # this 0.8, about the same as sqrt(252/365)
    }
    #
    # --------------------------------------------------
    max_date <- max(spx$d)
    min_date <- min(spx$d)
    if (!is.null(start.date)) {
        if (start.date > min_date) min_date <- start.date
    } 
    if (!is.null(end.date)) {
        if (end.date < max_date) max_date <- end.date
    } 
    I <- which(spx$d >= min_date & spx$d <= max_date)
    J <- which(vix$d >= min_date & vix$d <= max_date)
    
    # --------------------------------------------------
    vx <- vix$p*vix.adj.ratio
    y_max <- max(c(v[I], vx[J]))*1.5 # 1.5 to extend the top half empty for prices
    y_min <- min(c(v[I], vx[J]))*0.9
    
    # spx vol
    plot(spx$d[I], v[I], type="l", col="black",
         ylim=c(y_min, y_max), 
         xlab=sprintf("Date (%s to %s)", min_date, max_date), 
         ylab="V (volatility)",
         main="Expected Volatility from 10-State HMM vs Adjusted VIX")
    # vix
    lines(vix$d[J], vx[J], col=scales::alpha("magenta",0.7), lwd=2)
    abline(h=0, lty=1)
    
    # index px
    spxr <- spx$p[I]
    spxr <- log(spxr/spxr[1])
    # plot prices on the upper 1/3 of area
    if (is.null(px.origin)) px.origin <- y_max-(y_max-y_min)/3
    if (is.null(px.scale)) px.scale <- (y_max-px.origin)/max(abs(spxr))*0.5
    lines(spx$d[I], px.origin+spxr*px.scale, col="blue")
    
    ld_stats <- ldhmm.ld_stats(hss)
    for (i in 1:hss@m) {
        abline(h=ld_stats[i,2]*sqrt(days.pa)*100, lty=2, col="blue")
    }
    # these black dots should always on top of other plots
    points(spx$d[I], v[I], col="black", cex=0.4) 
    
    text(max_date, px.origin*1.05, pos=2, cex=0.8, col="blue",
         label="SPX (rescaled)")
    text(max_date, y_max*0.95, pos=2, cex=0.8, col="red",
         label=sprintf("VIX Adj. Ratio = %.2f", vix.adj.ratio))
    
    # insert
    # --------------------------------------------------
    # plot the volatility/return insert
    if (insert.plot) {
        if (is.null(insert.viewport)) insert.viewport <- grid::viewport(.23, .77, .23, .23)
        R  <- ld_stats[,1]*days.pa
        V <- ld_stats[,2]*sqrt(days.pa)*100
        qp  <- ggplot2::qplot(V, R, data=data.frame(R=R, V=V))
        print(qp, vp=insert.viewport)
    }
}
### <---------------------------------------------------------------------->
