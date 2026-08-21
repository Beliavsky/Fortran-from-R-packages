#' Plots Results of a Monte Carlo Simulation
#'
#' Plots the empirical cumulative distribution function of a \samp{mcnode} or a \samp{mc}
#' object (\samp{"0"} and \samp{"V"} nodes) or the empirical CDF of the estimate of a
#' \samp{mcnode} or \samp{mc} object (\samp{"U"} and \samp{"VU"} nodes).
#'
#' For \samp{"VU"} and \samp{"U"} nodes, quantiles are calculated using
#' \code{\link{quantile.mc}} within each of the \samp{nsu} simulations (by columns). The
#' medians (or means if \samp{stat="mean"}) are plotted, with \samp{lim} quantiles as
#' the envelope.
#'
#' @references Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment.
#'   Plenum Press, USA, pp. 81-155.
#'
#' @param x a \samp{mcnode} or a \samp{mc} object.
#' @param prec the precision of the plot. \samp{0.001} gives an ecdf from
#'   the 0.000, 0.001, ..., 1.000 quantiles.
#' @param stat the function used for estimates (2D \samp{mc} or \samp{mcnode}). By default the median.
#' @param lim a vector of numbers (between 0 and 1) indicating the envelope. May be \samp{NULL} or empty.
#' @param na.rm should NA values be discarded?
#' @param griddim a vector of two integers for the grid size. If \samp{NULL}, calculated automatically.
#' @param xlab vector of labels for the x-axis. If \samp{NULL}, use the name of the node.
#' @param ylab vector of labels for the y-axis.
#' @param main vector of main titles of the graph.
#' @param draw should the plot be drawn?
#' @param paint should the envelopes be filled?
#' @param xlim x coordinate range: a vector of length 2 or a list of such vectors.
#' @param ylim y coordinate range: a vector of length 2 or a list of such vectors.
#' @param ... further arguments to be passed to \samp{plot.stepfun}.
#' @return a \samp{plot.mc} object (list of quantiles used to draw the plot), invisibly.
#' @seealso \code{\link{ecdf}}, \code{\link{plot}}, \code{\link{quantile.mc}},
#'   \code{\link{ggplotmc}} for a ggplot2 version.
#' @keywords hplot
#' @name plot.mc
#' @examples
#' data(total)
#' plot(xVUM3)
#' ## Only one envelope for 0.025 and 0.975
#' plot(xVUM3, lim=c(0.025, 0.975))
#' @export
plot.mc <- function(x, prec=0.001, stat = c("median","mean"), lim = c(0.025, 0.25, 0.75, 0.975), na.rm=TRUE, griddim = NULL, xlab = NULL, ylab = "Fn(x)", main = "", draw = TRUE, paint=TRUE, xlim=NULL,ylim=NULL,...)
{
  if(inherits(x,"mc")==TRUE) {
    x <- quantile.mc(x, probs=seq(0,1,prec),lim = lim, na.rm=na.rm, lnames=xlab)
    }

  if(draw) {
	y <- x                           # for a correct return
    stat <- match.arg(stat)

	 beau <- function(n){
		nc <- round(sqrt(n))
		nr <- ceiling(n/nc)
		c(nc,nr)}

   noms <- names(rapply(x,function(x) 1))    #moche mais efficace
   if(is.null(xlab)) xlab <- noms
   n <- length(noms)

   if(!is.null(ylim) & ((is.list(ylim) & length(ylim)!= n)|(is.vector(ylim) & length(ylim)!= 2))) stop("ylim should be NULL, a vector of 2 elements or a list of length the number of nodes")
   if(!is.null(xlim) & ((is.list(xlim) & length(xlim)!= n)|(is.vector(xlim) & length(xlim)!= 2))) stop("xlim should be NULL, a vector of 2 elements or a list of length the number of nodes")

	 main <- rep(main,n)
	 xlab <- rep(xlab,n)
	 ylab <- rep(ylab,n)

  if(is.null(griddim)) griddim <- beau(n)
  if(prod(griddim) < n) op <- par(mfrow=griddim,ask=TRUE,mar=c(5,4,.2,.2))
     else op <- par(mfrow=griddim, mar=c(5,4,.2,.2))

  try({   #to restore par in case of error

  i <- 1
  env <- environment()

  LEPLOT <- function(y,...){
      if(nrow(y) != 1) {
        if(stat=="median") y <- y[-2,,drop=FALSE]
        else y <- y[-1,,drop=FALSE]}                                              #Retrait median or mean
  		nr <- nrow(y)
      i <- get("i",envir=env)
      xlima <- if(is.null(xlim)) range(y,na.rm=TRUE) else
		xlima <- if(is.list(xlim)) xlim[[i]] else xlim
  	  if(xlima[1]==xlima[2]) xlima <- xlima + c(-0.5,0.5)
      ylima <- if(is.null(ylim)) c(0,1) else
		ylima <- if(is.list(ylim)) ylim[[i]] else ylim
      x50 <- plot.stepfun(y[1,], main=main[i], xlim=xlima, ylim=ylima, ylab = ylab[i], verticals = TRUE, do.points = FALSE, xlab=xlab[i], lwd=3, ...)
      abline(h = c(0, 1), col =  "gray70", lty = 3)

      # Points for the polygon used to fill the envelope
      if(paint){
        ti.l <- x50$t[-length(x50$t)]
        ti.r <- x50$t[-1L]
        y50 <- x50$y
        thex50 <- rev(as.vector(rbind(ti.l,ti.r)))
        they50 <- rev(as.vector(rbind(y50, y50)))
      }

      if(nr > 1){
        rankplot <- 1 + order(-abs(lim-0.5))      # in order to draw over in the good order
        for(j in rankplot) {
          xp <- plot.stepfun(y[j,], verticals=TRUE, do.points=FALSE, add= TRUE, lty=3 ,col="gray30",...)
          if(paint){
            ti.lp <- xp$t[-length(xp$t)]
            ti.rp <- xp$t[-1L]
            yp <- xp$y
            thexp <- as.vector(rbind(ti.lp,ti.rp))
            theyp <- as.vector(rbind(yp, yp))
            color <- grey(abs(lim[j-1]-.5)+.25)
            polygon(c(thexp,thex50), c(theyp,they50), col= color)
          }
        }
      }
    assign("i",i+1,envir=env) }

  rapply(y,LEPLOT)

  })
	par(op)
  }
  class(x) <- "plotmc"
  return(invisible(x))}

#' @rdname plot.mc
#' @export
plot.mcnode <- function(x, ...)
{
  nom <- deparse(substitute(x), width.cutoff = 500L, nlines=1)
  x <- list(x)
  names(x) <- nom
  class(x) <- "mc"
  x <- plot.mc(x, ... )
  return(invisible(x))}

#' @rdname plot.mc
#' @export
plot.plotmc <- function(x, ...)
{
  x <- plot.mc(x, ... )
  return(invisible(x))}
