#' Draw a Tornado Chart
#'
#' Draws a Tornado chart as provided by \code{\link{tornado}} or \code{\link{tornadounc}}.
#'
#' A point is drawn at the estimate, and the segment reflects the uncertainty around it.
#'
#' @param x a \code{\link{tornado}} or \code{\link{tornadounc}} object.
#' @param which which output to plot (for multivariate output).
#' @param name vector of names for input variables. If \samp{NULL}, taken from the object.
#' @param stat the statistic of the output to be considered. For a \samp{tornado} object:
#'   \samp{"median"} or \samp{"mean"}. For a \samp{tornadounc} object: should match a row
#'   name, or use the row number.
#' @param xlab label of the x axis. If \samp{"method"}, use the correlation method.
#' @param ylab label of the y axis.
#' @param ... further arguments to be passed to the \samp{plot} function.
#' @return \samp{NULL} (called for side effects).
#' @seealso \code{\link{tornado}}
#' @keywords hplot
#' @name plot.tornado
#' @examples
#' data(ec)
#' x <- evalmcmod(ec$modEC2, nsv=100, nsu=100, seed=666)
#' tor <- tornado(x, 7)
#' plot(tor)
#' @export
plot.tornado <- function(x,which=1,name=NULL,stat=c("median","mean"),xlab="method",ylab="",...)
{
  val <- x$value[[which]]
  if(is.null(val)) stop("Invalid value for which")
  nc <- ncol(val)
  nr <- nrow(val)
  if(!is.null(name)) {colnames(val) <- (rep(name,nc))[1:nc]}

  if(xlab=="method") xlab <- c("Spearman's rho statistic","Kendall's tau statistic","Pearson correlation")[pmatch(x$method,c("spearman","kendall","pearson"))]

	plot(c(-1.5,1),c(1,nc),type="n",axes= FALSE, xlab=xlab,ylab=ylab,...)
  axis(1,at=c(-1,-0.5,0,0.5,1))

  stat <- match.arg(stat)
  stat <- ifelse(stat=="mean" && nr!=1, 2 ,1)
  val <- val[,order(abs(val[stat,])),drop=FALSE]
	if(nr==1){
		segments(0,1:nc,val,1:nc,lwd=2)
    points(val,1:nc, pch="|", lwd=2)
    }

  else {
		segments(0,1:nc,val[stat,],1:nc,lwd=2,col="grey")
    points(val[stat,],1:nc, pch="|", lwd=2)
    if(nr>3){
        val <- apply(val[3:nr,],2,range)
        segments(val[1,],1:nc,val[2,],1:nc,lwd=2)
        points(val[1,],1:nc, pch="|", lwd=2)
        points(val[2,],1:nc, pch="|", lwd=2)
            }
    }

		abline(v=0)
		text(-1.4,1:nc,labels=paste(colnames(val),sep=":"),cex=.8)
		}

#' @rdname plot.tornado
#' @export
plot.tornadounc <- function(x,which=1, stat="median", name=NULL, xlab="method", ylab="",...)
{
  statposs <- rownames(x$value[[which]])

  if(is.character(stat)) stat <- pmatch(stat, rownames(x$value[[which]]))
  if(is.na(stat)) stop("stat should match with: ",paste(statposs,collapse=", "))

  x$value <- list(x$value[[which]][stat,,drop=FALSE])
	plot.tornado(x,which=1, stat="median", name=name, xlab=xlab,ylab=ylab,...)
 }
