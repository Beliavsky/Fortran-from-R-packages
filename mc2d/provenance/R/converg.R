#' Graph of Running Statistics for Convergence Diagnostics
#'
#' Provides basic graphs to evaluate the convergence of a node of a \code{\link{mc}} or
#' a \code{\link{mccut}} object in the variability or in the uncertainty dimension.
#'
#' If the node is of type \samp{"V"}, the running mean, median and \samp{probs} quantiles
#' according to the variability dimension will be provided. If the node is of type \samp{"VU"}
#' and \samp{margin="var"}, this graph will be provided on one simulation in the uncertainty
#' dimension (chosen by \samp{iter}).
#'
#' If the node is of type \samp{"U"}, the running mean, median and \samp{lim} quantiles
#' according to the uncertainty dimension will be provided.
#'
#' If the node is of type \samp{"VU"} (with \samp{margin="unc"} or from a \samp{mccut} object),
#' one graph is provided for each of the mean, median and \samp{probs} quantiles calculated in
#' the variability dimension.
#'
#' @note This function may be used on a \samp{mccut} object only if a \samp{summary.mc}
#' function was used in the third block of the \code{\link{evalmccut}} call. The values
#' used as \samp{probs} in \samp{converg} should have been used in that \samp{summary.mc}.
#'
#' @param x a \code{\link{mcnode}}, \code{\link{mc}} or \code{\link{mccut}} object.
#' @param node the node to be considered in a \samp{mc} or \samp{mccut} object, either as
#'   order number or name. By default: the last node. Cannot be of type \samp{"0"}.
#' @param margin the margin used to plot the graph. Used only if the node is a \samp{"VU"} node.
#' @param nvariates the variate to be considered. Used only for multivariate nodes.
#' @param iter if \samp{margin=="var"} and the node is \samp{"VU"}, specifies which iteration
#'   in the uncertainty dimension to use.
#' @param probs the quantiles to be provided in the variability dimension.
#' @param lim the quantiles to be used in the uncertainty dimension.
#' @param griddim a vector of two integers for the grid size. If \samp{NULL}, calculated automatically.
#' @param log if \samp{TRUE}, the data will be log-transformed.
#' @return \samp{invisible(NULL)} (called for its side effect of drawing a plot).
#' @keywords hplot
#' @examples
#' data(total)
#' converg(xVU, margin="var")
#' converg(xVU, margin="unc")
#' @export
converg <- function(x, node=length(x), margin=c("var","unc"), nvariates=1, iter=1, probs=c(0.025,0.975), lim=c(0.025,0.975), griddim=NULL, log=FALSE)
{
  margin <- match.arg(margin)

  cmcc <- inherits(x,"mccut")
  cmc  <- inherits(x,"mc")
  if(!cmc && !inherits(x,"mcnode") && !cmcc) stop("converg function is not valid for this object")

  if(cmcc){     # selection of the summary
    x <- x[[which(sapply(x,inherits,what="summary.mccut"))[1]]]
    if(length(x)==0) stop("summary.mc was not evaluated in evalmcmodcut : impossible to produce a converg")
    }

	if(cmc || cmcc){   # selection of the node
     if(is.numeric(node)) node <- names(x)[node]
	   if(!(node %in% names(x)))	stop("node is not a valid value")
     x <- x[[node]]}  # selection of the node
  else node <-  deparse(substitute(x), width.cutoff = 500L, nlines=1) # for mcnode

cat("Convergence for ",node,", variates: ", nvariates,"\n",sep="")

  typen <- attr(x,"type")
  if(typen=="0") stop("converg is not valid for 0 mcnode")

  if(cmcc){   #selection of the quantiles
    if(typen=="V") stop("converg is not valid for V mcnode in a mccut object")
    if(is.list(x)) x <- x[[nvariates]]                                          # summary from multivariates nodes
    nvariates <- 1
    x <- aperm(x,c(3,2,1))
    quel <- c("mean",paste(c(50,probs*100),"%",sep=""))
    if(typen != "U") {
      if(!(all(quel %in% dimnames(x)[[1]]))) stop("impossible to provide those probs since they were not evaluated in the evalmcmodcut function")
      x <- x[quel,,1]
      }
    }

  if(log) x <- log(x)

  stat <- function(i,x,lesq){
      c(mean=mean(x[1:i],na.rm=TRUE),quantile(x[1:i],probs=c(0.5,lesq),na.rm=TRUE,names=FALSE))}

  graph <- function(x,typen,titre){
    if(typen=="U") {xlab <- "Iterations (Unc)" ; lty <- 2}
    else {xlab <- "Iterations (Var)" ; lty <- 1}
    ylim <- range(x)
    plot(x[1,],ylim=ylim,main=titre,xlab=xlab,ylab="",type="l",col="blue",lwd=2)
    lines(x[2,],col="red",lwd=2,lty=lty)
    for(i in 2:nrow(x)){
      lines(x[i,],col="red",lty=lty)
      }
    }

  if(typen=="VU" && !cmcc && margin=="var"){
      x <- x[,iter,nvariates,drop=FALSE]
      node <- paste(node,"[,",iter,"]",sep="")
      typen <- "V"
      cat("Iteration #",iter," in the uncertainty dimension\n",sep="")
      }

  if(typen %in% c("V","U")){
      x <- as.vector(x[,,nvariates])
      if(typen == "U") mesq <- lim else mesq <- probs
      x <- sapply(1:length(x),stat,x=x,lesq=mesq)
      graph(x,typen,node)
      }
  else {  # 2DVU with margin == "unc" and mccut
          if(!cmcc) x <- apply(x[,,nvariates,drop=FALSE],2,function(x) c(mean=mean(x,na.rm=TRUE),quantile(x,probs=c(0.5,probs),na.rm=TRUE)))
          ng <- nrow(x)
          nom <- paste(node,", stat:",rownames(x),sep="")
          typen <- "U"

          if(is.null(griddim)){
            nc <- round(sqrt(ng))
	          nr <- ceiling(ng/nc)
            griddim <- c (nc,nr)
          }
          if(prod(griddim) < ng) op <- par(mfrow=griddim,ask=TRUE) else op <- par(mfrow=griddim)

          cat("Each graph draws the evolution in the uncertainty dimension of a statistics calculated on the variability dimension\n")

          try({                 #to restore par in case of error
            for(i in 1:ng) {
              mesq <- lim
              tmp <- sapply(1:length(x[i,]),stat,x=x[i,],lesq=mesq)
              graph(tmp,"U",nom[i])
            }
          })
          par(op)
        }
nom <- paste(mesq*100,"%",sep="")
cat("Legend:\nthick blue = running mean\nthick red = running median\nthin red =",mesq,"running quantiles\n")
if(typen=="U") cat("... in the uncertainty dimension\n") else cat("... in the variability dimension\n")
return(invisible(NULL))}
