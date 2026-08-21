#' @title Spaghetti Plot of mc/mcnode Object
#' @description Use plot to draw spaghetti plots for the mc/mcnode objects.
#' @param x mc/mcnode object  
#' @param griddim a vector of two integers, indicating the size of the grid of the graph. If NULL, the grid is calculated to produce a "nice" graph.
#' @param xlab vector of labels for the x-axis. If NULL, use the name of the node.
#' @param ylab vector of labels for the y-axis.
#' @param main vector of main titles of the graph.
#' @param maxlines the maximum number of ecdf to draw.
#' @param ... further arguments to be passed to plot.stepfun()
#' @aliases spaghetti.mc
#' @aliases spaghetti.mcnode
#' @examples 
#' data(total)
#' spaghetti(mc(xVUM))
#' spaghetti(xVUM)
#' 
#' @rdname spaghetti
#' @export
spaghetti <- function(x,...){
  UseMethod("spaghetti")
}

#' @rdname spaghetti
#' @method spaghetti mc
#' @export
spaghetti.mc <- function(x, 
                         griddim = NULL, 
                         xlab = names(x),
                         ylab = "F(n)", 
                         main = "",  
                         maxlines = 100,
                         ...)
{
  beau <- function(n){
    nc <- round(sqrt(n))
    nr <- ceiling(n/nc)
    c(nc,nr)
  }
  
  l <- length(x)
  main <- rep(main,l)
  xlab <- rep(xlab,l)
  ylab <- rep(ylab,l)
  
  
  loutm <- lapply(x,attr,which="outm")
  dimm <- sapply(x,dim)
  n <- sum(dimm[3,]* (loutm=="each") + (loutm!="each" & loutm!="none"))
  
  
  if(is.null(griddim)) griddim <- beau(n)
  if(prod(griddim) < n) {
    op <- par(mfrow=griddim,ask=TRUE,mar=c(5,4,.2,.2))
  } else {
    op <- par(mfrow=griddim,mar=c(5,4,.2,.2))
  }
  
  
  try({  #to restore par in case of error
    
    for(i in 1:l){
      
      if(is.null(loutm[[i]])) loutm[[i]] <- "each"
      if(loutm[[i]][1] == "none") next                                             # Pass outm == none
      for(j in loutm[[i]]){                                                     # loop on the nb of stat, j is the name of the function
        
        if(j == "each"){
          nvar <- dim(x[[i]])[3]
          if(nvar==1) xlab2 <- xlab[i] else xlab2 <- paste(xlab[i],1:nvar,sep="")
        } else {
          func <- get(j,mode="function")                                        # apply the function
          x[[i]] <- apply(x[[i]],c(1,2),func)
          dim(x[[i]]) <- c(dim(x[[i]]),1)
          nvar <- 1                                                             #1 dimension now for this stat
          xlab2 <- paste(j,xlab[i])                                             #change the name with the name of the stat
        }
        
        if(is.logical(x[[i]])) x[[i]] <- unclass(x[[i]]) * 1                        # unclass to avoid Ops.mcnode
        
        for(k in 1:nvar){
          
          value <- as.matrix(unmc(x[[i]][,,k])) # nvariates data
          variable <- rep(1:dim(value)[2], each=dim(value)[1])
          ncol.value <- min(maxlines, dim(value)[2])
          
          dim(value) <- c(prod(dim(value)),1)
          
          
          data.df <- data.frame(variable=as.factor(variable),value=value)
          data.df <- data.df[which(!is.na(data.df$value)), ]
          
          plot.stepfun(data.df[which(data.df$variable==1),2], 
                       xlim=c(0,max(data.df[,2])),
                       main = main[i],
                       xlab = xlab2[k],
                       ylab = ylab[i],
                       do.points = FALSE)
          if(ncol.value>1){
            for (p in 2: ncol.value){
              plot.stepfun(data.df[which(data.df$variable==p),2], do.points = FALSE, add=T)}}}
      }
    }
  })
  par(op)
  return(invisible())
}

#' @rdname spaghetti
#' @method spaghetti mcnode
#' @export
spaghetti.mcnode <- function(x, ...)
{ 
  nom <- deparse(substitute(x), width.cutoff = 500L, nlines=1)
  x <- list(x)
  names(x) <- nom
  spaghetti.mc(x, ...)
}


