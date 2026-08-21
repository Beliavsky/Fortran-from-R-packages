#' @title Histogram of a Monte Carlo Simulation (ggplot version)
#' @description Shows histogram of a \samp{mcnode} or a \samp{mc} object by ggplot framework.
#' @param x  an `mc` or an `mcnode` object  
#' @param griddim A vector of two integers, indicating the size of the grid of the graph. 
#' If `NULL`, the grid is calculated to produce a "nice" graph.
#' @param xlab Vector of labels for the x-axis. If `NULL`, use the name of the node.
#' @param ylab Vector of labels for the y-axis.
#' @param main Vector of main titles of the graph
#' @param bins Number of bins. Defaults to 30.
#' @param which An argument used for a multivariate `mcnode`. 
#' Can specify which variate plot to display. 
#' When variates are more than one, the output will be saved in a plot list by 
#' default or use the number of which variate to display.
#' @param ... Further arguments to be passed to geom_histogram()
#' 
#' @aliases gghist.mc
#' @aliases gghist.mcnode
#'
#' @returns a ggplot object.
#' 
#' @author Yu Chen and Regis Pouillot
#' 
#' @seealso [hist.mc()]
#' 
#' @examples 
#' data(total)
#' # When mcnode has one variate
#' gghist(xV)
#' # When mcnode has two variates, the two plots will be saved in a list 
#' # if affected to a variable
#' gplots <- gghist(xVUM) 
#' # show the first variate plot of xVUM mcnode
#' gplots[[1]] 
#' # directly show the first variate plot of xVUM mcnode
#' gghist(xVUM, which = 1) #directly show the first variate plot of xVUM mcnode
#' # Post process
#' gplots[[1]] + ggplot2::geom_histogram(color = "red",fill="blue")
#' 
#' @rdname gghist
#' @export
gghist <- function(x,...){
  UseMethod("gghist")
}

#' @rdname gghist
#' @method gghist mcnode
#' @export
gghist.mcnode <- function(x, griddim = NULL, xlab = names(x),ylab = "Frequency", main = "", bins=30, which=NULL,...)
{
  
  nom <- deparse(substitute(x), width.cutoff = 500L, nlines=1)
  x <- list(x)
  names(x) <- nom
  # the function beau calculate a nice grid
  beau <- function(n){
    nr <- round(sqrt(n))          # plots layout
    nc <- round(sqrt(n))
    while (nr*nc<n) {
      if ((nr-nc)==0){
        nc <- nc+1
      } else {
        nr <- nr+1
      }
    }
    c(nr,nc)}
  
  l <- length(x)
  main <- rep(main,l)
  xlab <- rep(xlab,l)
  ylab <- rep(ylab,l)
  
  
  loutm <- lapply(x,attr,which="outm")
  dimm <- sapply(x,dim)
  n <- sum(dimm[3,]* (loutm=="each") + (loutm!="each" & loutm!="none"))
  
  
  if(is.null(griddim)) griddim <- beau(n)
  
  try({  #to restore par in case of error
    plotlist <- list()
    
    for(i in 1:l){
      
      if(is.null(loutm[[i]])) loutm[[i]] <- "each"
      if(loutm[[i]][1] == "none") next                                             # Pass outm == none
      
      
      for(j in loutm[[i]]){                                                     # loop on the nb of stat, j is the name of the function
        
        if(j == "each"){
          nvar <- dim(x[[i]])[3]
          if(nvar==1) xlab2 <- xlab[i] else xlab2 <- paste(xlab[i],1:nvar,sep="")
        }else{
          func <- get(j,mode="function")                                        # apply the function
          x[[i]] <- apply(x[[i]],c(1,2),func)
          dim(x[[i]]) <- c(dim(x[[i]]),1)
          nvar <- 1                                                             #1 dimension now for this stat
          xlab2 <- paste(j,xlab[i])                                             #change the name with the name of the stat
        }
        
        if(is.logical(x[[i]])) x[[i]] <- unclass(x[[i]]) * 1                        # unclass to avoid Ops.mcnode
        
        for(k in 1:nvar){
          df <- data.frame(x = as.vector(x[[i]][,,k])) # nvariates data
          plotlist[[k]] <- ggplot2::ggplot(df,aes(x=x))+
            ggplot2::geom_histogram(bins = bins,fill="grey",colour="black",...)+
            ggplot2::xlab(xlab2[k])+
            ggplot2::ylab(ylab[i])+
            ggplot2::ggtitle(main[i])+
            ggplot2::theme_bw()+
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                  panel.grid.minor = ggplot2::element_blank(),
                  plot.title = ggplot2::element_text(hjust = 0.5))
          
        }
      }
    }
    
    if(nvar>1){
      if(is.null(which)){
        cat("There are more than one variate in this mcnode. Output is a plotlist.\n")
        cat("Can achieve single variate plot by use list[[n]]")
        return(invisible(plotlist))
      } else {return(plotlist[[which]])}
    } else {return(plotlist[[1]])}
  })
}


#' @rdname gghist
#' @method gghist mc
#' @export
gghist.mc <- function(x, griddim = NULL, xlab = names(x),ylab = "Frequency", main = "", bins=30, ...)
{
  if(!inherits(x,"mc")) stop("Please input mc object data")
  # the function beau calculate a nice grid
  beau <- function(n){
    nr <- round(sqrt(n))          # plots layout
    nc <- round(sqrt(n))
    while (nr*nc<n) {
      if ((nr-nc)==0){
        nc <- nc+1
      } else {
        nr <- nr+1
      }
    }
    c(nr,nc)}
  
  l <- length(x)
  main <- rep(main,l)
  xlab <- rep(xlab,l)
  ylab <- rep(ylab,l)
  
  
  loutm <- lapply(x,attr,which="outm")
  dimm <- sapply(x,dim)
  n <- sum(dimm[3,]* (loutm=="each") + (loutm!="each" & loutm!="none"))
  
  
  if(is.null(griddim)) griddim <- beau(n)
  
  try({  #to restore par in case of error
    plotlist <- list()
    
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
        
        if(is.logical(x[[i]])) x[[i]] <- unclass(x[[i]]) * 1                    # unclass to avoid Ops.mcnode
        
        for(k in 1:nvar){
          df <- data.frame(x = as.vector(x[[i]][,,k])) # nvariates data
          ggp <- ggplot2::ggplot(df, ggplot2::aes(x=x))+
            ggplot2::geom_histogram(bins = bins,fill="grey",colour="black", ...)+
            ggplot2::xlab(xlab2[k])+
            ggplot2::ylab(ylab[i])+
            ggplot2::ggtitle(main[i])+
            ggplot2::theme_bw()+
            ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                  panel.grid.minor = ggplot2::element_blank(),
                  plot.title = ggplot2::element_text(hjust = 0.5))
          plotlist <- c(plotlist, list(ggp))
        }
      }
    }
    gga <- ggpubr::ggarrange(plotlist = plotlist,ncol=griddim[2],nrow=griddim[1])
    
  })
  return(gga)
}




