#' @title Spaghetti Plot of `mc` or `mcnode` Object
#' @description Use ggplot to draw spaghetti plots for the [mc] or [mcnode] objects.
#' @param x  an `mc` or an `mcnode` object  
#' @param griddim a vector of two integers, indicating the size of the grid of the graph. If `NULL`, the grid is calculated to produce a "nice" graph.
#' @param xlab vector of labels for the x-axis. If `NULL`, use the name of the node.
#' @param ylab vector of labels for the y-axis.
#' @param main vector of main titles of the graph
#' @param which An argument used for an `mcnode` with multivariates. 
#' Can specify which variate plot to display. When variates are more than one, 
#' the output will be saved in a plot list by default or use the number of which variate to display.
#' @param maxlines the maximum number of ecdf to draw.
#' @param ... further arguments to be passed to [ggplot2::stat_ecdf()]
#' @aliases ggspaghetti.mc
#' @aliases ggspaghetti.mcnode
#' @examples 
#' data(total)
#' # mcnode has one variate
#' 
#' @author Yu Chen and Regis Pouillot
#' 
#' @rdname ggspaghetti
#' @export
ggspaghetti <- function(x,...){
  UseMethod("ggspaghetti")
}

#' @rdname ggspaghetti
#' @method ggspaghetti mc
#' @export
ggspaghetti.mc <- function(x, 
                           griddim = NULL, 
                           xlab = names(x),
                           ylab = "F(n)", 
                           main = "",  
                           maxlines= 100,
                           ...)
{
  beau <- function(n){
    nc <- round(sqrt(n))
    nr <- ceiling(n/nc)
    c(nc,nr)}
  
  l <- length(x)
  
  main <- rep(main,l)
  xlab <- rep(xlab,l)
  ylab <- rep(ylab,l)
  
  loutm <- lapply(x,attr,which="outm")
  dimm <- sapply(x,dim)
  n <- sum(dimm[3,]* (loutm=="each") + (loutm!="each" & loutm!="none"))
  
  if(is.null(griddim)) griddim <- beau(n)
  
  try({ 
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
        
        if(is.logical(x[[i]])) x[[i]] <- unclass(x[[i]][,,k]) * 1                        # unclass to avoid Ops.mcnode
        
        for(k in 1:nvar){
          
          value <- as.matrix(unmc(x[[i]][,,k]) ) # nvariates data
          ncol.value <- pmin(dim(value)[2], maxlines)
          # Limite to the maxlines 
          value <- value[, 1:ncol.value, drop=FALSE]
          
          variable <- rep(1:dim(value)[2],each=dim(value)[1])
          
          
          dim(value) <- c(prod(dim(value)),1)
          data.df <- data.frame(variable=as.factor(variable),value=value)
          
          gg.temp <- ggplot2::ggplot(data.df,aes(value,col=variable))+
            ggplot2::scale_color_manual(values = rep("Black", ncol.value)) +
            ggplot2::stat_ecdf(geom = "step", pad = FALSE, show.legend = FALSE,...) +
            ggplot2::ggtitle(main[i])+
            ggplot2::theme_bw()+
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
            ggplot2::xlab(xlab2[k])+
            ggplot2::ylab(ylab[i])
        }
      }
      plotlist[[i]] <- gg.temp
    }
    gga <- ggpubr::ggarrange(plotlist = plotlist,ncol=griddim[2],nrow=griddim[1])
  })
  return(gga)
}

#' @rdname ggspaghetti
#' @method ggspaghetti mcnode
#' @export
ggspaghetti.mcnode <- function(x, 
                               griddim = NULL, 
                               xlab = names(x),
                               ylab = "F(n)", 
                               main = "",
                               which = NULL,
                               maxlines = 100,
                               ...)
{
  beau <- function(n){
    nc <- round(sqrt(n))
    nr <- ceiling(n/nc)
    c(nc,nr)}
  x <- mc(x)
  l <- length(x)
  main <- rep(main,l)
  xlab <- rep(xlab,l)
  ylab <- rep(ylab,l)
  
  loutm <- lapply(x,attr,which="outm")
  dimm <- sapply(x,dim)
  n <- sum(dimm[3,]* (loutm=="each") + (loutm!="each" & loutm!="none"))
  
  if(is.null(griddim)) griddim <- beau(n)
  
  try({ 
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
        
        if(is.logical(x[[i]])) x[[i]] <- unclass(x[[i]][,,k]) * 1               # unclass to avoid Ops.mcnode
        
        for(k in 1:nvar){
          
          value <- as.matrix(unmc(x[[i]][,,k]) ) # nvariates data
          ncol.value <- pmin(dim(value)[2],maxlines)
          
          value <- value[,1:ncol.value,drop=FALSE]
          
          variable <- rep(1:dim(value)[2],each=dim(value)[1])
          
          
          dim(value) <- c(prod(dim(value)),1)
          data.df <- data.frame(variable=as.factor(variable),value=value)
          
          gg.temp <- ggplot2::ggplot(data.df, ggplot2::aes(value,col=variable))+
            ggplot2::scale_color_manual(values = rep("Black", ncol.value)) +
            ggplot2::stat_ecdf(geom = "step", pad = FALSE, show.legend = FALSE,...) +
            ggplot2::ggtitle(main[i])+
            ggplot2::theme_bw()+
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                  panel.grid.major = ggplot2::element_blank(),
                  panel.grid.minor = ggplot2::element_blank())+
            ggplot2::xlab(xlab2[k])+
            ggplot2::ylab(ylab[i])
          plotlist[[k]] <- gg.temp
        }
      }
    }
    
    if(nvar>1){
      if(is.null(which)){
        warning("There are more than one variate in this mcnode. Output is a plotlist.
                You can get a single variate plot using list[[n]] (see example)",
                call.=FALSE)
        return(invisible(plotlist))
      } else {return(plotlist[[which]])}
    } else {return(plotlist[[1]])}
    
  })
}




