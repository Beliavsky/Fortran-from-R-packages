#' @title ggplotmc
#' @description Plots the empirical cumulative distribution function of a [mcnode] 
#' or a [mc] object ("`0`" and "`V`" nodes) or the empirical cumulative distribution 
#' function of the estimate of a [mcnode] or [mc] object ("`U`" and "`VU`" nodes) based on 
#' [ggplot2::ggplot] package.
#' @param x  and `mc` or an `mcnode` object  
#' @param prec the precision of the plot. 0.001 will provide an ecdf using the 0.000, 0.001, .002, ..., 1.000 quantiles.
#' @param stat the function used for estimates (2D `mc` or `mcnode`). By default the median.
#' @param lim a  vector of numbers (between 0 and 1) indicating the envelope (2D `mc` or `mcnode`) . Maybe NULL or empty.
#' @param na.rm Should `NA` values be discarded
#' @param griddim a vector of two integers, indicating the size of the grid of the graph. If NULL, the grid is calculated to produce a "nice" graph.
#' @param xlab vector of labels for the x-axis. If `NULL`, the name of the node is used.
#' @param ylab vector of labels for the y-axis.
#' @param main vector of main titles of the graph
#' @param paint Should the envelopes be filled?
#' @param xlim x coordinate range. `xlim` is either a vector of length 2, used for each graph, or a list of vectors of length 2, whose ith element is used for the ith graph. By default, the data range is used as xlim.
#' @param ylim y coordinate range. `ylim` is either a vector of length 2, used for each graph, or a list of vectors of length 2, whose ith element is used for the ith graph. By default, the data range is 0-1.
#' @param which An argument used for an `mcnode` with multivariates. Can specify which variate plot to display. When variates are more than one, the output will be saved in a plot list by default or use the number of which variate to display.
#' @param ... further arguments to be passed to [ggplot2::stat_ecdf()]
#' 
#' @aliases ggplotmc.mc
#' @aliases ggplotmc.mcnode
#' 
#' @returns a ggplot object.
#' 
#' @author Yu Chen and Regis Pouillot
#' 
#' @seealso [plot.mc()]
#' 
#' @examples
#' data(total)
#' # When mcnode has one variate
#' ggplotmc(xV)
#' # Post process
#' ggplotmc(xV) + ggplot2::ggtitle("post processed")
#' # When mcnode has two variates
#' gplots <- ggplotmc(xVUM) #will save two plots in a list
#' gplots[[1]] # show the first variate plot of xVUM mcnode
#' ggplotmc(xVUM, which = 1) #directly show the first variate plot of xVUM mcnode
#' 
#' @rdname ggplotmc
#' @export
ggplotmc <- function(x,...){
  UseMethod("ggplotmc")
}

#' @rdname ggplotmc
#' @method ggplotmc mcnode
#' @export
ggplotmc.mcnode <- function(x, prec=0.001, 
                 stat = c("median","mean"), 
                 lim = c(0.025, 0.25, 0.75, 0.975), 
                 na.rm=TRUE, griddim = NULL, 
                 xlab = NULL, ylab = "Fn(x)", 
                 main = "", paint=TRUE, 
                 xlim=NULL,ylim=NULL,which=NULL,...)
{
  if(inherits(x,"mcnode")==TRUE) {
    x <- quantile.mcnode(x, probs=seq(0,1,prec),lim = lim, na.rm=na.rm, lnames=xlab)
  }
  
    
  y <- x                           # for a correct return
  stat <- match.arg(stat)
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
  
  noms <- names(rapply(x,function(x) 1))    #moche mais efficace
  if(is.null(xlab)) xlab <- noms
  n <- length(noms)
  
  if(!is.null(ylim) & ((is.list(ylim) & length(ylim)!= n)|(is.vector(ylim) & length(ylim)!= 2))) stop("ylim should be NULL, a vector of 2 elements or a list of length the number of nodes") 
  if(!is.null(xlim) & ((is.list(xlim) & length(xlim)!= n)|(is.vector(xlim) & length(xlim)!= 2))) stop("xlim should be NULL, a vector of 2 elements or a list of length the number of nodes") 
  
  main <- rep(main,n)
  xlab <- rep(xlab,n)
  ylab <- rep(ylab,n)
  
  if(is.null(griddim)) griddim <- beau(n)
  
  try({   
    
    i <- 1
    env <- environment()
    plotlist <- list()
    
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
      
      #ggplot ecdf
      ggp <- ggplot2::ggplot()+
        ggplot2::stat_ecdf(ggplot2::aes(unlist(y[1,])),geom = "step",na.rm = FALSE,...)+
        ggplot2::xlim(xlima)+
        ggplot2::ylim(ylima)+
        ggplot2::xlab(xlab[i])+
        ggplot2::ylab(ylab[i])
      
      colorlist <- list()
      
      if(nr > 1){
        rankplot <- 1 + order(-abs(lim-0.5)) 
        for(j in rankplot) {
          colorlist <- c(colorlist,grey(abs(lim[j-1]-.5)+.25))
          df <- data.frame(x = unlist(y[j,]))
          ggp <- ggp +
            ggplot2::stat_ecdf(data=df, mapping=ggplot2::aes(x=x),
                               geom = "step",na.rm = FALSE,...)
        }
        if(paint){
          ggbuild <- ggplot2::ggplot_build(ggp)    #can achieve data
          ti.l <- ggbuild$data[[1]]$x[-length(ggbuild$data[[1]]$x)]  #Points for the polygon used to fill the envelope
          ti.r <- ggbuild$data[[1]]$x[-1L]
          y50 <- ggbuild$data[[1]]$y[-length(ggbuild$data[[1]]$y)]
          thex50 <- rev(as.vector(rbind(ti.l,ti.r)))
          they50 <- rev(as.vector(rbind(y50, y50)))
          for (j in 1:length(lim)){
            ti.lp <- ggbuild$data[[j+1]]$x[-length(ggbuild$data[[j+1]]$x)]
            ti.rp <- ggbuild$data[[j+1]]$x[-1L]
            yp <- ggbuild$data[[j+1]]$y[-length(ggbuild$data[[j+1]]$y)]
            thexp <- as.vector(rbind(ti.lp,ti.rp))
            theyp <- as.vector(rbind(yp, yp))
            df <- data.frame(x=c(thexp,thex50), y=c(theyp,they50))
            ggp<-ggp+
              ggplot2::geom_polygon(data=df, ggplot2::aes(x=x, y=y),
                           fill=colorlist[[j]],
                           col="gray30")
          }
        }
      }
      
      #Add some aesthetics
      ggp <- ggp +
        ggplot2::geom_hline(ggplot2::aes(yintercept=0),linetype="dashed")+
        ggplot2::geom_hline(ggplot2::aes(yintercept=1),linetype="dashed")+
        ggplot2::ggtitle(main[i])+
        ggplot2::theme_bw()+
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
              panel.grid.minor = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5))
      plotlist <<- c(plotlist, list(ggp))
      assign("i",i+1,envir=env)
      return(NULL)}
    
    rapply(y,LEPLOT,how="list")
    
    if(n>1){
      if (is.null(which)){
        warning("There are more than one variate in this mcnode. Output is a plotlist. You can get a specific variate plot by using list[[n]] (see Example)")
        return(invisible(plotlist))
      } else {return(plotlist[[which]])}
    } else {
      return(plotlist[[1]])
    }
    }
  )
  
}

#' @rdname ggplotmc
#' @method ggplotmc mc
#' @export
ggplotmc.mc <- function(x, prec=0.001, 
                     stat = c("median","mean"), 
                     lim = c(0.025, 0.25, 0.75, 0.975), 
                     na.rm=TRUE, griddim = NULL, 
                     xlab = NULL, ylab = "Fn(x)", 
                     main = "", paint=TRUE, 
                     xlim=NULL,ylim=NULL,...)
{
  if(inherits(x,"mc")==TRUE) {
    x <- quantile.mc(x, probs=seq(0,1,prec),lim = lim, na.rm=na.rm, lnames=xlab)
  }
  
  y <- x                           # for a correct return
  stat <- match.arg(stat)
  
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
  
  noms <- names(rapply(x,function(x) 1))    #moche mais efficace
  if(is.null(xlab)) xlab <- noms
  n <- length(noms)
  
  if(!is.null(ylim) & ((is.list(ylim) & length(ylim)!= n)|(is.vector(ylim) & length(ylim)!= 2))) stop("ylim should be NULL, a vector of 2 elements or a list of length the number of nodes") 
  if(!is.null(xlim) & ((is.list(xlim) & length(xlim)!= n)|(is.vector(xlim) & length(xlim)!= 2))) stop("xlim should be NULL, a vector of 2 elements or a list of length the number of nodes") 
  
  main <- rep(main,n)
  xlab <- rep(xlab,n)
  ylab <- rep(ylab,n)
  
  if(is.null(griddim)) griddim <- beau(n)
  
  try({   
    
    i <- 1
    env <- environment()
    
    plotlist <- list()
    
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
      
      #ggplot ecdf
      ggp <- ggplot2::ggplot()+
        ggplot2::stat_ecdf(ggplot2::aes(unlist(y[1,])),geom = "step",na.rm = FALSE,...)+
        ggplot2::xlim(xlima)+
        ggplot2::ylim(ylima)+
        ggplot2::xlab(xlab[i])+
        ggplot2::ylab(ylab[i])
      
      colorlist <- list()
      
      if(nr > 1){
        rankplot <- 1 + order(-abs(lim-0.5)) 
        for(j in rankplot) {
          colorlist <- c(colorlist,grey(abs(lim[j-1]-.5)+.25))
          df <- data.frame(x=unlist(y[j,]))
          ggp <- ggp+
            ggplot2::stat_ecdf(data=df, mapping=ggplot2::aes(x=x),geom = "step",na.rm = FALSE,...)
        }
        if(paint){
          ggbuild <- ggplot2::ggplot_build(ggp)    #can achieve data
          ti.l <- ggbuild$data[[1]]$x[-length(ggbuild$data[[1]]$x)]  #Points for the polygon used to fill the envelope
          ti.r <- ggbuild$data[[1]]$x[-1L]
          y50 <- ggbuild$data[[1]]$y[-length(ggbuild$data[[1]]$y)]
          thex50 <- rev(as.vector(rbind(ti.l,ti.r)))
          they50 <- rev(as.vector(rbind(y50, y50)))
          for (j in 1:length(lim)){
            ti.lp <- ggbuild$data[[j+1]]$x[-length(ggbuild$data[[j+1]]$x)]
            ti.rp <- ggbuild$data[[j+1]]$x[-1L]
            yp <- ggbuild$data[[j+1]]$y[-length(ggbuild$data[[j+1]]$y)]
            thexp <- as.vector(rbind(ti.lp,ti.rp))
            theyp <- as.vector(rbind(yp, yp))
            df <- data.frame(x=c(thexp,thex50), y=c(theyp,they50))
            ggp <- ggp +
              ggplot2::geom_polygon(data=df, ggplot2::aes(x=x, y=y),
                           fill=colorlist[[j]],
                           col="gray30")
          }
        }
      }
      
      #Add some aesthetics
      ggp <- ggp+
        ggplot2::geom_hline(ggplot2::aes(yintercept=0),linetype="dashed")+
        ggplot2::geom_hline(ggplot2::aes(yintercept=1),linetype="dashed")+
        ggplot2::ggtitle(main[i])+
        ggplot2::theme_bw()+
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
              panel.grid.minor = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5))
      
      plotlist <<- c(plotlist, list(ggp))
      assign("i",i+1,envir=env) 
      return(NULL)}
    
    rapply(y,LEPLOT,how="list")

    gga <- ggpubr::ggarrange(plotlist = plotlist,ncol=griddim[2],nrow=griddim[1])
  })

  return(gga)
}






