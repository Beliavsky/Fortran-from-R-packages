#' @title Draws a Tornado chart as provided by tornado (ggplot version).
#' @name ggtornado
#' @aliases ggtornado
#' @aliases ggtornadounc
#' @description Draws a Tornado chart as provided by tornado. 
#' @usage 
#' ## For class 'tornado'
#' ggtornado(x, 
#'   which=1, 
#'   name=NULL, 
#'   stat=c("median","mean"), 
#'   xlab="method", 
#'   ylab=""
#' )
#' @usage 
#' ## For class 'tornadounc'
#' ggtornadounc(x,
#'   which=1, 
#'   stat="median", 
#'   name=NULL, 
#'   xlab="method", 
#'   ylab=""
#' )
#' @param x A tornado object as provided by the \code{\link{tornado}} function.
#' @param which Which output to print -for multivariates output-.
#' @param name Vector of name of input variables. If NULL, the name will be given from the name of the elements.
#' @param stat The name of the statistics of the output to be considered. For a tornado object: "median" or "mean". For a tornadounc object: the value should match one row name of the tornadounc object. Alternatively, for a tornadounc object, the number of the row may be used.
#' @param xlab Label of the x axis. Default is to use the correlation method used in the tornado object.
#' @param ylab Label of the y axis. Default is empty.
#' @seealso \code{\link{tornado}}
#' @examples
#' data(ec)
#' x <- evalmcmod(ec$modEC2, nsv=100, nsu=100, seed=666)
#' tor <- tornado(x, 7)
#' ggtornado(tor)
#' data(total)
#' ggtornado(tornadounc(total, 10, use="complete.obs"), which=1)
#' @export 
ggtornado <- function(x,which=1,name=NULL,stat=c("median","mean"),xlab="method",ylab="")
{
  val <- x$value[[which]]
  if(is.null(val)) stop("Invalid value for which")
  nc <- ncol(val)
  nr <- nrow(val)
  if(!is.null(name)) {colnames(val) <- (rep(name,nc))[1:nc]}
  
  if(xlab=="method") xlab <- c("Spearman's rho statistic","Kendall's tau statistic","Pearson correlation")[pmatch(x$method,c("spearman","kendall","pearson"))]
  
  #create blank plot
  ggp <- ggplot2::ggplot() +
    ggplot2::ylab(xlab)+
    ggplot2::xlab(ylab)+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(linewidth=0.5, colour = "black"),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size=8))
  
  stat <- match.arg(stat)
  stat <- ifelse(stat=="mean" && nr!=1, 2 ,1)
  val <- val[,order(abs(val[stat,])),drop=FALSE]
  
  
  if(nr==1){
    df <- as.data.frame(t(val))
    ggp <- ggp + 
      ggplot2::geom_bar(data = df, ggplot2::aes(x = rownames(df), y = df[,1]),
                        stat = "identity", 
                        position = "identity",
                        width=0.2)+
      ggplot2::coord_flip()+
      ggplot2::scale_y_continuous(limits = c(-1,1),
                                  breaks=seq(-1,1,0.5),
                                  labels = c("-1","-0.5","0","0.5","1"))+
      ggplot2::scale_x_discrete(limits = rownames(df))+
      ggplot2::geom_point(aes(x=1:nc, y=val[1,]),shape="|", size = 8.5)+
      ggplot2::geom_segment(aes(x=0,y=0,xend=nrow(df)+0.2,yend=0))
    
  } else {
    df <- as.data.frame(t(val))
    ggp <- ggp + 
      ggplot2::geom_bar(data = df, 
                        aes(x = rownames(df), y = df[,stat]),
                        stat = "identity", 
                        position = "identity",
                        width=0.2)+
      ggplot2::coord_flip()+
      ggplot2::scale_y_continuous(limits = c(-1,1),
                                  breaks=seq(-1,1,0.5),
                                  labels = c("-1","-0.5","0","0.5","1"))+
      ggplot2::scale_x_discrete(limits = rownames(df))
    
    nrow.df <- nrow(df) #assit in drawing the x=0 line
    my <- val[stat,] #instore median or mean val
    
    if(nr>3){
      xmax <- xmin <- ymax <- ymin <- NULL
      val <- apply(val[3:nr,],2,range)
      df <- as.data.frame(t(val))
      df[,"xminus"] <- c(1:nc-0.1)
      df[,"xplus"] <- c(1:nc+0.1)
      colnames(df) <- c("ymin","ymax","xmin","xmax")
      ggp <- ggp +
        ggplot2::geom_rect(data = df, 
                           ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+
        ggplot2::geom_point(aes(x=1:nc, y =val[1,]), shape="|", size=8.5)+               #left        
        ggplot2::geom_point(aes(x=1:nc, y =val[2,]), shape="|", size=8.5)+               #right
        ggplot2::geom_point(aes(x=1:nc, y=my),shape="|", size=8.5)                       #mid
    }
    ggp <- ggp + ggplot2::geom_segment(ggplot2::aes(x=0,y=0,xend=nrow.df+0.2,yend=0))
  }
  #ggp <- ggp + geom_text(aes(x=1:nc, y =-1.4),label = rownames(df),size=2.5)  # add text
  return(ggp)
}

#' @rdname ggtornado
#' @export
ggtornadounc <- function(x,which=1, stat="median", name=NULL, xlab="method", ylab="")

{
  statposs <- rownames(x$value[[which]])
  
  if(is.character(stat)) stat <- pmatch(stat, rownames(x$value[[which]]))
  if(is.na(stat)) stop("stat should match with: ",paste(statposs,collapse=", ")) 
  
  x$value <- list(x$value[[which]][stat,,drop=FALSE])
  ggtornado(x,which=1, stat="median", name=name, xlab=xlab,ylab=ylab)
}








