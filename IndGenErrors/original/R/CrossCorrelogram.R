#'@title Cross-correlogram
#'
#'@description This function, used in crosscor_2series and crosscor_3series plots the graphs of the cross-correlation statistics.
#'
#'@param object   List of the output (statistics, pvalues) from crosscor_2series and crosscor_3series
#'@param comb           Name (string) of series, e.g.,  comb="(x,y)"
#'@param rot      Rotation of labels (default=0)
#'@return \item{Output}{No values are returned; only the graph is printed}
#'
#'@references  Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'@examples
#'  #Romano-Siegel's example #
#'data(romano_ex)
#'outr = crosscor_3series(romano_ex$x,romano_ex$y,romano_ex$z,5,2)
#'CrossCorrelogram(outr$out123,"{x,y,z}",rot=90)
#'
#'@export

CrossCorrelogram = function(object,comb,rot=0)
{



  subsets = object$subsets
  m=length(subsets)
  n = object$n
  Z = object$stat


  z = 1.96/sqrt(n);


  mycol = c("blue","red")
  x=c(1:m)

  Sig = as.factor(as.numeric(abs(Z) > z))

  x1.framed <- data.frame(x, Z,factor=Sig)

  new_theme <- theme_update(axis.text.x  = element_text(angle=rot, vjust=0.5, size=8))
  new_theme <- theme_update(axis.text.x  = element_text(angle=rot, vjust=0.5, size=8))

  #############################################################
  p1 <- ggplot(x1.framed, aes(x, Z))
  p1 <- p1 + new_theme
  p1 <- p1 + ggtitle(paste("Cross-correlogram for",comb)) + ylab("r_{A,n}")+xlab("Lags")
  p1 <- p1 + geom_point(stat = "identity")
  p1 <- p1+geom_point(aes(color=Sig))+ scale_color_manual(values = mycol)



  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,Z[i])
    lines <- data.frame(xs,ys)
    p1 <- p1 + geom_path(data = lines,  aes(xs,ys))
  }

  p1 <- p1 + geom_hline(yintercept= z, color = "red",lty=3)
  p1 <- p1 + geom_hline(yintercept=-z, color = "red",lty=3)
  p1 <- p1 + geom_hline(yintercept=0, color = "black")
  p1 <- p1 + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  p1 <- p1+scale_fill_manual(values = mycol)
  print(p1)


}


