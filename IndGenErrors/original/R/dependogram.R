#'@title Dependogram for Cramer-von Mises statistics
#'
#'@description This function, used in cvm_2series and cvm_3series draws the P-values of the Moebius Cramer-von Mises statistics.
#'
#'@param object   List of the output (statistics, pvalues) from cvm_2series and cvmr_3series
#'@param stat     Name (string) of statistics to be used
#'@param rot      Rotation of labels  (default=0)
#
#'@return \item{Output}{No values are returned; only the graph is printed}
#'
#'@references  Duchesne, Ghoudi & Remillard  (2012). On Testing for independence between the innovations of several time series. CJS, vol. 40, 447-479.
#'
#'@examples
#' #Romano-Siegel's example #
#'data(romano_ex)
#'out = cvm_3series(romano_ex$x,romano_ex$y,romano_ex$z,5,2)
#'dependogram(out$out123,"{x,y,z}",rot=90)
#'
#'@export


dependogram = function(object,stat,rot=0)
{

  pval = object$pvalue$cvm

  m = length(pval)
  subsets = object$subsets



  mycol = c("blue","red")
  x=c(1:m)
  Sig = as.factor(as.numeric(pval<0.05))
  x.framed <- data.frame(x, pval,factor=Sig)
  new_theme <- theme_update(axis.text.x  = element_text(angle=rot, vjust=0.5, size=8))
  new_theme <- theme_update(axis.text.x  = element_text(angle=rot, vjust=0.5, size=8))
  plot <- ggplot(x.framed, aes(x, pval))
  plot <- plot + new_theme
  plot <- plot + ggtitle(paste("Dependogram of CVM",stat, "tests")) + ylab("P-value of CVM statistics")+xlab("Lags")
  plot <- plot + geom_point(stat = "identity")
  plot <- plot+geom_point(aes(color=Sig))+ scale_color_manual(values = mycol)

  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,pval[i])
    lines <- data.frame(xs,ys)
    plot <- plot + geom_path(data = lines,  aes(xs,ys))
  }

  plot <- plot + geom_hline(yintercept=0.05, color = "red",lty=3)
  plot <- plot + geom_hline(yintercept=0, color = "black")
  plot <- plot + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  plot <- plot + scale_fill_manual(values = mycol)+ theme(legend.position = "none")
  print(plot)


}


