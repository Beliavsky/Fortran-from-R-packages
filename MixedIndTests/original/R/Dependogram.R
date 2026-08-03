#'@title Dependogram for Cramer-von Mises statistics
#'
#'@description This function, used in EstDep, TestIndCopula and TestIndSerCopula, draws the P-values of the Moebius Cramer-von Mises statistics from the multilinear copula and their combination for a tests of randomness for k consectives values X(1), ..., X(k) or for a test of independence between random variables.
#'
#'@param out      List of the output from EstDep, EstDepSerial, TestIndCopula or TestIndSerCopula (P-values, subsets)
#'@param stat     Name of statistics to be used (default is "CVM")
#
#'
#'@references Genest, Neslehova, Remillard & Murphy (2019). Testing for independence in arbitrary distributions
#'@examples
#' x <- matrix(rnorm(250),ncol=5)
#' out <-TestIndCopula(x)
#' Dependogram(out)
#'
Dependogram = function(out,stat="CVM")
{

  if(stat=="CVM"){
    pval = out$pvalue$cvm
  }else{
    pval = out$pvalues
  }
  m = length(pval)
  subsets = out$subsets



  mycol = c("black","red")
  x=c(1:m)
  Sig = as.factor(as.numeric(pval<5))
  x.framed <- data.frame(x, pval,factor=Sig)
  new_theme <- theme_update(
    axis.text.x  = element_text(angle=90, vjust=0.5, size=8))

  plot <- ggplot(x.framed, aes(x, pval))
  plot <- plot + new_theme
  plot <- plot + ggtitle(paste("Dependogram of ", stat, "tests")) + ylab(paste("P-value (%) of", stat))+xlab("Subsets")
  plot <- plot + geom_point(stat = "identity")
  plot <- plot+geom_point(aes(color=Sig))+ scale_color_manual(values = mycol)

  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,pval[i])
    lines <- data.frame(xs,ys)
    plot <- plot + geom_path(data = lines,  aes(xs,ys))
  }

  plot <- plot + geom_hline(yintercept=5, color = "red",lty=3)
  plot <- plot + geom_hline(yintercept=0, color = "blue")
  plot <- plot + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  scale_fill_manual(values = mycol)
  print(plot)


}


