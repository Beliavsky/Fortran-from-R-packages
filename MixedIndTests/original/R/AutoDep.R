#'@title Dependogram for Kendall's tau and Spearman's rho
#'
#'@description This function, used in EstDepSerial, draws the P-values of Kendall's tau and Spearman's rho for a given number of lags.
#'
#'@param out      List of the output of EstDepSerial (P-values, subsets)
#'
#
#'@references B.R Nasri (2021). Tests of serial dependence for arbitrary distributions
#'
#'@examples
#' out <-EstDepSerial(SimAR1Poisson(c(5,0.4),100),10)
#' AutoDep(out)
#'
AutoDep = function(out)
{
  pval = out$pvaltau

  m = length(pval)
  subsets = out$subsets

  lags=c(1:m)


  mycol = c("black","red")
  Sig = as.factor(as.numeric(pval<5))

  x.framed <- data.frame(lags, pval,factor=Sig)
  new_theme <- theme_update(
    axis.text.x  = element_text(angle=90, vjust=0.5, size=8))

  plot <- ggplot(x.framed, aes(lags, pval))
  plot <- plot + new_theme
  plot <- plot + ggtitle("Dependogram of Kendall's taus") + ylab("P-value (%) of tau")+xlab("Lag")
  plot <- plot + geom_point(stat = "identity")
  plot <- plot+geom_point(aes(color=Sig))+ scale_color_manual(values = mycol)

  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,pval[i])
    lines <- data.frame(xs,ys)
    plot <- plot + geom_path(data = lines,  aes(xs, ys))
  }

  plot <- plot + geom_hline(yintercept=5, color = "red",lty=3)
  plot <- plot + geom_hline(yintercept=0, color = "blue")
  plot <- plot + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  scale_fill_manual(values = mycol)
  print(plot)
  #

  pval = out$pvaltau
  Sig = as.factor(as.numeric(pval<5))

  x.framed <- data.frame(lags, pval,factor=Sig)
  new_theme <- theme_update(
    axis.text.x  = element_text(angle=90, vjust=0.5, size=8))

  plot <- ggplot(x.framed, aes(lags, pval))
  plot <- plot + new_theme
  plot <- plot + ggtitle("Dependogram of Spearman's rhos") + ylab("P-value (%) of rho")+xlab("Lag")
  plot <- plot + geom_point(stat = "identity")
  plot <- plot+geom_point(aes(color=Sig))+ scale_color_manual(values = mycol)

  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,pval[i])
    lines <- data.frame(xs,ys)
    plot <- plot + geom_path(data = lines,  aes(xs, ys))
  }

  plot <- plot + geom_hline(yintercept=5, color = "red",lty=3)
  plot <- plot + geom_hline(yintercept=0, color = "blue")
  plot <- plot + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  scale_fill_manual(values = mycol)
  print(plot)
}
