#'@title Dependogram for Moebius correlations
#'
#'@description This function, used in EstDepMoebius and EstDepSerialMoebius plot the graphs of the correlation statistics of Spearman, van der Waerden and Savage as a function of the subsets for tests of randomness or test of independence between random variables. Under the null hypothesis, the statistics should be independent N(0,1).
#'
#'@param out      List of the output from EstDep, EstDepSerial, TestIndCopula or TestIndSerCopula (P-values, subsets)
#'@param n        Number of observations
#'
#'
#'@references Nasri & Remillard (2023). Tests of independence and randomness for arbitrary data using copula-based covariances
#'@examples
#' x <- matrix(rnorm(250),ncol=5)
#' out <-EstDepMoebius(x)
#' DependogramZ(out,50)
#'
DependogramZ = function(out,n)
{


  m = length(out$card)
  subsets = out$subsets

  Z1 = out$stat$spearman;
  Z2 = out$stat$vdw;
  Z3 = out$stat$savage;

  z1 = 1.96/sqrt(n);
  z2 = -z1;
  #library(ggplot2)
  mycol = c("black","red")
  x=c(1:m)

  Sig1 = as.factor(as.numeric(abs(Z1) > z1))
  Sig2 = as.factor(as.numeric(abs(Z2) > z1))
  Sig3 = as.factor(as.numeric(abs(Z3) > z1))
  x1.framed <- data.frame(x, Z1,factor=Sig1)
  x2.framed <- data.frame(x, Z2,factor=Sig2)
  x3.framed <- data.frame(x, Z3,factor=Sig3)
  new_theme <- theme_update(
    axis.text.x  = element_text(angle=90, vjust=0.5, size=8))

  #############################################################
  p1 <- ggplot(x1.framed, aes(x, Z1))
  p1 <- p1 + new_theme
  p1 <- p1 + ggtitle("Dependogram of Spearman's rho") + ylab("r_{A,n}")+xlab("Subsets")
  p1 <- p1+geom_point(aes(color=Sig1))+ scale_color_manual(values = mycol)
  p1 <- p1 + geom_point(stat = "identity")

  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,Z1[i])
    lines <- data.frame(xs,ys)
    p1 <- p1 + geom_path(data = lines,  aes(xs,ys))
  }

  p1 <- p1 + geom_hline(yintercept= z1, color = "red",lty=3)
  p1 <- p1 + geom_hline(yintercept=-z1, color = "red",lty=3)
  p1 <- p1 + geom_hline(yintercept=0, color = "black")
  p1 <- p1 + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  scale_fill_manual(values = mycol)
  print(p1)


  #############################################################
  p2 <- ggplot(x2.framed, aes(x, Z2))
  p2 <- p2 + new_theme
  p2 <- p2 + ggtitle("Dependogram of van der Waerden's coefficent") + ylab("r_{A,n}")+xlab("Subsets")
  p2 <- p2+geom_point(aes(color=Sig2))+ scale_color_manual(values = mycol)
  p2 <- p2 + geom_point(stat = "identity")
  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,Z2[i])
    lines <- data.frame(xs,ys)
    p2 <- p2 + geom_path(data = lines,  aes(xs,ys))
  }

  p2 <- p2 + geom_hline(yintercept= z1, color = "red",lty=3)
  p2 <- p2 + geom_hline(yintercept= z2, color = "red",lty=3)
  p2 <- p2 + geom_hline(yintercept=0, color = "black")
  p2 <- p2 + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  scale_fill_manual(values = mycol)

  print(p2)

  #############################################################

  #############################################################
  p3 <- ggplot(x3.framed, aes(x, Z3))
  p3 <- p3 + new_theme
  p3 <- p3 + ggtitle("Dependogram of Savage's coefficent") + ylab("r_{A,n}")+xlab("Subsets")
  p3 <- p3 + geom_point(stat = "identity")
  p3 <- p3+geom_point(aes(color=Sig3))+ scale_color_manual(values = mycol)
  p3 <- p3 + geom_point(stat = "identity")
  for(i in 1:m){
    xs=c(i,i)
    ys=c(0,Z3[i])
    lines <- data.frame(xs,ys)
    p3 <- p3 + geom_path(data = lines,  aes(xs,ys))
  }

  p3 <- p3 + geom_hline(yintercept= z1, color = "red",lty=3)
  p3 <- p3 + geom_hline(yintercept= z2, color = "red",lty=3)
  p3 <- p3 + geom_hline(yintercept=0, color = "black")
  p3 <- p3 + scale_x_continuous(breaks = 1:m, labels = as.vector(subsets))
  scale_fill_manual(values = mycol)
  print(p3)

  #############################################################

}


