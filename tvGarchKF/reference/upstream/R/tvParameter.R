#' @title Structure of the Time-Varying GARCH(1,1) Parameters
#' @description This function performs an exploratory analysis to uncover the dynamic structure of the time-varying GARCH(1,1) parameters. Specifically, the observation domain \eqn{\{1, \ldots,T\}} is partitioned into \eqn{M} overlapping blocks, each of length \eqn{N}, with a constant shift of size \eqn{S} between consecutive blocks. The relation between these quantities satisfies \eqn{T = S(M-1)+N}. The midpoint of the \eqn{j}-th block, for \eqn{j=1,\ldots,M}, is denoted \eqn{t_j=S(j-1)+N/2}. For each block, a local estimation of the stationary GARCH(1,1) model is performed using the observations within that block. The resulting sequence of local estimates, evaluated across all blocks, provides an empirical trajectory that reflects the underlying evolution of the time-varying parameters. This trajectory can serve as a guide for selecting flexible function classes capable of capturing their temporal variation.
#' @param data Represents the financial return series employed to investigate the temporal evolution of the parameters in the tv-GARCH(1,1) model.
#' @param S The number of observations by which the analysis window is shifted to define the starting point of the next block; also known as the step size or shift parameter.
#' @param N The total number of observations contained within each data block, representing the block or window length over which local model estimation is performed.
#' @param plot A Boolean flag indicating whether a graphical representation of the estimation results should be generated.
#' @return Data frame who contains omega, alpha, beta of GARCH(1,1) model and midpoint each block.
#' @importFrom fGarch garchFit
#' @importFrom graphics abline
#' @importFrom graphics par
#' @export
#' @examples
#' ipsa<-diff(log(indipsa))*100
#' S = 100
#' N = 800
#' tv <- tvParameter(ipsa,S,N)

tvParameter<-function(data, S,N, plot = TRUE){

omega <- vector("numeric")
alpha <- vector("numeric")
beta <- vector("numeric")

y = data
T. = length(y)
S  = S
N  = N
M  = trunc((T.-N)/S+1)
estacionario <- garchFit(~garch(1,1), data=y, include.mean=F, trace = FALSE)
G <- estacionario@fit$coef

for(j in 1:M){
  model= garchFit(~garch(1,1),data=y[(S*(j-1)+1):(S*(j - 1)+N)], include.mean=F, trace = FALSE)
  omega[j]=model@fit$matcoef[1,1]
  alpha[j]=model@fit$matcoef[2,1]
  beta[j]=model@fit$matcoef[3,1]
}

t = S * (1:M - 1) + N/2
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))
if (plot == TRUE) {
  par(mfrow=c(1,3))
  plot(t,omega, xlim=c(0,T.),ylim =c(min(omega)-0.1, max(omega)+0.1), bty = "n", main = expression(~c), las = 1, lwd = 3, ylab = "", cex.axis = 1.2,xlab = "", cex.main = 1.5, pch = 20)
  abline(h=G[1],lwd=2,lty=2)

  plot(t,alpha, xlim=c(0,T.),ylim = c(min(alpha)-0.1, max(alpha)+0.1), bty = "n", main = expression(~alpha), las = 1, lwd = 3, ylab = "", cex.axis = 1.2, xlab = "", cex.main = 1.5, pch = 20)
  abline(h=G[2],lwd=2,lty=2)

  plot(t,beta, xlim=c(0,T.),ylim = c(min(beta)-0.1, max(beta)+0.1), bty = "n", main =expression(~beta), las = 1, lwd = 3, ylab = "", cex.axis = 1.2, xlab = "", cex.main = 1.5, pch = 20)
  abline(h=G[3],lwd=2,lty=2)
}

data <- data.frame(t=t, omega=omega, alpha=alpha, beta=beta)

return(data)
}
