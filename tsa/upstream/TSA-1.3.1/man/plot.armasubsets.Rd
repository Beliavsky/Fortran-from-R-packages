\name{plot.armasubsets}
\alias{plot.armasubsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the Best Subset ARMA models}
\description{
This function is adapted from the plot.regsubsets function of the leaps
package, and its main use is to plot the output from the
armasubsets function.}
\usage{
\method{plot}{armasubsets}(x, labels = obj$xnames, main = NULL, 
scale = c("BIC", "AICc", "AIC", "Cp", "adjR2", "R2"), 
col = gray(c(seq(0.4, 0.7, length = 10), 0.9)), draw.grid = TRUE, 
axis.at.3 = TRUE, ...)
}

\arguments{
  \item{x}{an object of class armasubsets}
  \item{labels}{variable names}
  \item{main}{title for plot}
  \item{scale}{which summary statistic to use for ordering plots}
  \item{col}{the last color should be close to but distinct from white}
  \item{draw.grid}{a logical argument; if it is true (default), gray grid lines are superimposed on the graph.}
  \item{axis.at.3}{a logical argument; if if it is true (default), the x-labels are drawn on the upper
 horizontal axis.}
  \item{\dots}{other arguments}
}

\value{Plot the few best subset ARMA models.}
\author{Kung-Sik Chan, based on previoud work by Thomas Lumley and Merlise Clyde}
\seealso{armasubsets}
\examples{
set.seed(53331)
test=arima.sim(model=list(ar=c(rep(0,11),.8),ma=c(rep(0,11),0.7)),n=120)
res=armasubsets(y=test,nar=14,nma=14,y.name='test',ar.method='ols')
plot(res)
}
\keyword{methods}
