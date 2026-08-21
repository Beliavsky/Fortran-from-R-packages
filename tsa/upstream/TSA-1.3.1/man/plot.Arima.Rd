\name{plot.Arima}
\alias{plot.Arima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute and Plot the Forecasts Based on a Fitted Time Series Model}
\description{ Plots the time series data and its predictions with 95\% prediction bounds.}
\usage{
\method{plot}{Arima}(x, n.ahead = 12, col = "black", ylab = object$series, 
lty = 2, n1, newxreg, transform, Plot=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a fitted arima model}
  \item{n.ahead}{number of prediction steps ahead (default=12)}
  \item{col}{color of the prediction bounds}
  \item{ylab}{label of the y-axis}
  \item{lty}{line type of the point predictor; default=dashed lines}
  \item{n1}{starting time point of the plot (default=earliest time point)}
  \item{newxreg}{a matrix of covariate(s) over the period of prediction}
  \item{transform}{function used to transform the forecasts and their prediction
bounds; if missing, no transformation will be carried out. This option is useful if the model was fitted to the transformed data and it is desirable to obtain the forecasts on the original scale. For example, if the model was fitted with the logarithm of the data, then transform = exp will plot the forecasts and their
prediction bounds on the original scale.} 
  \item{Plot}{Plotting will be suppressed if Plot is set to be FALSE; default is TRUE}
  \item{\dots}{additional parameters passed to the plot function}
}
\value{
Side effects of the function: plot the forecasts and their 95\% prediction
bounds, unless Plot is set to be FALSE.
The part of the observed series is plotted with all data plotted as
open circles and linked by a smooth line. By default the predicted values
are plotted as open circles joined up by a dashed line. The plotting style of
the predicted values can be altered by supplying relevant plotting options, e.g
specifying the options type='o', pch=19 and lty=1 will plot the predicted
values as solid circles that are overlaid on the connecting smooth solid line.
The prediction limits are plotted as dotted lines, with default color being
black. However, the prediction limits can be drawn in 
other colors.
For example, setting col='red' paints the prediction limits in red.  
An interesting use of the col argument is setting col=NULL which has the effect
of not drawing the prediction limits.

The function returns 
an invisible list containing the following components.
\item{pred}{the time series of predicted values}
\item{lpi}{the corresponding lower 95\% prediction limits}
\item{upi}{the corresponding upper 95\% prediction limits}
} 
\author{Kung-Sik Chan}
\examples{
data(oil.price)
oil.IMA11alt=arima(log(oil.price),order=c(0,1,1),
# create the design matrix of the covariate for prediction
xreg=data.frame (constant=seq(oil.price)))
n=length(oil.price)
n.ahead=24
newxreg=data.frame(constant=(n+1):(n+n.ahead))
# do the prediction and plot the results
plot(oil.IMA11alt,n.ahead=n.ahead,newxreg=newxreg,
ylab='Log(Oil Price)',xlab='Year',n1=c(2000,1))
# do the same thing but on the orginal scale
plot(oil.IMA11alt,n.ahead=n.ahead,newxreg=newxreg,
ylab='Oil Price',xlab='Year',n1=c(2000,1),transform=exp,pch=19, lty=1,type='o')
# Setting pch=19 plots the predicted values as solid circles.
res=plot(oil.IMA11alt,n.ahead=n.ahead,newxreg=newxreg,
ylab='Oil Price',xlab='Year',n1=c(2000,1),transform=exp,pch=19,col=NULL)
# Setting col=NULL will make the prediction bands invisible. Try col='red'. 
res
# prints the predicted values and their 95\% prediction limits.

}
\keyword{methods}
