\name{tsdiag.Arimax}
\alias{tsdiag.Arimax}
\title{Model Diagnostics for a Fitted ARIMAX Model}
\description{This function is modified from the tsdiag function of the
stats package.  
}
\usage{
\method{tsdiag}{Arimax}(object, gof.lag, tol = 0.1, col = "red", omit.initial = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted ARIMAX model}
  \item{gof.lag}{maximum lag used in ACF and Ljung-Box tests for the residuals}
  \item{tol}{tolerance (default=0.1); see below}
  \item{col}{color of some warning lines in the figures (default=red)}
  \item{omit.initial}{suppress the initial (d+Ds) residuals if true}
  \item{...}{other arguments to be passed to the acf function} 
}
\value{Side effects: 
Plot the time plot of the standardized residuals. Red dashed line at
       +/-qnorm(0.025/no of data) are added to the plot. Data beyond these lines are deemed outliers,
       based on the Bonferronni rule.
       The ACF of the standardized residuals is plotted.
       The p-values of the Ljung-Box test are plotted using a variety of 
       the first K residuals. K starts at the lag on and beyond which the 
       moving-average weights (in the MA(infinity) representation) 
are less than tol.}
\author{Kung-Sik Chan, based on the tsdiag function of the stats pacakage}
\examples{
data(color)
m1.color=arima(color,order=c(1,0,0))
tsdiag(m1.color,gof=15,omit.initial=FALSE)
}
\keyword{methods}

