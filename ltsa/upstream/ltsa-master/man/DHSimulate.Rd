\name{DHSimulate}
\alias{DHSimulate}

\title{ Simulate General Linear Process }
\description{
Uses the Davies-Harte algorithm to simulate a Gaussian time
series with specified autocovariance function.
}
\usage{
DHSimulate(n, r, ReportTestOnly = FALSE, rand.gen = rnorm, ...)
}
\arguments{
  \item{n}{ length of time series to be generated }
  \item{r}{ autocovariances at lags 0,1,...}
  \item{ReportTestOnly}{ FALSE -- Run normally so terminates with an error if Davies-Harte condition
  does not hold. Othewise if TRUE, then output is TRUE if the Davies-Harte condition holds
  and FALSE if it does not.}
  \item{rand.gen}{ random number generator to use. It is assumed to have mean zero and
  variance one.}
  \item{\dots}{optional arguments passed to \code{rand.gen} }
}
\details{
The method uses the FFT and so is most efficient if the series length, n,
is a power of 2.
The method requires that a complicated non-negativity condition be satisfed.
Craigmile (2003) discusses this condition in more detail and shows for 
anti-persistent time series this condition will always be satisfied.
Sometimes, as in the case of fractinally differenced white noise with
parameter d=0.45 and n=5000, this condition fails and the algorithm doesn't
work.  
In this case, an error message is generated and the function halts.
}
\value{
Either a vector of length containing the simulated time series if Davies-Harte condition
holds and ReportTestOnly = FALSE. 
If argument ReportTestOnly is set to TRUE, then output is logical variable indicating
if Davies-Harte condition holds, TRUE, or if it does not, FALSE.
}

\author{ A.I. McLeod}

\references{ 

Craigmile, P.F. (2003).
Simulating a class of stationary Gaussian processes using the Davies-Harte algorithm, 
with application to long memory processes.
Journal of Time Series Analysis, 24, 505-511.

Davies,  R. B. and Harte, D. S. (1987).
Tests for Hurst Effect. Biometrika 74,  95--101.

McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}

\seealso{ \
\code{\link{DLSimulate}} ,
\code{\link{SimGLP}},
\code{\link{arima.sim}} 
}

\examples{
#simulate a process with autocovariance function 1/(k+1), k=0,1,...
# and plot it
n<-2000
r<-1/sqrt(1:n)
z<-DHSimulate(n, r)
plot.ts(z)

#simulate AR(1) and produce a table comparing the theoretical and sample
# autocovariances and autocorrelations
phi<- -0.8
n<-4096
g0<-1/(1-phi^2)
#theoretical autocovariances
tacvf<-g0*(phi^(0:(n-1)))
z<-DHSimulate(n, tacvf)
#autocorrelations
sacf<-acf(z, plot=FALSE)$acf
#autocovariances
sacvf<-acf(z, plot=FALSE,type="covariance")$acf
tacf<-tacvf/tacvf[1]
tb<-matrix(c(tacvf[1:10],sacvf[1:10],tacf[1:10],sacf[1:10]),ncol=4)
dimnames(tb)<-list(0:9, c("Tacvf","Sacvf","Tacf","Sacf"))
tb

#Show the Davies-Harte condition sometimes hold and sometimes does not
#   in the case of fractionally differenced white noise
#
#Define autocovariance function for fractionally differenced white noise
`tacvfFdwn` <-
function(d, maxlag)
{
    x <- numeric(maxlag + 1)
    x[1] <- gamma(1 - 2 * d)/gamma(1 - d)^2
    for(i in 1:maxlag) 
        x[i + 1] <- ((i - 1 + d)/(i - d)) * x[i]
    x
}
#Build table to show values of d for which condition is TRUE when n=5000
n<-5000
ds<-c(-0.45, -0.25, -0.05, 0.05, 0.25, 0.45)
tb<-logical(length(ds))
names(tb)<-ds
for (kd in 1:length(ds)){
    d<-ds[kd]
    r<-tacvfFdwn(d, n-1)
    tb[kd]<-DHSimulate(n, r, ReportTestOnly = TRUE)
    }
tb

}

\keyword{ ts }
\keyword{ datagen }
