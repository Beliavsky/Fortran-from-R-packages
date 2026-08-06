\name{maryo}
\alias{maryo}
\docType{data}
\title{ Marona and Yohai Artificial Data}
\description{
 Simple artificial data set generated according the example by 
Marona and Yohai (1998). 
The data set consists of 20 bivariate normal observations generated
with zero means, unit variances and correlation 0.8. The sample 
correlation is 0.81. Two outliers are introduced (i.e. these are 
10\% of the data) in the following way: two points are modified 
by interchanging the largest (observation 19) and smallest 
(observation 9) value of the first coordinate. The sample correlation 
becomes 0.05. 
This example  provides a good example of the fact that a 
multivariate outlier need not be an outlier in any of its 
coordinate variables.
}
\usage{data(maryo)}
\format{
  A data frame with 20 observations on 2 variables. To introduce the outliers
  x[9,1] with x[19,1] are interchanged.

%  \describe{
%    \item{V1}{x[1]}
%    \item{V2}{x[2]}
%  }
}
%\details{}
\source{
    R. A. Marona and V. J. Yohai (1998) Robust estimation of multivariate 
    location and scatter. In
    \emph{Encyclopedia of Statistical Sciences, Updated Volume 2} 
    (Eds. S.Kotz, C.Read and D.Banks). Wiley, New York p. 590
}
%\references{}
\examples{

data(maryo)
getCorr(CovClassic(maryo))          ## the sample correlation is 0.81

## Modify 10%% of the data in the following way:
##  modify two points (out of 20) by interchanging the 
##  largest and smallest value of the first coordinate
imin <- which(maryo[,1]==min(maryo[,1]))        # imin = 9
imax <- which(maryo[,1]==max(maryo[,1]))        # imax = 19
maryo1 <- maryo
maryo1[imin,1] <- maryo[imax,1]
maryo1[imax,1] <- maryo[imin,1]

##  The sample correlation becomes 0.05
plot(maryo1)
getCorr(CovClassic(maryo1))         ## the sample correlation becomes 0.05
getCorr(CovMcd(maryo1))      ## the (reweighted) MCD correlation is 0.79

}
\keyword{datasets}
