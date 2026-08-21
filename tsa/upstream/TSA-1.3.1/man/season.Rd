\name{season}
\alias{season}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract the season info from a time series}
\description{
Extract the season info from a equally spaced 
time series and create a vector of the
season info. For example for monthly data, the function outputs a
vector containing the months of the data. 
}
\usage{
season(x, labels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a time series}
  \item{labels}{ the user supplied labels for the seasons}
}
\details{
The time series must have frequency greater than 1, otherwise the 
function will stop and issue an error message.
If labels is missing, labels will be set as follows:
It is set to be c("1Q","2Q","3Q","4Q)
if the frequency of x equals 4, c("January",...,"December")
if the frequency equals 12, and c("Monday",...,"Sunday") if
frequency equals 7. Otherwise, it is set to be c("S1",...)
}
\value{
An invisible vector containing the seasons of the data 
}
\author{Kung-Sik Chan}
\seealso{ \code{\link{harmonic}}}
\examples{
data(tempdub)
month.=season(tempdub) # the period sign is included to make the printout from
# the commands two line below clearer; ditto below.
model2=lm(tempdub~month.-1) # -1 removes the intercept term 
summary(model2)
}
\keyword{methods}
