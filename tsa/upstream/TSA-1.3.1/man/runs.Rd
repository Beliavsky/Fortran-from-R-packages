\name{runs}
\alias{runs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Runs test}
\description{
Test the independence of a sequence of random variables by checking
whether there are too many or too few 
runs above (or below) the median.  
}
\usage{
runs(x,k=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time series}
  \item{k}{the value above or below which runs are counted; default is zero, so
data is assumed to have zero median}
}
\details{
The runs
test examines the data 
in sequence to look for patterns
that would give
evidence against independence. Runs above or below k 
are counted. A small number
of runs would indicate that neighboring values are 
positively dependent and tend to hang together over time. 
On the other hand, too many runs would indicate that the
data 
oscillate back and forth across their median of zero. 
Then neighboring residuals
are negatively dependent. 
So either too few or too many runs lead us to reject independence.
When applied to residuals, the runs test is useful for model diagnostics.
}
\value{
  \item{pvalue}{p-value of the test}
  \item{observed.runs}{observed number of runs}
 \item{expected.runs}{expected number of runs}
 \item{n1}{number of data less than or equal to k}
 \item{n2}{number of data above k}
}
\author{Kung-Sik Chan}
\examples{
data(tempdub)
month.=season(tempdub) # the period sign is included to make the printout from
# the following command clearer.
model3=lm(tempdub~month.) # intercept is automatically included so one month (Jan) is dropped
summary(model3)
runs(rstudent(model3))
}
\keyword{methods}
