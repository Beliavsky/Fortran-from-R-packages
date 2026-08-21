\name{tlrt}
\alias{tlrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood ratio test for threshold nonlinearity}
\description{
Carry out the likelihood ratio test for threshold nonlinearity, with
the null hypothesis being a normal AR process and the alternative
hypothesis being a TAR model with homogeneous, normally distributed errors. 
}
\usage{
tlrt(y, p, d = 1, transform = "no", a = 0.25, b = 0.75,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{time series}
  \item{p}{working AR order}
  \item{d}{delay}
  \item{transform}{available transformations: "no" (i.e. use raw data), "log", "log10" and
"sqrt"}
  \item{a}{lower percent; the threshold is searched over the interval defined by the
a*100 percentile to the b*100 percentile of the time-series variable} 
  \item{b}{upper percent}
  \item{...}{other arguments to be passed to the ar function which determines the 
Ar order, if p is missing}
}
\details{
The search for the threshold parameter may be narrower than that defined by the
user as the function attempts to ensure adequate sample size in each 
regime of the TAR model.
The p-value of the test is based on large-sample approximation and also 
is more reliable for small p-values. 
}
\value{
  \item{p.value}{p-value of the test}
  \item{test.statistic}{likelihood ratio test statistic}
  \item{a}{the actual lower fraction that defines the interval of search
for the threshold; it may differ from the a specified by the user} 
  \item{b}{the actual upper fraction that defines the interval of search
for the threshold} 
}
\references{
Chan, K.S. (1990).  Percentage points of likelihood ratio tests for
   threshold autoregression.   Journal of Royal Statistical  Society, B  53, 3,  691-696.}
\author{Kung-Sik Chan}
\seealso{\code{\link{Keenan.test}}, \code{\link{Tsay.test}}}
\examples{
data(spots)
pvaluem=NULL
for (d in 1:5){
res=tlrt(sqrt(spots),p=5,d=d,a=0.25,b=0.75)
pvaluem= cbind( pvaluem, round(c(d,signif(c(res$test.statistic, 
res$p.value))),3))
}
rownames(pvaluem)=c('d','test statistic','p-value')
pvaluem
}
\keyword{methods}
