\name{methodRuleTS1}
\alias{methodRuleTS1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rule for determining method for permTS}
\description{
This is the default function which determines which method to use in \code{\link{permTS}}.  
}
\usage{
methodRuleTS1(x, group, exact, Nbound = c(1000, 200, 100, 50, 16))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of response scores  }
  \item{group}{group membership vector}
  \item{exact}{logical, TRUE=exact method chosen, FALSE=pclt}
  \item{Nbound}{vector of bounds (see details)}
}
\details{
This function determines which of several methods will be used in \code{\link{permTS}}; 
see that help for description of methods. 

When exact=FALSE then returns 'pclt'. When exact=TRUE then returns either 'exact.network' if the 
estimated time of calculation is not too large or 'exact.mc' otherwise. When exact=NULL then returns 
either 'exact.network' if the estimated time is not too large or 'pclt' otherwise. The estimation of the 
calculation time is as follows: if the smallest number of unique values in one of the two groups is equal to 
kmin, then calculation time is large if the sample size <=Nbound[kmin-1], if Nbound[kmin-1] exists, or 
is large if the sample size <= min(Nbound) otherwise.
}
\value{
a character vector with one of the following values: "pclt","exact.network","exact.mc"
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{ \code{\link{permTS}}}
\examples{
N<-100
set.seed(1)
methodRuleTS1(x=sample(1:2,N,replace=TRUE),group=sample(c(0,1),N,replace=TRUE),exact=NULL)
N<-100
methodRuleTS1(sample(1:500,N,replace=TRUE),sample(c(0,1),N,replace=TRUE),TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
