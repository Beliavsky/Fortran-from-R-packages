\name{methodRuleTREND1}
\alias{methodRuleTREND1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rule for determining method for permTREND}
\description{
This is the default function which determines which method to use in \code{\link{permTREND}}.  
}
\usage{
methodRuleTREND1(x, y, exact, Nbound = c(20))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of response scores  }
  \item{y}{group membership vector}
  \item{exact}{logical, TRUE=exact method chosen, FALSE=pclt}
  \item{Nbound}{gives 'pclt' if length(x) > Nbound}
}
\details{
This function determines which of two methods will be used in \code{\link{permTREND}}; 
see that help for description of methods. 

When exact=FALSE then returns 'pclt'. When exact=TRUE then returns  'exact.mc'. When exact=NULL then returns 
either 'exact.mc' if length(x) is less than or equal to Nbound, otherwise 
returns 'pclt'. 
}
\value{
a character vector with one of the following values: "pclt","exact.mc"
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{ \code{\link{permKS}}}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
