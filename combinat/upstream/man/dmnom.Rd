\name{dmnom}
\alias{dmnom}
\alias{fact}
\alias{logfact}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{  
density of multinomial, and support
functions
        }
\description{
density of multinomial
}
\usage{
dmnom(x, size=sum(x), prob=stop("no prob arg"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ vector }
  \item{size}{ total }
  \item{prob}{ parameter vector (sums to 1) }
}
%\details{
%}
%\value{
%}
\author{ code by Scott Chasalow, R pack and maint by VJ Carey <stvjc@channing.harvard.edu>}
\examples{
 dmnom(c(1,1,4,4),10,c(.2,.2,.3,.3))
}
\keyword{ models }% at least one, from doc/KEYWORDS
%\keyword{ combinatorics }% __ONLY ONE__ keyword per line


