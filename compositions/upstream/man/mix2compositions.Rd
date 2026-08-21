\name{transformations from 'mixtures' to  'compositions' classes}
\alias{mix.2aplus}
\alias{mix.2acomp}
\alias{mix.2rcomp}
\alias{mix.2rplus}
\alias{mix.2rmult}
\title{Transformations from 'mixtures' to  'compositions' classes}
\description{
  Transformations from 'mixtures' of the "mixR" library to 'compositions'
  classes 'aplus', 'acomp', 'rcomp', 'rplus'
  and 'rmult'.
  }
\usage{
   mix.2aplus(X)
   mix.2acomp(X)
   mix.2rcomp(X)
   mix.2rplus(X)
   mix.2rmult(X)
   }
\arguments{
  \item{X}{mixture object to be converted}
  }
\value{
    A 'compositions' object of the class 'aplus', 'acomp', 'rcomp', 'rplus' or 'rmult'.
 }
\details{
    A 'compositions' object is obtained from the mixtute object m, having the same data matrix
    as mixture object m i.e. \code{m$mat}.
   }
\seealso{
  \code{\link{aplus}}
  \code{\link{acomp}}
  \code{\link{rcomp}}
  \code{\link{rplus}}
  \code{\link{rmult}}  
}
\examples{
\dontrun{
m <- mix.Read("Glac.dat")      # reads the Glacial data set from Aitchison (1986)
m <- mix.Extract(m,c(1,2,3,4))   # mix object with closed four parts subcomposition
ap <- mix.2aplus(m)   # ap is a 'compositions' object of the aplus class
ac <- mix.2acomp(m)   # ac is a 'compositions' object of the acomp class
}

}
\keyword{classes}