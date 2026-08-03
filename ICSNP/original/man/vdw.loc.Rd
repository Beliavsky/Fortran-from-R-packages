\name{vdw.loc}
\alias{vdw.loc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Van der Waerden Estimator of Location}
\description{
  Iterative algorithm to compute the location estimator based on van der Waerden scores (sometimes also referred to as normal scores).
  }
\usage{
vdw.loc(x, int.diff = 10, maxiter = 1000, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector.}
  \item{int.diff}{number of observations in internal interval when the estimate is searched. }
  \item{maxiter}{maximum number of iterations.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  The algorithm searches among the observations and all Walsh averages for the  two points nearest around the root of the van der Waerden score criterion. Since the criterion function
  is monotone first the \code{int.diff} of the sorted data points are searched that contain the root. After then determining there the two points of question a linear interpolation is used as an estimate.    
}
\value{
  the van der Waerden score estimator of location.
}
\references{ \cite{Hettmansperger, T.P. and McKean, J.W. (1998), \emph{Robust Nonparametric Statistical Methods},
London, Arnold. }}
\author{Klaus Nordhausen}

\examples{
set.seed(1)
x <- rt(100, df = 3)
vdw.loc(x)
rm(.Random.seed)
}
\keyword{nonparametric}
\keyword{robust}
\keyword{univar}

