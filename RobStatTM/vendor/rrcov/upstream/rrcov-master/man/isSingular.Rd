\name{isSingular-methods}
\docType{methods}
\alias{isSingular}
\alias{isSingular-methods}
\alias{isSingular,ANY-method}
\alias{isSingular,Cov-method}
\title{
    Check if a covariance matrix (object of class 'Cov') is singular 
}
\description{
  Returns TRUE if the covariance matrix contained in a \code{\link{Cov-class}} object (or derived from) is singular.
}
\usage{
    \S4method{isSingular}{Cov}(obj)
}
\arguments{
  \item{obj}{an object of class (derived from) \code{"Cov"}.}
}

\section{Methods}{
\describe{

    \item{isSingular}{\code{signature(x = Cov)}: Check if a covariance matrix 
    (object of class \code{\link{Cov-class}}) is singular.}
}}

\seealso{
    \code{\link{Cov-class}},
    \code{\link{CovClassic}},
    \code{\link{CovRobust-class}}.
}
\examples{

data(hbk)
cc <- CovClassic(hbk)
isSingular(cc)
}
\keyword{multivariate}

