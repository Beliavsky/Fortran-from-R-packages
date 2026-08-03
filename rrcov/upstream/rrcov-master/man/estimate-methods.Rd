\name{restimate-methods}
\docType{methods}
\alias{restimate}
\alias{restimate-methods}
\alias{restimate,CovControlMest-method}
\title{ Methods for Function estimate in Package 'rrcov'}
\description{
    Each concrete control class, like \code{CovControlMest}, \code{CovControlOgk},
    etc., should implement an \code{restimate} method which will call the correponding
    (constructor)-function and will return the obtained S4 class, derived from
    \code{CovRobust}.
}
\usage{
\S4method{restimate}{CovControlMest}(obj, x, ...)
}
\arguments{
 \item{obj}{an object of class \code{"CovControlEstimate"}}
  \item{x}{ Data frame or matrix containing the data }.
  \item{... }{ other parameters to be passed through to the estimation function. }
}


\section{Methods}{
\describe{

\item{obj = "CovControlMcd"}{ Compute the MCD estimates of multivariate location and
    scatter by calling\code{CovMcd} }
\item{obj = "CovControlMest"}{ Compute the constrained M-estimates of multivariate location and
    scatter by calling\code{CovMest} }
\item{obj = "CovControlOgk"}{ Compute the Ortogonalized Gnanadesikan-Kettenring (OGK) estimates of multivariate location and
    scatter  by calling\code{CovOgk} }
}}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
