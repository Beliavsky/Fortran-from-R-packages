\name{CovControl-class}
\docType{class}
\alias{CovControl-class}

\title{Class "CovControl" is a VIRTUAL base control class }
\description{  The class "CovControl" is a VIRTUAL base control class
    for the derived classes representing the control parameters
    for the different robust methods }
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\arguments{
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in
    \code{\link{mahalanobis}}.}
}
\section{Methods}{
No methods defined with class "CovControl" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
%\note{}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\seealso{}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
