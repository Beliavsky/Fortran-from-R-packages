\name{QdaRobust-class}
\docType{class}
\alias{QdaRobust-class}

\title{Class "QdaRobust" is a virtual base class for all robust QDA classes }
\description{ The class \code{QdaRobust} searves as a base class for deriving all other 
classes representing the results of robust Quadratic Discriminant Analysis methods}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
    \item{\code{call}:}{The (matched) function call.}
    \item{\code{prior}:}{Prior probabilities used, default to group proportions}
    \item{\code{counts}:}{number of observations in each class}
    \item{\code{center}:}{the group means}
    \item{\code{cov}:}{the group covariance matrices}
    \item{\code{covinv}:}{the inverse of the group covariance matrices}
    \item{\code{covdet}:}{the determinants of the group covariance matrices}
    \item{\code{method}:}{a character string giving the estimation method used}
    \item{\code{X}:}{the training data set (same as the input parameter x of the constructor function)}
    \item{\code{grp}:}{grouping variable:  a factor specifying the class for each observation.} 
    \item{\code{control}:}{Object of class \code{"CovControl"} specifying which estimate to 
    use for the group means and covariances }    
  }
}
\section{Extends}{
Class \code{"\linkS4class{Qda}"}, directly.
}
\section{Methods}{
No methods defined with class "QdaRobust" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{Qda-class}}, \code{\link{QdaClassic-class}}, 
}

\examples{
showClass("QdaRobust")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
