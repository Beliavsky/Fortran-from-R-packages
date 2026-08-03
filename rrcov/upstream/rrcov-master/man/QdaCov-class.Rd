\name{QdaCov-class}
\docType{class}
\alias{QdaCov-class}
\title{Class "QdaCov" - Robust methods for Quadratic Discriminant Analysis }
\description{ 
    Robust quadratic discriminant analysis is performed by replacing the classical 
    group means and withing group covariance matrices by their robust equivalents. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("QdaCov", ...)} but the 
    usual way of creating \code{QdaCov} objects is a call to the function
    \code{QdaCov} which serves as a constructor.
}
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
Class \code{"\linkS4class{QdaRobust}"}, directly.
Class \code{"\linkS4class{Qda}"}, by class "QdaRobust", distance 2.
}

\section{Methods}{No methods defined with class "QdaCov" in the signature.}

\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{QdaRobust-class}}, \code{\link{Qda-class}}, \code{\link{QdaClassic}}, \code{\link{QdaClassic-class}}
}

\examples{
showClass("QdaCov")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
