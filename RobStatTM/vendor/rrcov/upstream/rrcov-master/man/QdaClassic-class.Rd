\name{QdaClassic-class}
\Rdversion{1.1}
\docType{class}
\alias{QdaClassic-class}

\title{Class "QdaClassic" - Quadratic Discriminant Analysis }
\description{ Contains the results of classical Quadratic Discriminant Analysis }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("QdaClassic", ...)} but the 
    usual way of creating \code{QdaClassic} objects is a call to the function
    \code{QdaClassic} which serves as a constructor.
}\section{Slots}{
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
    \item{\code{control}:}{Object of class \code{"CovControl"} inherited from class \code{Qda} 
    specifying which estimate and with what estimation options to use for the group 
    means and covariances. It is always \code{NULL} for classical discriminant analysis.}    
  
  }
}
\section{Extends}{
Class \code{"\linkS4class{Qda}"}, directly.
}
\section{Methods}{
No methods defined with class "QdaClassic" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{QdaRobust-class}}, \code{\link{Qda-class}}, \code{\link{QdaClassic}} 
}

\examples{
showClass("QdaClassic")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
