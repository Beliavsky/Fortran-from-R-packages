\name{LdaRobust-class}
\docType{class}
\alias{LdaRobust-class}

\title{Class "LdaRobust" is a virtual base class for all robust LDA classes }
\description{ The class \code{LdaRobust} searves as a base class for deriving all other 
classes representing the results of the robust Linear Discriminant Analysis methods}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
    \item{\code{call}:}{The (matched) function call.}
    \item{\code{prior}:}{Prior probabilities used, default to group proportions}
    \item{\code{counts}:}{number of observations in each class}
    \item{\code{center}:}{the group means}
    \item{\code{cov}:}{the common covariance matrix}
    \item{\code{ldf}:}{a matrix containing the linear discriminant functions}
    \item{\code{ldfconst}:}{a vector containing the constants of each linear discriminant function}
    \item{\code{method}:}{a character string giving the estimation method used}
    \item{\code{X}:}{the training data set (same as the input parameter x of the constructor function)}
    \item{\code{grp}:}{grouping variable:  a factor specifying the class for each observation.}   
  }
}
\section{Extends}{
Class \code{"\linkS4class{Lda}"}, directly.
}
\section{Methods}{
No methods defined with class "LdaRobust" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{Lda-class}}, \code{\link{LdaClassic-class}}, 
}

\examples{
showClass("LdaRobust")
}
\keyword{robust}
\keyword{multivariate}
