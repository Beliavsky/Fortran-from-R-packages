\name{Linda-class}
\docType{class}
\alias{Linda-class}
\title{Class "Linda" - Robust method for LINear Discriminant Analysis }
\description{ 
    Robust linear discriminant analysis is performed by replacing the classical 
    group means and withing group covariance matrix by robust equivalents 
    based on MCD. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Linda", ...)} but the 
    usual way of creating \code{Linda} objects is a call to the function
    \code{Linda} which serves as a constructor.
}
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
    \item{\code{l1med}:}{wheather L1 median was used to compute group means.}   
  }
}
\section{Extends}{
Class \code{"\linkS4class{LdaRobust}"}, directly.
Class \code{"\linkS4class{Lda}"}, by class "LdaRobust", distance 2.
}

\section{Methods}{No methods defined with class "Linda" in the signature.}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{LdaRobust-class}}, \code{\link{Lda-class}}, \code{\link{LdaClassic}}, \code{\link{LdaClassic-class}}
}

\examples{
showClass("Linda")
}
\keyword{robust}
\keyword{multivariate}
