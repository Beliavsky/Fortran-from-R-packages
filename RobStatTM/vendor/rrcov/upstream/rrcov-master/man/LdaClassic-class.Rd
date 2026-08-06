\name{LdaClassic-class}
\docType{class}
\alias{LdaClassic-class}

\title{Class "LdaClassic" - Linear Discriminant Analysis }
\description{ Contains the results of a classical Linear Discriminant Analysis }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LdaClassic", ...)} but the 
    usual way of creating \code{LdaClassic} objects is a call to the function
    \code{LdaClassic} which serves as a constructor.
}\section{Slots}{
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
No methods defined with class "LdaClassic" in the signature.
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{LdaRobust-class}}, \code{\link{Lda-class}}, \code{\link{LdaClassic}} 
}

\examples{
showClass("LdaClassic")
}
\keyword{robust}
\keyword{multivariate}
