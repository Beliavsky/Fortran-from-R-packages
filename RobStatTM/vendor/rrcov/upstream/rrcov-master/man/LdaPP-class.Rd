\name{LdaPP-class}
\docType{class}
\alias{LdaPP-class}
\alias{predict,LdaPP-method}
\title{Class "LdaPP" - Robust method for Linear Discriminant Analysis by Projection-pursuit}
\description{ 
    The class \code{LdaPP} represents an algorithm for robust linear discriminant 
    analysis by projection-pursuit approach. The objects of class \code{LdaPP} contain the results
    of the robust linear discriminant analysis by projection-pursuit approach.
    
    }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LdaPP", ...)} but the 
    usual way of creating \code{LdaPP} objects is a call to the function
    \code{LdaPP} which serves as a constructor.
}
\section{Slots}{
     \describe{
    \item{\code{call}:}{The (matched) function call.}
    \item{\code{prior}:}{Prior probabilities used, default to group proportions}
    \item{\code{counts}:}{number of observations in each class}
    \item{\code{center}:}{the group means}
    \item{\code{cov}:}{the common covariance matrix}
    \item{\code{raw.ldf}:}{a matrix containing the raw linear discriminant functions - see Details in \code{\link{LdaPP}}}
    \item{\code{raw.ldfconst}:}{a vector containing the raw constants of each raw linear discriminant function - see Details in \code{\link{LdaPP}}}
    \item{\code{ldf}:}{a matrix containing the linear discriminant functions}
    \item{\code{ldfconst}:}{a vector containing the constants of each linear discriminant function}
    \item{\code{method}:}{a character string giving the estimation method used}
    \item{\code{X}:}{the training data set (same as the input parameter x of the constructor function)}
    \item{\code{grp}:}{grouping variable:  a factor specifying the class for each observation.}   
  }
}
\section{Extends}{
    Class \code{"\linkS4class{LdaRobust}"}, directly.
    Class \code{"\linkS4class{Lda}"}, by class "LdaRobust", distance 2.
}

\section{Methods}{
\describe{
    \item{predict}{\code{signature(object = "LdaPP")}: calculates prediction using the results in 
        \code{object}. An optional data frame or matrix in which to look for variables with which 
        to predict. If omitted, the training data set is used. If the original fit used a formula or 
        a data frame or a matrix with column names, newdata must contain columns with the 
        same names. Otherwise it must contain the same number of columns, 
        to be used in the same order. If the argument \code{raw=TRUE} is set the raw 
        (obtained by the first approximation algorithm) linear discriminant 
        function and constant will be used.}
    }
}
\references{
    Pires, A. M. and A. Branco, J. (2010)
    Projection-pursuit approach to robust linear discriminant analysis
    \emph{Journal Multivariate Analysis}, Academic Press, Inc., \bold{101}, 2464--2485.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at} and
    Ana Pires \email{apires@math.ist.utl.pt} 
}

\seealso{
     \code{\link{LdaRobust-class}}, \code{\link{Lda-class}}, 
     \code{\link{LdaClassic}}, \code{\link{LdaClassic-class}}, 
     \code{\link{Linda}}, \code{\link{Linda-class}}
}

\examples{
showClass("LdaPP")
}
\keyword{robust}
\keyword{multivariate}
