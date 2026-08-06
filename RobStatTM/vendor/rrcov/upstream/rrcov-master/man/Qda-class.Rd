\name{Qda-class}
\Rdversion{1.1}
\docType{class}
\alias{Qda-class}
\alias{predict,Qda-method}
\alias{show,Qda-method}
\alias{summary,Qda-method}

\title{Class "Qda" - virtual base class for all classic and robust QDA classes 
}
\description{ The class \code{Qda} serves as a base class for deriving 
    all other classes representing the results of classical 
    and robust Quadratic Discriminant Analisys methods }
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
  
    \item{\code{call}:}{the (matched) function call.}
    \item{\code{prior}:}{prior probabilities used, default to group proportions}
    \item{\code{counts}:}{number of observations in each class}
    \item{\code{center}:}{the group means}
    \item{\code{cov}:}{the group covariance matrices}
    \item{\code{covinv}:}{the inverse of the group covariance matrices}
    \item{\code{covdet}:}{the determinants of the group covariance matrices}
    \item{\code{method}:}{a character string giving the estimation method used}
    \item{\code{X}:}{the training data set (same as the input parameter x of the constructor function)}
    \item{\code{grp}:}{grouping variable:  a factor specifying the class for each observation.}   
    \item{\code{control}:}{object of class \code{"CovControl"} specifying which estimate 
    and with what estimation options to use for the group means and covariances 
    (or \code{NULL} for classical discriminant analysis)}    
  }
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "Qda")}: calculates prediction using the results in 
        \code{object}. An optional data frame or matrix in which to look for variables with which 
        to predict. If omitted, the scores are used. If the original fit used a formula or 
        a data frame or a matrix with column names, newdata must contain columns with the 
        same names. Otherwise it must contain the same number of columns, 
        to be used in the same order. }
    \item{show}{\code{signature(object = "Qda")}: prints the results }
    \item{summary}{\code{signature(object = "Qda")}: prints summary information }
   }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{QdaClassic}}, \code{\link{QdaClassic-class}}, \code{\link{QdaRobust-class}}
}
\examples{
showClass("Qda")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
