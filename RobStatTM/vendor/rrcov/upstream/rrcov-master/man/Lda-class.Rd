\name{Lda-class}
\docType{class}
\alias{Lda-class}
\alias{predict,Lda-method}
\alias{show,Lda-method}
\alias{summary,Lda-method}

\title{Class "Lda" - virtual base class for all classic and robust LDA classes 
}
\description{ The class \code{Lda} serves as a base class for deriving 
    all other classes representing the results of classical 
    and robust Linear Discriminant Analisys methods }
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
  
    \item{\code{call}:}{the (matched) function call.}
    \item{\code{prior}:}{prior probabilities used, default to group proportions}
    \item{\code{counts}:}{number of observations in each class}
    \item{\code{center}:}{the group means}
    \item{\code{cov}:}{the common covariance matrix}
    \item{\code{ldf}:}{a matrix containing the linear discriminant functions}
    \item{\code{ldfconst}:}{a vector containing the constants of each linear discriminant function}
    \item{\code{method}:}{a character string giving the estimation method used}
    \item{\code{X}:}{the training data set (same as the input parameter x of the constructor function)}
    \item{\code{grp}:}{grouping variable:  a factor specifying the class for each observation.}   
    \item{\code{covobj}:}{object of class \code{"Cov"} containing the estimate 
        of the common covariance matrix of the centered data. It is not NULL 
        only in case of method "B".} 
    \item{\code{control}:}{object of class \code{"CovControl"} specifying which estimate 
    and with what estimation options to use for the group means and common covariance 
    (or \code{NULL} for classical linear discriminant analysis)}   
  }
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "Lda")}: calculates prediction using the results in 
        \code{object}. An optional data frame or matrix in which to look for variables with which 
        to predict. If omitted, the training data set is used. If the original fit used a formula or 
        a data frame or a matrix with column names, newdata must contain columns with the 
        same names. Otherwise it must contain the same number of columns, 
        to be used in the same order. }
    \item{show}{\code{signature(object = "Lda")}: prints the results }
    \item{summary}{\code{signature(object = "Lda")}: prints summary information }
    
   }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{LdaClassic}}, \code{\link{LdaClassic-class}}, \code{\link{LdaRobust-class}}
}
\examples{
showClass("Lda")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
