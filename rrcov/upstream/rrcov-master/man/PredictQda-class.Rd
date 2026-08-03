\name{PredictQda-class}
\docType{class}
\alias{PredictQda-class}
\alias{show,PredictQda-method}

\title{Class "PredictQda" - prediction of "Qda" objects }
\description{ The prediction of a "Qda" object }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PredictQda", ...)}
but most often by invoking 'predict' on a "Qda" object. They contain values 
meant for printing by 'show'}
\section{Slots}{
     \describe{
    \item{\code{classification}:}{a factor variable containing the classification of each object }
    \item{\code{posterior}:}{a matrix containing the posterior probabilities }
    \item{\code{x}:}{ matrix with the discriminant scores }
    \item{\code{ct}:}{re-classification table of the training sample}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "PredictQda")}: prints the results }
     }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
%\note{}
\seealso{
  \code{\link{Qda-class}}
}
\examples{
showClass("PredictQda")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
