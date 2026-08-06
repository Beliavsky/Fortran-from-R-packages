\name{SummaryCovRobust-class}
\docType{class}
\alias{SummaryCovRobust-class}
\alias{isClassic,SummaryCovRobust-method}
\alias{show,SummaryCovRobust-method}

\title{Class "SummaryCovRobust" - summary of "CovRobust" objects }
\description{    Summary information for CovRobust objects meants for printing by 'show' }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummaryCovRobust", ...)}, 
but most often by invoking 'summary' on an "Cov" object. They contain values 
meant for printing by 'show'.
}
\section{Slots}{
     \describe{
    \item{\code{covobj}:}{Object of class \code{"Cov"}}
    \item{\code{evals}:}{Eigenvalues of the covariance or correlation matrix}
  }
}
\section{Extends}{
Class \code{"SummaryCov"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SummaryCovRobust")}: ... }
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
  \code{\link{CovRobust-class}},   \code{\link{SummaryCov-class}}
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
cv <- CovMest(hbk.x)
cv
summary(cv)
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
