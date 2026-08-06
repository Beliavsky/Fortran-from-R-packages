\name{SummaryPca-class}
\docType{class}
\alias{SummaryPca-class}
\alias{show,SummaryPca-method}

\title{Class "SummaryPca" - summary of "Pca" objects  }

\description{The "Pca" object plus some additional summary information}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummaryPca", ...)},
but most often by invoking 'summary' on a "Pca" object. They contain values 
meant for printing by 'show'.
}
\section{Slots}{
    \describe{
        \item{\code{pcaobj}:}{Object of class \code{"Pca"}}
        \item{\code{importance}:}{matrix with additional information: importance of components }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SummaryPca")}: display the object }
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
  \code{\link{Pca-class}}
}
\examples{
showClass("SummaryPca")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
