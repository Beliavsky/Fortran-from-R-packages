\name{SummaryQda-class}
\docType{class}
\alias{SummaryQda-class}
\alias{show,SummaryQda-method}

\title{ Class "SummaryQda" - summary of "Qda" objects   }
\description{ 
Summary information about a \code{Qda} - Quadratic Discriminant Analysis object }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummaryQda", ...)},
but most often by invoking 'summary' on an "Qda" object. They contain values 
meant for printing by 'show'.
}
\section{Slots}{
    \describe{
        \item{\code{qdaobj}:}{Object of class \code{"Qda"}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SummaryQda")}: display the object }
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
showClass("SummaryQda")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
