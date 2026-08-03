\name{SummaryLda-class}
\docType{class}
\alias{SummaryLda-class}
\alias{show,SummaryLda-method}

\title{ Class "SummaryLda" - summary of "Lda" objects   }
\description{ 
Contains summary information about an \code{Lda} object - Linear Discriminant Analysis object }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummaryLda", ...)},
but most often by invoking 'summary' on an "Lda" object. They contain values 
meant for printing by 'show'.
}
\section{Slots}{
    \describe{
        \item{\code{ldaobj}:}{Object of class \code{"Lda"}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SummaryLda")}: display the object }
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
  \code{\link{Lda-class}}
}
\examples{
showClass("SummaryLda")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
