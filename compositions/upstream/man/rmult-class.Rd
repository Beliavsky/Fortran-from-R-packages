\name{rmult-class}
\Rdversion{1.1}
\docType{class}
\alias{rmult-class}
\alias{coerce,rmult,data.frame-method}
\alias{coerce,rmult,structure-method}
\alias{coerce<-,rmult,data.frame-method}

\title{Class \code{"rmult"}}
\description{
The S4-version of the data container "rmult" for compositional data. More information in 
\code{\link{rmult}}
}
\section{Objects from the Class}{A virtual Class: No objects may be directly created from it. 
This is provided to ensure that rmult objects behave as data.frame or structure under certain circumstances. Use \code{\link{rmult}} to create these objects.}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} containing the data itself }
    \item{\code{names}:}{Object of class \code{"character"} with column names }
    \item{\code{row.names}:}{Object of class \code{"data.frameRowLabels"} with row names }
    \item{\code{.S3Class}:}{Object of class \code{"character"} with the class string }
  }
}
\section{Extends}{
Class \code{"\linkS4class{data.frame}"}, directly.
Class \code{"\linkS4class{compositional}"}, directly.
Class \code{"\linkS4class{list}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{oldClass}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{vector}"}, by class "data.frame", distance 3.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "rmult", to = "data.frame")}: to generate a data.frame }
    \item{coerce}{\code{signature(from = "rmult", to = "structure")}: to generate a structure (i.e. a vector, matrix or array) }
    \item{coerce<-}{\code{signature(from = "rmult", to = "data.frame")}: to overwrite a composition with a data.frame}
	 }
}
\references{
see \code{\link{rmult}}
}
\author{
Raimon Tolosana-Delgado
}
\note{
see \code{\link{rmult}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
see \code{\link{rmult}}

}
\examples{
showClass("rmult")
}
\keyword{classes}
