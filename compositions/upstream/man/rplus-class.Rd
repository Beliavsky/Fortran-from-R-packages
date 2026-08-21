\name{rplus-class}
\Rdversion{1.1}
\docType{class}
\alias{rplus-class}
\alias{coerce,rplus,data.frame-method}
\alias{coerce,rplus,structure-method}
\alias{coerce<-,rplus,data.frame-method}

\title{Class \code{"rplus"}}
\description{
The S4-version of the data container "rplus" for compositional data. More information in 
\code{\link{rplus}}
}
\section{Objects from the Class}{A virtual Class: No objects may be directly created from it. 
This is provided to ensure that rplus objects behave as data.frame or structure under certain circumstances. Use \code{\link{rplus}} to create these objects.}
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
    \item{coerce}{\code{signature(from = "rplus", to = "data.frame")}: to generate a data.frame }
    \item{coerce}{\code{signature(from = "rplus", to = "structure")}: to generate a structure (i.e. a vector, matrix or array) }
    \item{coerce<-}{\code{signature(from = "rplus", to = "data.frame")}: to overwrite a composition with a data.frame}
	 }
}
\references{
see \code{\link{rplus}}
}
\author{
Raimon Tolosana-Delgado
}
\note{
see \code{\link{rplus}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
see \code{\link{rplus}}

}
\examples{
showClass("rplus")
}
\keyword{classes}
