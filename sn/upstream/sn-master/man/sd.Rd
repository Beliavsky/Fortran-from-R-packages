\name{sd}
\title{Standard deviation}
\alias{sd}
\alias{sd.default}
 

\description{
The \code{sd} function from the \pkg{stats} is replaced by a new \code{method} 
in order to introduce a separate method to deal with objects of class
\code{SECdistrUv}. The function \code{sd.default} is an alias of the original
function \code{\link[stats]{sd}}.
}

\usage{
  sd(x, \dots)
  \method{sd}{default}(x, na.rm = FALSE, \dots)
}

\arguments{
  \item{x}{a numeric vector, matrix or data frame.}
  \item{na.rm}{logical. Should missing values be removed?}
  \item{\dots}{further arguments passed to or from other methods.}
}

\seealso{
  \code{\link[stats]{sd}}, \code{\linkS4class{SECdistrUv}} 
}


\keyword{univar}
