\name{fitted.statfacmodBurSt}
\alias{fitted.statfacmodBurSt}
\title{ Variance Matrix From Statistical Factor Model }
\description{
Takes a statistical factor model object and produces the variance matrix
that it represents.
}
\usage{
\method{fitted}{statfacmodBurSt}(object, output = "full", ...)
}
\arguments{
  \item{object}{ 
	an object of class \code{statfacmodBurSt}.
	}
  \item{output}{ 
	a character string that must partially match one of:
	\code{"full"}, \code{"systematic"}, \code{"specific"}.
  }
  \item{...}{ 
	not used, but allows methods for inheriting objects to have
	additional arguments.
  }
}
\value{
a numeric matrix of the variance represented by the factor model, or the
systematic or specific variance portion of the model.

This has two additional attributes:
\item{number.of.factors }{
	the number of factors used in the model.
}
\item{timestamp }{
	the date and time at which the model was originally created.
}
}

\section{Revision }{
This help was last revised 2012 February 12.
}

\author{ Burns Statistics }

\seealso{ 
\code{\link{factor.model.stat}}.
}

\examples{
\dontrun{
varfac <- factor.model.stat(retmat, zero=TRUE, output="factor")
# perhaps modify loadings and uniquenesses for missing values
varian2 <- fitted(varfac) # get variance matrix from factor model

specif <- fitted(varfac, output="specific") # diagonal matrix
}

}
\keyword{multivariate}
