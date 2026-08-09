\name{deviance_o}
\alias{deviance_o}
\encoding{latin1}

\title{Deviance of the orthogonal residuals}

\description{
Returns the deviance \eqn{\sum_{i=1}^N \| D_i \|^2} of the orthogonal residuals from the fitted \code{\link{onls}} model.
}

\usage{
deviance_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\value{
The deviance of the orthogonal residuals.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## See 'onls'.
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

