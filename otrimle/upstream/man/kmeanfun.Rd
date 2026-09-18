\name{kmeanfun}

\alias{kmeanfun}
\alias{ksdfun}


\title{Mean and standard deviation of unimodality statistic}

\description{
  These functions approximate the mean and standard deviation of the
  unimodality statistic computed by \code{\link{kerndensmeasure}}
  assuming standard Gaussian data dependent on the number of
  observations\code{n}. They have been chosen based on a simulation involving 74
  different values of \code{n}. Used for standardisation in
  \code{\link{kerndensp}}.
} 

\usage{
  kmeanfun(n)
  ksdfun(n)
}


\arguments{
  \item{n}{integer. Number of observations.}
}

\value{
  The resulting mean (\code{kmeanfun}) or standard deviation
  (\code{ksdfun}).
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{kerndensp}} 
}

\examples{
   kmeanfun(50)
   ksdfun(50)
}

\keyword{htest}


