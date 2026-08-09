\name{logLik_o}
\alias{logLik_o}
\encoding{latin1}

\title{Log-Likelihood for the orthogonal residuals}

\description{
Returns the log-likelihood as calculated from the orthogonal residuals obtained by \code{\link{residuals_o}}.
}

\usage{
logLik_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\value{
The log-likelihood as calculated from the orthogonal residuals.
}

\note{
\code{logLik_o} has no other generic functions on top, so for calculating \code{\link{AIC}}, one has to apply \code{AIC(logLik_o(model))}, see 'Examples'. The usual \code{\link{logLik}} applies to the vertical residuals of the orthogonal model.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
logLik_o(mod)

## compare AIC of vertical versus orthogonal residuals.
AIC(mod)
AIC(logLik_o(mod))
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

