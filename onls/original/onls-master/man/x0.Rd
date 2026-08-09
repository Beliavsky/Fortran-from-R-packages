\name{x0}
\alias{x0}
\alias{y0}
\encoding{latin1}

\title{x0/y0-values from orthogonal nonlinear least squares regression}

\description{
Returns the \eqn{{\color{red} x_{0i}/y_{0i}}} values as obtained from minimizing the Euclidean distance \deqn{\min\| D_i \| \equiv \min\sqrt{(x_i - {\color{red}x_{0i}})^2 + (y_i - {\color{red}y_{0i}})^2}} 
}

\usage{
x0(object) 
y0(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\value{
A vector of \eqn{\color{red}x_{0i}} or \eqn{\color{red}y_{0i}} values.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
x0(mod)
y0(mod)
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

