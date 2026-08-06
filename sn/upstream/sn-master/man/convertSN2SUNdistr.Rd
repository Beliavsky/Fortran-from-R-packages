\name{convertSN2SUNdistr}
\alias{convertSN2SUNdistr}
 
\title{Convert a SN distribution into a SUN}

\description{An object of \code{SECdistrMv-class} or \code{SECdistrUv-class} 
  representing a \acronym{SN} or \acronym{ESN}  distribution is converted  
  into a \code{SUNdistr-class} object representing the same distribution.}

\usage{convertSN2SUNdistr(object, HcompNames = "h", silent = FALSE)}

\arguments{
  \item{object}{an object of \code{SECdistrMv-class} with \code{family} of type
    \code{SN} or \code{ESN}.}
  \item{HcompNames}{an optional character string for the hidden component} 
  \item{silent}{a logical value which controls the behaviour if the supplied
    \code{object} is not suitable. If \code{silent = FALSE}  (detault value)
    an error message is generated, otherwise a \code{NULL} is silently
    returned.}
}

%\details{%%  ~~ If necessary, more details than the description above ~~}

\value{an object of \code{SUNdistr-class}}

%\references{%% ~put references to the literature/web site here ~}

\author{Adelchi Azzalini}

% \note{%%  ~~further notes~~}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{SUNdistr-class}}, \code{\link{SECdistrMv-class}}, 
  \code{\link{SECdistrUv-class}}}

\examples{
esn <- makeSECdistr(dp=c(0, 1, 2, 0.5), family="ESN")
sun <- convertSN2SUNdistr(esn)
mean(sun) - mean(esn)
vcov(sun) - sd(esn)^2
#
dp0 <- list(xi=1:2, Omega=diag(3:4), alpha=c(3, -5))
f10 <- makeSECdistr(dp=dp0, family="SN", name="SN-2d", compNames=c("u1", "u2"))
sun10 <- convertSN2SUNdistr(f10)
mean(sun10) - mean(f10)
vcov(sun10) - vcov(f10)
}

\keyword{distribution}
\keyword{multivariate} 
\concept{SUN distribution}
\concept{Unified Skew-Normal distribution}
