\name{plot.loadings}
\alias{plot.loadings}
\title{ Plot loadings }
\description{
 Plots loadings as a separate barplot for each factor.
}
\usage{
\method{plot}{loadings}(x, n = 5, k = ncol(x), mfrow = c(k, 1), ...)
}
\arguments{
  \item{x}{A \code{loadings} object.}
  \item{n}{Number of components of each factor to plot.}
  \item{k}{Number of factors to plot.}
  \item{mfrow}{Passed to \code{par(mfrow = ...)} if \code{k > 1}.}
  \item{\dots}{Other arguments passed to \code{barplot}.}
}
\details{
  The top \code{n} components of each of the top \code{k} factors
  are displayed in a separate \code{barplot}.  
}
\value{
  Return value is a list of the return values from each \code{barplot}
  invocation.
}
\seealso{\code{\link{barplot}}}
\examples{
data(m.barra.9003)
rtn <- m.barra.9003
stat.fac <- factanal(rtn, factors = 3)
m.barra.loadings <- loadings(stat.fac)
plot(m.barra.loadings)
}
\keyword{hplot}
