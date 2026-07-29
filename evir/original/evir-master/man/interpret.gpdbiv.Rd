\name{interpret.gpdbiv}
\alias{interpret.gpdbiv}
\title{Interpret Results of Bivariate GPD Fit}
\description{
  Interprets the results of a bivariate GPD model fitted
  using the bivariate POT method. 
}
\usage{
interpret.gpdbiv(out, x, y)
}
\arguments{
\item{out}{a \code{gpdbiv} object}
\item{x}{a scalar value greater than first threshold}
\item{y}{a scalar value greater than second threshold}
}
\value{
 A vector of probabilities is invisibly returned, in printed order. 
}
\section{Side Effects}{
  A simple interpretation of the fit in terms of exceedance 
  probabilities for the point (x,y) is printed. 
}
\details{
  First marginal probabilities of exceeding the points x and y
  are calculated, and then joint and conditional probabilities. 
}
\seealso{\code{\link{gpdbiv}}, \code{\link{plot.gpdbiv}}}

\examples{
data(bmw) ; data(siemens)
out <- gpdbiv(-bmw, -siemens, ne1 = 100, ne2 = 100) 
interpret.gpdbiv(out, 0.05, 0.05) 
# probabilities of 5\% falls in BMW and Siemens stock prices 
}
\keyword{htest}
