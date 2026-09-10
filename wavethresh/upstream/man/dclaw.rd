\name{dclaw}
\alias{dclaw}
\alias{rclaw}
\alias{pclaw}
\title{Claw distribution}
\usage{
rclaw(n)
dclaw(x)
pclaw(q)
}
\arguments{
\item{n}{Number of draws from \code{rclaw} distribution}
\item{x}{Vector of ordinates}
\item{q}{Vector of quantiles}
}
\description{
Random generation, density and cumulative probability for the claw distribution.
}
\details{
The claw distribution is a normal mixture distribution, introduced in Marron & Wand (1992).
Marron, J.S. & Wand, M.P. (1992). Exact Mean Integrated Squared Error.
\emph{Ann. Stat.}, \bold{20}, 712--736.
}
\value{Random samples (rclaw), density (dclaw) or probability (pclaw)
of the claw distribution.
}
\examples{
# Plot the claw density on the interval [-3,3]
x <- seq(from=-3, to=3, length=500)
\dontrun{plot(x, dclaw(x), type="l")}
}
\author{David Herrick}
\keyword{smooth}
