\name{Lasym}
\alias{Lasym}

\title{Lorenz Asymmetry Coefficient}

\description{
Coefficient of asymmetry in the Lorenz curve.
}

\usage{
Lasym(x, n = rep(1, length(x)), interval = FALSE, na.rm = TRUE)
}

\arguments{
 \item{x}{a vector containing non-negative elements.}
 \item{n}{a vector of frequencies, must be same length as \code{x}.}
 \item{interval}{logical. In the case where there are observations
   exactly equal to the mean, either an interval of asymmetry coefficients
   can be returned or their midpoint.}
 \item{na.rm}{logical. Should missing values (\code{NA}s) be removed
   prior to computations? If set to \code{FALSE} the computations yield
   \code{NA}.}
}

\details{
Damgaard and Weinter (2000) have suggested an additional measure for
comparing inequality in in distributions (specifically for describing
plant size or fecundity distributions) to accompany the Lorenz curve
and Gini coefficient. It assesses the asymmetry in the Lorenz curve
of the distributions.
}

\references{
C Damgaard, J Weiner: Describing Inequality in Plant Size or Fecundity, 2000.
Ecology 81(4), 1139--1142.
}

\seealso{\code{\link{Lc}}, \code{\link{Gini}}}

\examples{
## Examples from Damgaard & Weiner (2000)

## Figure 2
x <- rep(c(50/9, 50), c(9, 1))
y <- rep(c(2, 18), c(5, 5))
plot(table(x))
plot(table(y))

## statistics
mean(x)
mean(y)
Gini(x, corr = TRUE)
Gini(y, corr = TRUE)
Lasym(x)
Lasym(y)

## Figure 3
plot(Lc(x))
lines(Lc(y), col = "slategray")
abline(1, -1, lty = 2)

}

\keyword{misc}
