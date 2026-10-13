\name{qcochran}
\alias{qcochran}
\alias{pcochran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Critical values and p-values for Cochran outlying variance test }
\description{

This functions calculates quantiles (critical values) and reversively p-values
for Cochran test for outlying variance.

}
\usage{
qcochran(p, n, k)
pcochran(q, n, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of probabilities. }
  \item{q}{ vector of quantiles. }
  \item{n}{ number of values in each group (if not equal, use arithmetic mean). }
  \item{k}{ number of groups.}
}
\value{

Vector of p-values or critical values.

}
\references{

Snedecor, G.W., Cochran, W.G. (1980). Statistical Methods (seventh edition).
		 Iowa State University Press, Ames, Iowa.


}
\author{ Lukasz Komsta }
\seealso{ \code{\link{cochran.test}} }
\examples{
qcochran(0.05,5,5)
pcochran(0.293,5,5)
}
\keyword{distribution}
