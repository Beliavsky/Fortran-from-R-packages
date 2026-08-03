\name{hl.loc}
\alias{hl.loc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hodges - Lehmann Estimator of Location}
\description{
  Function to compute the Hodges - Lehmann estimator of location in the one sample case.
}
\usage{
hl.loc(x, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  The Hodges - Lehmann estimator is the median of the combined data points and Walsh averages.
  It is the same as the Pseudo Median returned as a by-product of the function \code{\link{wilcox.test}}.
}
\value{
  the Hodges - Lehmann estimator of location.
}
\references{ \cite{Hettmansperger, T.P. and McKean, J.W. (1998), \emph{Robust Nonparametric Statistical Methods},
London, Arnold. }

\cite{Hodges, J.L., and Lehmann, E.L. (1963), Estimates of location based on rank tests. \emph{The Annals of Mathematical Statistics}, \bold{34}, 598--611. }
}


\author{ Klaus Nordhausen }

\seealso{\code{\link{wilcox.test}}}
\examples{
set.seed(1)
x <- rt(100, df = 3)
hl.loc(x)
# same as
wilcox.test(x,  conf.int = TRUE)$estimate
rm(.Random.seed)
}
\keyword{nonparametric}
\keyword{robust}
\keyword{univar}
