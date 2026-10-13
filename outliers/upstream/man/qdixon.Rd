\name{qdixon}
\alias{qdixon}
\alias{pdixon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ critical values and p-values for Dixon tests}
\description{

Approximated quantiles (critical values) and distribution function (giving p-values)
for Dixon tests for outliers.

}
\usage{
qdixon(p, n, type = 10, rev = FALSE)
pdixon(q, n, type = 10)
}
\arguments{
  \item{p}{ vector of probabilities. }
  \item{q}{ vector of quantiles. }
  \item{n}{ length of sample. }
  \item{type}{ integer value: 10, 11, 12, 20, or 21. For description see \code{dixon.test}. }
  \item{rev}{ function \code{qdixon} with this parameter set to TRUE acts as \code{pdixon}.}
}
\details{

This function is based on tabularized Dixon distribution, given by Dixon (1950) and corrected
by Rorabacher (1991). Continuity is reached due to smart interpolation using \code{\link{qtable}} function.
By now, numerical procedure to obtain these values for n>3 is not known.

}
\value{

Critical value or p-value (vector).

}
\references{ 
Dixon, W.J. (1950). Analysis of extreme values. Ann. Math. Stat. 21, 4, 488-506.

Dixon, W.J. (1951). Ratios involving extreme values. Ann. Math. Stat. 22, 1, 68-78.

Rorabacher, D.B. (1991). Statistical Treatment for Rejection of Deviant Values: Critical Values of Dixon Q Parameter and Related Subrange Ratios at the 95 percent Confidence Level. Anal. Chem. 83, 2, 139-146.

}

\author{ Lukasz Komsta }
\seealso{ \code{\link{qtable}}, \code{\link{dixon.test}}}
\keyword{htest}

