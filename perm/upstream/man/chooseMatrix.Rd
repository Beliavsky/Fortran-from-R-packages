\name{chooseMatrix}
\alias{chooseMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create matrix with choose(n,m) unique rows with m ones in each row the rest 0.}
\description{
Create a \code{choose(n,m)} by n matrix. The matrix has unique rows with m ones in each row and the rest zeros.
}
\usage{
chooseMatrix(n, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{an integer}
  \item{m}{an integer<=n}
}
%\details{}
\value{
A matrix with choose(n,m) rows n columns. 
The matrix has unique rows with m ones in each row and the rest zeros.
}
%\references{ ~put references to the literature/web site here ~ }
\author{M.P.Fay}
\note{Used for complete enumeration when method='exact.ce' in \code{\link{permTS}}}
\seealso{ \code{\link{permTS}} }
\examples{
chooseMatrix(5,2)
}
\keyword{array}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
