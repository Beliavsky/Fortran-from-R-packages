\name{powm}
\alias{powm}
\title{Exponentiation function}
\description{
  This function return \eqn{x ^ y mod n}.

  This function return \eqn{x ^ y mod n}{x ^ y mod n}
  pow.bigz do the same when modulus is set.

}
\usage{
powm(x, y, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Integer or big integer - possibly a vector}
  \item{y}{Integer or big integer - possibly a vector}
  \item{n}{Integer or big integer - possibly a vector}
}
% \details{
%   This function return \eqn{x ^ y mod n}{x ^ y mod n}

%   pow.bigz do the same when modulus is set.
% }
\value{
  A bigz class representing the parameter value.
}

\author{A. L.}

\seealso{\code{\link{pow.bigz}}}
\examples{
powm(4,7,9)


x = as.bigz(4,9)
x ^ 7

}


\keyword{arith}

