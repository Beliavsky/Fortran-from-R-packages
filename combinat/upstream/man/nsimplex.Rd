\name{nsimplex}
\alias{nsimplex}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Computes the number of points on a (p, n)-simplex lattice }
\description{
Computes the number of points on a (p, n)-simplex lattice; that is, the
number of p-part compositions of n. This gives the number of points in
the support space of a Multinomial(n, q) distribution, where
p == length(q).

Arguments p and n are replicated as necessary to have the length of the
longer of them.
}
\usage{
nsimplex(p, n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{p}{ vector of integers }
  \item{n}{ vector of integers }
}
%\details{
%}
\value{
integer
}
\examples{
nsimplex(3,5)
}
\keyword{ models }

