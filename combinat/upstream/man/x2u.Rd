\name{x2u}
\alias{x2u}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{  Convert an x-encoded simplex-lattice point to a u-encoded simplex-lattice point }
\description{
Convert an x-encoded simplex-lattice point to a u-encoded
simplex-lattice point  (equivalently,  "untabulate" bin counts)
}
\usage{
x2u(x, labels=seq(along = x))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ 
x:    A numeric vector.  x[i] is interpreted as the count in bin i.}
  \item{labels}{ 
A vector.  Interpreted as the bin labels;  default value is
        seq(along = x), which causes return of a u-encoded simplex-lattice
        point.  Other values of labels cause return of the result of
        subscripting labels with the u-encoded simplex-lattice point that
        would have been obtained if the default value of labels were used.
 }
}
%\details{
%}
\value{
        rep(labels, x), a vector of length sum(x).  If labels = seq(along = x)
        (the default),  value is the u-encoded translation of the simplex
        lattice point, x.  Equivalently,  value gives the bin numbers,
        in lexicographic order,  for the objects represented by the counts in
        x.  For other values of argument "labels", value gives the bin labels
        for the objects represented by the counts in x (equivalent to
        labels[x2u(x)]).
}
%\references{ }

\seealso{ tabulate, rep }

%\examples{
%}
\keyword{ models }

