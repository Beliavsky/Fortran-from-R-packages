\name{hcube}
\alias{hcube}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generate all points on a hypercuboid lattice. }
\description{
Generate all points on a hypercuboid lattice.
}
\usage{
hcube(x, scale, translation)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ 
         Argument x is an integer vector giving the extent of each dimension;
         the number of dimensions is length(x).
}
  \item{scale}{ 
         Argument scale is a vector of real numbers giving an amount by which
         to multiply the points in each dimension;  it will be replicated as
         necessary to have the same length as x.
 }
  \item{translation}{ 
         Argument translate is a vector of real numbers giving an amount to
         translate (from the "origin", rep(1,length(x))) the points in each
         dimension;  it will be replicated as necessary to have the same
         length as x.  To use rep(0,length(x)) as the origin,  use
         translation = -1.  Scaling,  if any,  is done BEFORE translation.
 }
}
%\details{
%}
\value{
         A prod(x) by length(x) numeric matrix;  element (i,j) gives the
         location of point i in the jth dimension.  The first column
         (dimension) varies most rapidly.

}
\references{ ~put references to the literature/web site here ~ }
\author{ code by Scott Chasalow, R pack and maint by VJ Carey <stvjc@channing.harvard.edu>}

\seealso{ 
fac.design, expand.grid }

%\examples{
%}
%\\keyword{ combinatorics }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line

