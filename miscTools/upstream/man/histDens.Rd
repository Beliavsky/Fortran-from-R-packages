\name{histDens}
\alias{histDens}
\title{Histogram with a Density Line}

\description{
   Plot a histrogram and add a kernel density line.
}

\usage{
histDens( x, breaks = "Sturges", ... )
}

\arguments{
   \item{x}{values of the variable.}
   \item{breaks}{passed to \code{\link[graphics]{hist}}.}
   \item{\dots}{further arguments are passed to \code{\link[graphics]{hist}}.}
}

\author{Arne Henningsen}

\examples{
   set.seed( 123 )
   x <- rnorm( 100 )
   histDens( x )
   histDens( x, 20 )
   histDens( x, 20, main = "My Title" )
}

\keyword{models}
