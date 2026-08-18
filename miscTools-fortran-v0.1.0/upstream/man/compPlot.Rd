\name{compPlot}
\alias{compPlot}
\title{Scatterplot to Compare two Variables}

\description{
   Plot a scatterplot to compare two variables.
}

\usage{
compPlot( x, y, lim = NULL, ... )
}

\arguments{
   \item{x}{values of the first variable (on the X axis).}
   \item{y}{values of the second variable (on the Y axis).}
   \item{lim}{optional vector of two elements specifying the limits
      of both axes).}
   \item{\dots}{further arguments are passed to \code{\link[graphics]{plot}}.}
}

\author{Arne Henningsen}

\examples{
   set.seed( 123 )
   x <- runif( 25 )
   y <- 2 + 3 * x + rnorm( 25 )
   ols <- lm( y ~ x )
   
   compPlot( y, fitted( ols ) )
   compPlot( y, fitted( ols ), lim = c( 0, 10 ) )
   compPlot( y, fitted( ols ), pch = 20 )
   compPlot( y, fitted( ols ), xlab = "observed", ylab = "fitted" )
   compPlot( y, fitted( ols ), log = "xy" )
}

\keyword{models}
