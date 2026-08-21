\name{PlotsML}
\Rdversion{2.1}
\alias{PlotsML}
\alias{plotsML.gb2}
\alias{saveplot}
\title{
Cumulative Distribution Plot and Kernel Density Plot for the Fitted GB2
}
\description{
Function \code{plotsML.gb2} produces two plots. The first is a plot of the empirical cumulative distribution function versus the fitted cumulative distibution function.  
The second is a plot of the kernel density versus the fitted GB2 density.  Function \code{saveplot} saves locally the produced plot.
}
\usage{
plotsML.gb2(z, shape1, scale, shape2, shape3, w=rep(1,length(z))) 
saveplot(name, pathout)
}
\arguments{
  \item{z}{numeric; vector of data values.}
  \item{w}{numeric; vector of weights. Must have the same length as \code{z}. By default \code{w} is a vector of 1.}
  \item{shape1}{numeric; positive parameter.}
  \item{scale}{numeric; positive parameter.}
  \item{shape2, shape3}{numeric; positive parameters of the Beta distribution.}
  \item{name}{string; the name of the plot.}
  \item{pathout}{string; the path of the folder or device where the plot will be saved.}
}
\details{
The used kernel is "Epanechnikov" (see \code{\link{plot}}).  
}
\author{
Monique Graf and Desislava Nedyalkova
}
\keyword{distribution}

