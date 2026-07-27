\name{interpol1d}
\alias{interpol1d}
\alias{interpolation1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear interpolation function. 
}
\description{ Interpolates linearly a function given at equally spaced points on the interval [minS,maxS].}

\usage{
interpol1d(x,F0,minS,maxS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{Point at which the function is interpolated.}
\item{minS}{Minimum value of the grid.}
\item{maxS}{Maximum value of the grid.}
\item{F0}{Value of the function at m equally spaced points on the grid.}
}


\value{
\item{interpol}{Linear interpolation of the function at point x.}
}
\references{
Chapter 3 of 'Statistical Methods for Financial Engineering, B. Remillard, CRC Press, (2013).
}
\author{
Bruno Remillard
}


\examples{
F0 = c(1:10)
minS = 1;
maxS = 10;
out = interpol1d(2.45,F0,1,10)
out #since the function is the identity, the answer should be 2.45!
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linear interpolation}
% __ONLY ONE__ keyword per line
