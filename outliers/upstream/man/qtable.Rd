\name{qtable}
\alias{qtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Interpolate tabularized distribution }
\description{

This function calculates critical values or p-values which cannot
be obtained numerically, and only tabularized version is available.

}
\usage{
qtable(p, probs, quants)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of probabilities. }
  \item{probs}{ vector of given probabilities. }
  \item{quants}{ vector of given corresponding quantiles. }
}
\details{

This function is internal routine used to obtain Grubbs and Dixon critical values.
It fits linear or cubical regression to closests values of its argument, then uses
obtained function to obtain quantile by interpolation.

}
\value{

A vector of interpolated values

}
\author{ Lukasz Komsta }
\note{ 

You can simply do "reverse" interpolation (p-value calculating) by reversing probabilities
and quantiles (2 and 3 argument).

}
\keyword{distribution}
