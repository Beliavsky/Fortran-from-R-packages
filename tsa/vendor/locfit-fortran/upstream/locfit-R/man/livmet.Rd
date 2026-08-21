\name{livmet}
\alias{livmet}
\title{liver Metastases dataset}
\usage{data(livmet)}
\format{
  Data frame with survival times (\code{t}), censoring indicator
  (\code{z}) and a number of covariates.
}
\description{
  Survival times for 622 patients diagnosed with Liver Metastases.

  Beware, the censoring variable
  is coded as 1 = uncensored, so use \code{cens=1-z} in
  \code{\link{locfit}()} calls.
}
\source{
  Haupt and Mansmann (1995)
}
\references{
  Haupt, G. and Mansmann, U. (1995)
  CART for Survival Data.
  Statlib Archive.
}
\keyword{datasets}
