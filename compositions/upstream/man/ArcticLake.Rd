\name{ArcticLake}
\docType{data}
\alias{ArcticLake}
\alias{Data05}
\title{Artic lake sediment samples of different water depth}
\description{
 Sand, silt and clay compositions of 39 sediment samples of different water depth
 in an Arctic lake. 
 }
\usage{
data(ArcticLake)
}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name ARCTIC.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}

\details{
Sand, silt and clay compositions of 39 sediment samples at different water depth (in meters)
 in an Arctic lake. The additional feature is a concomitant variable or \emph{covariate,}
 water depth, which may account for some of the variation in the compositions.
 In statistical terminology we have a multivariate regression problem with sediment composition 
 as regressand and water depth as regressor.
 
 All row percentage sums to 100, except for rounding errors.
}

\references{
   Aitchison: The Statistical Analysis of Compositional Data, 1986, Data 5, pp5.
 }
\note{Courtesy of J. Aitchison}
\keyword{datasets}
