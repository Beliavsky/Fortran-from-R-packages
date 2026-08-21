\name{Sediments}
\docType{data}
\alias{Data20}
\alias{Sediments}
\title{Proportions of sand, silt and clay in sediments specimens}
\description{
Data provide sand, silt and clay compositions of 21  sediments specimens,
10 of which are identified as \emph{offshore,} 7 as \emph{near shore} and 
4 new samples. 
}
\usage{
data(Sediments)
}
\format{
  \describe{
    \item{sand}{numeric: the portion of sand}
    \item{silt}{numeric: the protion of silt}
    \item{clay}{numeric: the portion of clay}
    \item{type}{numeric: 1 for offshore, 2 for near shore and 3 for new
      samples}
  }
}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name YATQUAD.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}
\details{
The data comprise 21 cases: 10 offshore, 7 near shore and 4 new samples,
and 4 variables: sand, silt and clay proportions and in addition the type of  
sediments specimens -- 1 for offshore, 2 for near shore and 3 for new samples. 

All 3-part compositions sum to one.
}

\references{
   Aitchison, J. (1986) The Statistical Analysis of Compositional Data, (Data 20) pp17. 
 }
\note{Courtesy of J. Aitchison}
\keyword{datasets}