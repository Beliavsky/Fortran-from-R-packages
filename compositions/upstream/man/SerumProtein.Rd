\name{SerumProtein}
\docType{data}
\alias{Data16}
\alias{SerumProtein}
\title{Serum Protein compositions of blood samples}
\description{
Data recording the proportions of the 4 serum proteins from blood samples
of 30 patients, 14 with known disease A, 16 with known disease B, and
6 new cases.
}
\usage{
data(SerumProtein)
}
\format{
  \describe{
    \item{a}{numeric a protein type}
    \item{b}{numeric a protein type}
    \item{c}{numeric a protein type}
    \item{d}{numeric a protein type}
    \item{Type}{1 deasease A, 2 disease B, 3 new cases}
  }
}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name SERPROT.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}
\details{
The data consist of 36 cases: 14 with known disease A, 16 with known disease B, and
6 new cases and 5 v variables: a, b, c, and d for 4 serum proteins and Type
for the diseases: 1 for  disease A, 2 for  disease B, and 3 for new cases. 
 All row serum proteins proportions sums to 1 except some rounding errors.
}
\references{
   Aitchison, J. (1986) The Statistical Analysis of Compositional Data (Data 16) pp20.
 }
\note{Courtesy of J. Aitchison}
\keyword{datasets}
