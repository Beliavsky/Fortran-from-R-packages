\name{WhiteCells}
\docType{data}
\alias{Data11}
\alias{WhiteCells}
\title{White-cell composition of 30 blood samples by two different methods}
\description{
In 30 blood samples portions of three kinds of white cells
\describe{
  \item{G:}{granulocytes,}
  \item{L:}{lymphocytes,}
  \item{M:}{monocytes,}
}
were determined with two methods, time-consuming microscopic and
automatic image analysis. The resulting 30 pairs of 3-part compositions
are recorded.
}
\usage{
data(WhiteCells)
}
\format{A 30x6  matrix}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name WCELLS.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}
\details{
In an experiment each of 30 blood samples was halved, one half being assigned
randomly to one method, the other half to the other method.
We have 60 cases of 3-part compositions but these are
essentially 30 pairs of related compositions. All 3-part portions sums to one, 
except for some rounding errors.
}
\note{Courtesy of J. Aitchison}
\references{
   Aitchison, J. (1986) The Statistical Analysis of Compositional Data, 1986 (Data 9) pp16. 
 }
\keyword{datasets}
