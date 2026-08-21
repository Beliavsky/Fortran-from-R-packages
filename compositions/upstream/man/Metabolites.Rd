\name{Metabolites}
\docType{data}
\alias{Data09}
\alias{Metabolites}
\title{Steroid metabolite patterns in adults and children}
\description{
Data shows the urinary excretion (mg/24 hours) of 37 normal adults and 30 normal children
of
total cortisol meatbolites,
total corticosterone meatbolites,
total pregnanetriol and \eqn{\Delta}{Delta}-5-pregnentriol.
}
\usage{
data(Metabolites)
}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name METABOL.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}

\details{
There are 67 cases for 37 adults and 30 children, and 5 columns:
Case no., met1, met2, met3 and Type, 1 for adults, $-1$ for children.
No sum constraint is placed on this data set: since the urinary excretion 
in mg for 24 hours are given.
}

\references{
   Aitchison, J. (1986) The Statistical Analysis of Compositional Data, (Data 9), pp14.
 }
\note{Courtesy of J. Aitchison}
\keyword{datasets}
