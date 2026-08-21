\name{Firework}
\docType{data}
\alias{Firework}
\alias{Data13}
\title{Firework mixtures}
\description{
Data show two measured properties, brilliance and vorticity, of 81 girandoles
composed of different mixtures of five ingredients: a -- e.
Of these ingredients, a and b are the primary light-producting, 
c is principal propellant, and d and e are binding agents for c. 
}
\usage{
data(Firework)
}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name YATQUAD.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}

\details{
The data consist of 81 cases and 7 variables: ingredients a, b, c, d, and e, 
and the two measured properties \emph{brilliance} and \emph{vorticity}.
The 81 different mixtures form a special experiment design. First the 81 
possible quadruples formed from the three values -1, 0, 1 were arranged 
in ascending order. Then for each such quadruple z, the corresponding mixture
x(z)=(a,b,c,d,e)=alrInv(z) is computed. Thus the No. 4 girandole corresponds 
to z=(-1,-1,0,-1) and so is composed of a mixture x=(.12,.12,.32,.12,.32) 
of the five ingredients. All 5-part mixtures sum up to one.}

\references{
   Aitchison: The Statistical Analysis of Compositional Data, 1986, Data 13, pp17.  
 }
\note{Courtesy of J. Aitchison}
\keyword{datasets}

