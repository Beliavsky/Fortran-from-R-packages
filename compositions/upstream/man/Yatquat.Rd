\name{Yatquat}
\docType{data}
\alias{Data12}
\alias{Yatquat}
\title{Yatquat fruit evaluation}
\description{
The quality of yatquat tree fruit is assessed in terms of relative proportions 
by volume of flesh, skin and stone.
In an experiment an arboriculturist uses 40 trees, randomly allocated 20
to the hormone treatment and leaves untreated the remaining 20 trees.
Data provides fruit compositions of the present season and the preceding season, 
as well as the treatment: 1 for the treated trees, -1 for untreated trees.
}
\usage{
data(Yatquat)
}
\format{A 40x7 data matrix}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name YATQUAT.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}
\details{
The yatquat tree produces each season a single large fruit.
Data provides fruit compositions of the present season, the compositions of the fruit
of the same 40 trees for the preceding season when none of the trees were treated, 
and in addition the Type: 1 for the treated trees, -1 for untreated trees. 
For each of the 40 cases we have two 3-part composition on flesh, skin and stone.
The column names are:
\tabular{rll}{
 \tab prFL\tab  portion of fruit flesh in the present season,\cr
 \tab prSK\tab  portion of fruit skin in the present season,\cr
 \tab prST\tab  portion of fruit stone in the present season,\cr
 \tab Type\tab  1 for treated,  $-1$ for untreated trees,\cr
 \tab paFL\tab  portion of fruit flesh in the preceding season,\cr
 \tab paSK\tab  portion of fruit skin in the preceding season,\cr
 \tab paST\tab  portion of fruit stone in the preceding season,
}
All 3-part compositions sum to one.
}
\examples{
#data(Yatquat)
#plot(acomp(Yatquat[,1:3]),col=as.numeric(Yatquat[,4])+2)
#plot(acomp(Yatquat[,5:7]),col=as.numeric(Yatquat[,4])+2)
}
\note{Courtesy of J. Aitchison}
\references{
   Aitchison, J. (1986) The Statistical Analysis of Compositional Data, (Data 12) pp17. 
 }
\keyword{datasets}
