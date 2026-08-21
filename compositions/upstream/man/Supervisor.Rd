\name{Supervisor}
\docType{data}
\alias{Data07}
\alias{Supervisor}
\title{Proportions of supervisor's statements assigned to different categories}
\description{
The results of a study of a single supervisor in his relationship to three supervisee
are recorded. Instructions in a technical subject took place in sessions of one hour 
and with only one supervisee at the time. Each supervisee attended six sessions
(once every two weeks in a twelve-week period).
All of 18 sessions were recorded and for each session the 'statements'
of the supervisor were classified into four categories.
Thus for each session the proportion of statements in the four categories
are set out in a two-way table according to the fortnight (6) and the supervisee (3).  
}
\usage{
data(Supervisor)
}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name SUPERVIS.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}
\format{A 18x13  matrix}
\details{
For each session the 'statements' of the supervisor were classified into
four categories
\describe{
  \item{C:}{commanding, posing a specific instruction to the supervisee,}
  \item{D:}{demanding, posing a specific question to the supervisee,}
  \item{E:}{exposing, providing the supervisee with an explanation,}
  \item{F:}{faulting, pointing out faulty technique to the supervisee.}
}
Thus for each session the proportion of statements in the four categories
are set out in a two-way table according to the fortnight (6) and the supervisee (3).
The  C, D, E, F values in the rows sum mostly to 1, except for some rounding errors.
}
\note{Courtesy of J. Aitchison}
\references{
   Aitchison, J. (1986): The Statistical Analysis of Compositional Data, (Data 7) pp12.
 }
\keyword{datasets}

