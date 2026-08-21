\name{DiagnosticProb}
\docType{data}
\alias{DiagnosticProb}
\alias{Data17}
\title{Diagnostic probabilities}
\description{
Data record the probabilities assigned by subjective diagnostic of 15
clinicians and 15 statisticians.
}
\usage{
data(DiagnosticProb)
}
\source{
  Aitchison: CODA microcomputer statistical package, 1986, the file name DIAGPROB.DAT,
  here included under the GNU Public Library Licence Version 2 or newer.
}

\details{
The data consist of 30 cases: 15 diagnostics probabilities assigned by clinicians,
15 diagnostics probabilities assigned by statisticians, and
4 variables: probabilities A, B, and C, and type i.e. 1 for clinicians, 
2 for statisticians.
 
 In the study of subjective performance in inferential task the subject is faced
 with the finite set of mutually exclusive and exhaustive hypothesis, 
 and the basis of specific information presented to him/her is required to
 divide the available unit of probability among these probabilities.
 In this study the task is presented as a problem of differential diagnosis of
 three mutually exclusive and exhaustive diseases of students, known under the
 generic title of 'newmath syndrome',
 
 \tabular{rll}{
 \tab A \tab - algebritis,  \cr
 \tab B \tab - bilateral paralexia,  \cr
 \tab C \tab - calculus deficiency.  
 }
The subject, playing the role of diagnostician, is informed that the three diseases
types are equally common and is shown the results of 10 diagnostic tests on 60
previous cases of known diagnosis, 20 of each type. The subject is then shown the results
of the 10 tests for a new undiagnosed cases and asked to assign diagnostic probabilities
to the three possible disease types.

Data record the subjective assessments of 15
clinicians and 15 statisticians for the same case. For this case the objective
diagnosis probabilities are known to be $(.08, .05, .87).$

 All row probabilities sum to 1, except for some rounding errors.
}
\references{
   Aitchison: The Statistical Analysis of Compositional Data, 1986, Data 17, pp20.
 }
\note{Courtesy of J. Aitchison}
\keyword{datasets}
