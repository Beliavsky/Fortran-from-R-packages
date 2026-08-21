\name{eeg}
\alias{eeg}
\docType{data}
\title{EEG Data}
\description{
An electroencephalogram (EEG) is a noninvasive test used to detect
and record the electrical activity generated in the brain. These data
were measured at a frequency of 256 per second and came from a
patient suffering a seizure. This a portion of a series on the website
of Professor Richard Smith, University of North Carolina. His
source: Professors Mike West and Andrew Krystal, Duke University.}
\usage{data(eeg)}
\format{
  The format is:
 ts [1:13000, 1]  -3.08 -20.15 -45.05 -69.95 -94.57 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "eeg"
 - attr(*, "tsp")= num [1:3]  2001 15000     1
}
\source{
\code{http://www.stat.unc.edu/faculty/rs/s133/Data/datadoc.html}}
\examples{
data(eeg)
## maybe str(eeg) ; plot(eeg) ...
}
\keyword{datasets}
