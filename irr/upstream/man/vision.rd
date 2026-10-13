\name{vision}
\alias{vision}
\docType{data}
\title{Eye-testing case records}
\description{
Case records of the eye-testing of N=7477 female employees in Royal Ordnance factories 
between 1943 and 1946. Data were primarily used by Stuart (1953) to illustrate the 
the estimation and comparison of strengths of association in contingency tables.
}
\usage{data(anxiety)}
\format{
  A data frame with 7477 observations (eye testing results with levels
  1st grade, 2nd grade, 3rd grade, 4th Grade) on the following 2 variables.
  \describe{
    \item{r.eye}{unaided distance vision performance of the right eye}
    \item{l.eye}{unaided distance vision performance of the left eye}
  }
}
\source{
Stuart, A. (1953). The Estimation and Comparison of Strengths of Association in 
Contingency Tables. Biometrika, 40, 105-110.
}
\references{
Stuart, A. (1953). The Estimation and Comparison of Strengths of Association in 
Contingency Tables. Biometrika, 40, 105-110.
}
\examples{
data(vision)
table(vision$r.eye, vision$l.eye)
}
\keyword{datasets}
