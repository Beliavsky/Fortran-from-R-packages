\name{anxiety}
\alias{anxiety}
\docType{data}
\title{Anxiety ratings by different raters}
\description{
The data frame contains the anxiety ratings of 20 subjects, rated by 3 raters.
Values are ranging from 1 (not anxious at all) to 6 (extremely anxious).
}
\usage{data(anxiety)}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{rater1}{ratings of the first rater}
    \item{rater2}{ratings of the second rater}
    \item{rater3}{ratings of the third rater}
  }
}
\source{
artificial data
}
\examples{
data(anxiety)
apply(anxiety,2,table)
}
\keyword{datasets}
