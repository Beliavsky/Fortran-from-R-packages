\name{video}
\alias{video}
\docType{data}
\title{Different raters judging the credibility of videotaped testimonies}
\description{
The data frame contains the credibility ratings of 20 subjects, rated by 4 raters.
Judgements could vary from 1 (not credible) to 6 (highly credible). Variance between and within raters is low.
}
\usage{data(video)}
\format{
  A data frame with 20 observations on the following 4 variables.
  \describe{
    \item{rater1}{ratings of rater 1}
    \item{rater2}{ratings of rater 2}
    \item{rater3}{ratings of rater 3}
    \item{rater4}{ratings of rater 4}
  }
}
\source{
artificial data
}
\examples{
data(video)
apply(video,2,table)
}
\keyword{datasets}
