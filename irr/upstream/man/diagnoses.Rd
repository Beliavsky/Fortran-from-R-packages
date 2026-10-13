\name{diagnoses}
\alias{diagnoses}
\docType{data}
\title{Psychiatric diagnoses provided by different raters}
\description{
Psychiatric diagnoses of n=30 patients provided by different sets of m=6 raters.
Data were used by Fleiss (1971) to illustrate the computation of Kappa for m raters.
}
\usage{data(diagnoses)}
\format{
  A data frame with 30 observations (psychiatric diagnoses with levels 1. Depression, 2. Personality Disorder, 3. Schizophrenia, 4. Neurosis, 5. Other) on 6 variables representing different raters.
  \describe{
    \item{rater1}{a factor including the diagnoses of rater 1 (levels see above)}
    \item{rater2}{a factor including the diagnoses of rater 2 (levels see above)}
    \item{rater3}{a factor including the diagnoses of rater 3 (levels see above)}
    \item{rater4}{a factor including the diagnoses of rater 4 (levels see above)}
    \item{rater5}{a factor including the diagnoses of rater 5 (levels see above)}
    \item{rater6}{a factor including the diagnoses of rater 6 (levels see above)}
  }
}
\source{
Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin, 76, 378-382.
}
\references{
Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin, 76, 378-382.
}
\examples{
data(diagnoses)
table(diagnoses[,1])
}
\keyword{datasets}
