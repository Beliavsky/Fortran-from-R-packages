\name{iota}
\alias{iota}
\title{iota coefficient for the interrater agreement of multivariate observations}
\description{
Computes iota as an index of interrater agreement of quantitative or nominal multivariate observations.
}
\usage{
iota(ratings, scaledata = c("quantitative","nominal"),
     standardize = FALSE)
}
\arguments{
  \item{ratings}{list of n*m matrices or dataframes with one list element for each variable, n subjects m raters.}
  \item{scaledata}{a character string specifying if the data is '"quantitative"' (default) or '"nominal"'. If the data is organized in factors, '"nominal"' is chosen automatically. You can specify just the initial letter.}
  \item{standardize}{a logical indicating whether quantitative data should be z-standardized within each variable before the computation of iota.}
}
\details{
Each list element must contain observations for each rater and subject without missing values.\cr
In case of one categorical variable (only one list element), iota reduces to the Fleiss exact kappa coefficient, which was proposed by Conger (1980).
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{value of iota.}
  \item{$detail}{a character string specifying if the values were z-standardized before the computation of iota.}
}
\references{
Conger, A.J. (1980). Integration and generalisation of Kappas for multiple raters. Psychological Bulletin, 88, 322-328.\cr\cr
Janson, H., & Olsson, U. (2001). A measure of agreement for interval or nominal multivariate observations. Educational and Psychological Measurement, 61, 277-289.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{icc}},
  \code{\link[irr]{kappam.fleiss}}
}
\examples{
data(diagnoses)
iota(list(diagnoses))  # produces the same result as...
kappam.fleiss(diagnoses, exact=TRUE)

# Example from Janson & Olsson (2001), Table 1
photo <- list()
photo[[1]] <- cbind(c( 71, 73, 86, 59, 71),  # weight ratings
                    c( 74, 80,101, 62, 83),
                    c( 76, 80, 93, 66, 77))
photo[[2]] <- cbind(c(166,160,187,161,172),  # height rating
                    c(171,170,174,163,182),
                    c(171,165,185,162,181))
iota(photo)
iota(photo, standardize=TRUE) # iota over standardized values
}
\keyword{univar}
