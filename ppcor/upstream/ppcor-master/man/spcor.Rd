\name{spcor}
\alias{spcor}
\title{
Semi-partial (part) correlation
}
\description{
The function \code{spcor} can calculate the pairwise semi-partial (part) correlations for each pair of variables given others. In addition, it gives us the p value as well as statistic for each pair of variables.
}
\usage{
spcor(x, method = c("pearson", "kendall", "spearman"))
}
\arguments{
  \item{x}{
a matrix or data fram. 
}
  \item{method}{
a character string indicating which semi-partial (part) correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman" can be abbreviated.
}
}
\details{
Semi-partial correlation is the correlation of two variables with variation from a third or more other variables removed only from the second variable. When the determinant of variance-covariance matrix is numerically zero, Moore-Penrose generalized matrix inverse is used. In this case, no \code{p-value} and \code{statistic} will be provided if the number of variables are greater than or equal to the sample size.
}
\value{
\item{estimate}{a matrix of the semi-partial (part) correlation coefficient between two variables}
\item{p.value}{a matrix of the p value of the test}
\item{statistic}{a matrix of the value of the test statistic}
\item{n}{the number of samples}
\item{gn}{the number of given variables}
\item{method}{the correlation method used}
}
\references{
Kim, S. (2015) ppcor: An R Package for a Fast Calculation to Semi-partial Correlation Coefficients. Communications for Statistical Applications and Methods,  22(6), 665-674. 
}
\author{
Seongho Kim <\email{biostatistician.kim@gmail.com}>
}
\note{
Missing values are not allowed.
}
\seealso{
\code{\link[ppcor]{spcor.test}}, \code{\link[ppcor]{pcor}}, \code{\link[ppcor]{pcor.test}}
}
\examples{
# data
y.data <- data.frame(
				hl=c(7,15,19,15,21,22,57,15,20,18),
				disp=c(0.000,0.964,0.000,0.000,0.921,0.000,0.000,1.006,0.000,1.011),
				deg=c(9,2,3,4,1,3,1,3,6,1),
				BC=c(1.78e-02,1.05e-06,1.37e-05,7.18e-03,0.00e+00,0.00e+00,0.00e+00
              ,4.48e-03,2.10e-06,0.00e+00)
			)

# semi-partial (part) correlation
spcor(y.data) 
}
\keyword{htest}
