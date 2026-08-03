\name{pcor}
\alias{pcor}
\title{
Partial correlation
}
\description{
The function \code{pcor} can calculate the pairwise partial correlations for each pair of variables given others. In addition, it gives us the p value as well as statistic for each pair of variables.
}
\usage{
pcor(x, method = c("pearson", "kendall", "spearman"))
}
\arguments{
  \item{x}{
a matrix or data fram. 
}
  \item{method}{
a character string indicating which partial correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman" can be abbreviated.
}
}
\details{
Partial correlation is the correlation of two variables while controlling for a third or more other variables. When the determinant of variance-covariance matrix is numerically zero, Moore-Penrose generalized matrix inverse is used. In this case, no \code{p-value} and \code{statistic} will be provided if the number of variables are greater than or equal to the sample size.
}
\value{
\item{estimate}{a matrix of the partial correlation coefficient between two variables}
\item{p.value}{a matrix of the p value of the test}
\item{statistic}{a matrix of the value of the test statistic}
\item{n}{the number of samples}
\item{gn}{the number of given variables}
\item{method}{the correlation method used}
}
\references{
Kim, S. (2015) ppcor: An R Package for a Fast Calculation to Semi-partial Correlation Coefficients. Communications for Statistical Applications and Methods, 22(6), 665-674.
}
\author{
Seongho Kim <\email{biostatistician.kim@gmail.com}>
}
\note{
Missing values are not allowed.
}
\seealso{
\code{\link[ppcor]{pcor.test}}, \code{\link[ppcor]{spcor}}, \code{\link[ppcor]{spcor.test}}
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
# partial correlation
pcor(y.data) 
}
\keyword{htest}
