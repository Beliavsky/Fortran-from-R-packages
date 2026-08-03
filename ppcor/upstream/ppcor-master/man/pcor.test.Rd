\name{pcor.test}
\alias{pcor.test}
\title{
Partial correlation for two variables given a third variable.
}
\description{
The function \code{pcor.test} can calculate the pairwise partial correlations between two variables. In addition, it gives us the p value as well as statistic.
}
\usage{
pcor.test(x, y, z, method = c("pearson", "kendall", "spearman"))
}
\arguments{
  \item{x}{
a numeric vector.
}
  \item{y}{
a numeric vector.
}
  \item{z}{
a numeric vector.
}
  \item{method}{
a character string indicating which partial correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman" can be abbreviated.
}
}
\details{
Partial correlation is the correlation of two variables while controlling for a third variable. When the determinant of variance-covariance matrix is numerically zero, Moore-Penrose generalized matrix inverse is used. In this case, no \code{p-value} and \code{statistic} will be provided if the number of variables are greater than or equal to the sample size.
}
\value{
\item{estimate}{the partial correlation coefficient between two variables}
\item{p.value}{the p value of the test}
\item{statistic}{the value of the test statistic}
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
Missing values are not allowed
}
\seealso{
\code{\link[ppcor]{pcor}}, \code{\link[ppcor]{spcor}}, \code{\link[ppcor]{spcor.test}}
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

# partial correlation between "hl" and "disp" given "deg" and "BC"
pcor.test(y.data$hl,y.data$disp,y.data[,c("deg","BC")])
pcor.test(y.data[,1],y.data[,2],y.data[,c(3:4)])
pcor.test(y.data[,1],y.data[,2],y.data[,-c(1:2)])
}
\keyword{ htest }
