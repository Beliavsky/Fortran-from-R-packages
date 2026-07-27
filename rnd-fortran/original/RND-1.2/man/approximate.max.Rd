\name{approximate.max}
\alias{approximate.max}
\title{Max Function Approximation}
\description{\code{approximate.max} gives a smooth approximation to the max function.}
\usage{
approximate.max(x, y, k = 5)
}
\arguments{
  \item{x}{the first argument for the max function}
  \item{y}{the second argument fot the max function}
  \item{k}{a tuning parameter.  The larger this value, the closer the function output to a true max function.}
}
\details{
 \code{approximate.max} approximates the max of x, and y as follows: 
 \deqn{
 g(x,y) = \frac{1}{1 + \exp(-k(x-y))}, \ \
 \max(x,y) \approx x g(x,y) + y(1 - g(x,y))
 }
}
\value{approximate maximum of x and y}
\references{
  Melick, W. R. and Thomas, C.P. (1997)
  Recovering an asset's implied pdf from option proces:
  An application to crude oil during the gulf crisis.
  \emph{Journal of Financial and Quantitative Analysis}, 32(1), 91-115
}
\author{Kam Hamidieh}
\examples{
#
# To see how the max function compares with approximate.max,
# run the following code.
#

i = seq(from = 0, to = 10, by = 0.25)
y = i - 5
max.values = pmax(0,y)
approximate.max.values = approximate.max(0,y,k=5)
matplot(i, cbind(max.values, approximate.max.values), lty = 1, type = "l", 
 col=c("black","red"), main = "Max in Black, Approximate Max in Red")
}
\keyword{ maximum }
\keyword{ numerical tricks }

