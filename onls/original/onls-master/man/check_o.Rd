\name{check_o}
\alias{check_o}
\encoding{latin1}

\title{Check the fit for orthogonality of all points}

\description{
Checks for orthogonality of all points by calculating the angle between the slope of the tangent at \eqn{(x_{0i}, y_{0i})} and the slope of the Euclidean vector \eqn{\|\vec{D}_i\|} to \eqn{(x_i, y_i)}, which should be \eqn{90^{\circ}} if the Euclidean distance has been minimized. See 'Details'.
}

\usage{
check_o(object, plot = TRUE) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
\item{plot}{logical. If \code{TRUE}, the \eqn{\alpha}-values are plotted for a quick overview of orthogonality for all points.}
}

\value{
A dataframe containing \eqn{x_i}, \eqn{x_{0i}}, \eqn{y_i}, \eqn{y_{0i}}, \eqn{\alpha_i}, \eqn{\frac{df}{dx}} and a logical for \eqn{89.95^\circ < \alpha_i < 90.05^\circ}. If \code{plot = TRUE}, a plot of the \eqn{\alpha}-values in black if orthogonal, or red otherwise. 
}

\details{
This is a validation method for checking the orthogonality between all \eqn{(x_{0i}, y_{0i})} and \eqn{(x_i, y_i)}. The function calculates the angle between the slope \eqn{\mathrm{m}_i} of the tangent obtained from the first derivative at \eqn{(x_{0i}, y_{0i})} and the slope \eqn{\mathrm{n}_i} of the \code{\link{onls}}-minimized Euclidean distance between \eqn{(x_{0i}, y_{0i})} and \eqn{(x_i, y_i)}:
\deqn{\tan(\alpha_i) = \left|\frac{\mathrm{m}_i - \mathrm{n}_i}{1 + \mathrm{m}_i \cdot \mathrm{n}_i}\right|, \,\, \mathrm{m}_i = \frac{df(x, \beta)}{dx_{0i}}, \,\, \mathrm{n}_i = \frac{y_i - y_{0i}}{x_i - x_{0i}}}
\deqn{=> \alpha_i[^{\circ}] = \tan^{-1} \left( \left|\frac{\mathrm{m}_i - \mathrm{n}_i}{1 + \mathrm{m}_i \cdot \mathrm{n}_i}\right| \right) \cdot \frac{360}{2\pi}}
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Compare 'data range extended' orthogonal model with a
## 'data range restricted' model by setting "extend = c(0, 0)"
## => some x may not be orthogonal!
x <- 1:20
y <- 10 + 3*x^2
y <- sapply(y, function(a) rnorm(1, a, 0.1 * a))
DAT <- data.frame(x, y)

mod1 <- onls(y ~ a + b * x^2, data = DAT, start = list(a = 1, b = 1))
check_o(mod1)

mod2 <- onls(y ~ a + b * x^2, data = DAT, start = list(a = 1, b = 1),
             extend = c(0, 0)) 
check_o(mod2)
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

