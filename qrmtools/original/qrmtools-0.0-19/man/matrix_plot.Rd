\name{matrix_plot}
\alias{matrix_plot}
\title{Graphical Tool for Visualizing Matrices}
\description{
  Plot of a matrix.
}
\usage{
matrix_plot(x, ran = range(x, na.rm = TRUE), ylim = rev(c(0.5, nrow(x) + 0.5)),
            xlab = "Column", ylab = "Row",
            scales = list(alternating = c(1,1), tck = c(1,0),
                          x = list(at = pretty(1:ncol(x)), rot = 90),
                          y = list(at = pretty(1:nrow(x)))),
            at = NULL, colorkey = NULL, col = c("royalblue3", "white", "maroon3"),
            col.regions = NULL, ...)
}
\arguments{
  \item{x}{\code{\link{matrix}}-like object.}
  \item{ran}{range (can be used to enforce (-1,1), for example).}
  \item{ylim}{y-axis limits in reverse order (for the rows to appear 'top down').}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{scales}{see \code{\link[lattice]{levelplot}()}; if \code{\link{NULL}},
    labels and ticks are omitted.}
  \item{at}{see \code{\link[lattice]{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on the given values in \code{x}.}
  \item{colorkey}{see \code{\link[lattice]{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on \code{at}.}
  \item{col}{\code{\link{vector}} of length two (if all values of
    \code{x} are non-positive or all are non-negative; note that also a
    vector of length three is allowed in this case) or three (if
    \code{x} contains negative and positive values) providing the color
    key's default colors.}
  \item{col.regions}{see \code{\link[lattice]{levelplot}()}. If \code{\link{NULL}}, a
    useful default is computed based on \code{at}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link[lattice]{levelplot}()}.}
}
\value{
  The plot, a Trellis object.
}
\details{
  Plot of a matrix.
}
\author{Marius Hofert}
\examples{
## Generate a random correlation matrix
d <- 50
L <- diag(1:d)
set.seed(271)
L[lower.tri(L)] <- runif(choose(d,2)) # random Cholesky factor
Sigma <- L %*% t(L)
P <- cor(Sigma)

## Default
matrix_plot(P)
matrix_plot(P, ran = c(-1, 1)) # within (-1, 1)
matrix_plot(abs(P)) # if nonnegative
L. <- L
diag(L.) <- NA
matrix_plot(L.) # Cholesky factor without diagonal

## Default if nonpositive
matrix_plot(-abs(P))

## Changing colors
matrix_plot(P, ran = c(-1, 1),
            col.regions = grey(c(seq(0, 1, length.out = 100),
                                 seq(1, 0, length.out = 100))))

## An example with overlaid lines
library(lattice)
my_panel <- function(...) {
    panel.levelplot(...)
    panel.abline(h = c(10, 20), v = c(10, 20), lty = 2)
}
matrix_plot(P, panel = my_panel)
}
\keyword{hplot}