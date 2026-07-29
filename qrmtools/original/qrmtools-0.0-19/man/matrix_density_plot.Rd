\name{matrix_density_plota}
\alias{matrix_density_plot}
\title{Density Plot of the Values from a Lower Triangular Matrix}
\description{
  Density plot of all values in the lower triangular part of a matrix.
}
\usage{
matrix_density_plot(x, xlab = "Entries in the lower triangular matrix",
                    main = "", text = NULL, side = 4, line = 1, adj = 0, ...)
}
\arguments{
  \item{x}{\code{\link{matrix}}-like object.}
  \item{xlab}{x-axis label.}
  \item{main}{title.}
  \item{text}{see \code{\link{mtext}()}. The \code{text = ""}, it is omitted.}
  \item{side}{see \code{\link{mtext}()}.}
  \item{line}{see \code{\link{mtext}()}.}
  \item{adj}{see \code{\link{mtext}()}.}
  \item{\dots}{additional arguments passed to the underlying \code{\link{plot}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  \code{matrix_density_plot()} is typically used for symmetric matrices
  (like correlation matrices, matrices of pairwise Kendall's tau or tail
  dependence parameters) to check the distribution of their off-diagonal
  entries.
}
\author{Marius Hofert}
\examples{
## Generate a random correlation matrix
d <- 50
L <- diag(1:d)
set.seed(271)
L[lower.tri(L)] <- runif(choose(d,2))
Sigma <- L %*% t(L) # ... real, symmetric, positive definite (covariance) matrix Sigma
P <- cor(Sigma)
## Density of its lower triangular entries
matrix_density_plot(P)
}
\keyword{hplot}