\name{NA_plot}
\alias{NA_plot}
\title{Graphical Tool for Visualizing NAs in a Data Set}
\description{
  Plot NAs in a data set.
}
\usage{
NA_plot(x, col = c("black", "white"), xlab = "Time", ylab = "Component",
        text = "Black: NA; White: Available data",
        side = 4, line = 1, adj = 0, ...)
}
\arguments{
  \item{x}{matrix (ideally an \code{xts} object).}
  \item{col}{bivariate vector containing the colors for missing and
    available data, respectively.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{text}{see \code{\link{mtext}()}. The \code{text = ""}, it is omitted.}
  \item{side}{see \code{\link{mtext}()}.}
  \item{line}{see \code{\link{mtext}()}.}
  \item{adj}{see \code{\link{mtext}()}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{image}()}.}
}
\value{
  \code{\link{invisible}()}.
}
\details{
  Indicate \code{\link{NA}}s in a data set.
}
\author{Marius Hofert}
\examples{
## Generate data
n <- 1000 # sample size
d <- 100 # dimension
set.seed(271) # set seed
x <- matrix(runif(n*d), ncol = d) # generate data

## Assign missing data
k <- ceiling(d/4) # fraction of columns with some NAs
j <- sample(1:d, size = k) # columns j with NAs
i <- sample(1:n, size = k) # 1:i will be NA in each column j
X <- x
for(k. in seq_len(k)) X[1:i[k.], j[k.]] <- NA # put in NAs

## Plot NAs
NA_plot(X) # indicate NAs
}
\keyword{hplot}