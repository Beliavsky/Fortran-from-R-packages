\name{sortDataFrame}
\alias{sortDataFrame}

\title{
Sorting a data frame by a key
}
\description{
A simple function to sort a data frame based on a certain keyword. This function was posted by r-fanatic at a dzone forum (the webpage is not available anymore).
}
\usage{
sortDataFrame(x, key, ...)
}

\arguments{
  \item{x}{
  The data frame to be sorted.
}
  \item{key}{
  The key by which the data frame will be sorted.
}
  \item{\dots}{
  Further arguments to be passed to the order function.
}
}
\value{
The sorted data frame.
}
\references{
The original webpage where r-fanatic posted the code is not available as of 3rd April 2017.
}
\author{
r-fanatic
}

\examples{
# Create a simple data frame and order using the "B" key
ppp <- data.frame(A=1:10,B=10:1)
ppp
sortDataFrame(ppp,"B")
ppp
}

\keyword{ misc }
\keyword{ methods }
