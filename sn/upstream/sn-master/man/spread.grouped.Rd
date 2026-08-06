\name{spread.grouped}
\alias{spread.grouped}
\concept{grouped data}

\title{Spreading grouped data over intervals}
\description{
Assuming that \code{counts} represents the frequencies of observations
falling into intervals identified by \code{breaks}, the function returns
a vector of values obtained by uniformly spreading each group of data 
over the pertaining interval.   
}

\usage{spread.grouped(breaks, counts, shift = "centre")}

\arguments{
  \item{breaks}{A numeric vector of strictly increasing finite values 
     which identify a set of contiguous intervals on the real line.}

  \item{counts}{A vector of non-negative integers representing the number 
     of observations falling in the intervals specified by \code{breaks};
     it is then required that \code{length(counts)+1=length(breaks)}. }
     
  \item{shift}{a character string which regulates the positioning of the 
    constructed points within a given interval, with possible values
    \code{"left"}, \code{"center"} (default choice) and \code{"right"},
    possibly abbreviated.}
}

\value{A numeric vector of length \code{sum(counts)} of values 
  within \code{range(breaks)}.}

\author{Adelchi Azzalini}

\seealso{\code{\link{fitdistr.grouped}}, \code{\link{fitdistr.grouped-class}} }

\examples{
breaks <- c(10, 12, 15, 20)
counts <- c(3, 2, 4)
spread.grouped(breaks, counts)
spread.grouped(breaks, counts, "l")
}
\keyword{manip}

