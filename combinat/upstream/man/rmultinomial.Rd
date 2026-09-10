\name{rmultinomial}
\alias{rmultinomial}
\alias{rmultz2}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generate random samples from multinomial distributions }
\description{
 rmultinomial: Generate random samples from multinomial distributions, where both n
 and p may vary among distributions
 
 rmultz2: fixed p case
}
\usage{
rmultinomial(n, p, rows=max(c(length(n), nrow(p))))
rmultz2(n, p, draws=length(n))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{ vector of sizes }
  \item{p}{ vector or probs }
  \item{rows}{ numeric giving desired number rows to be output }
  \item{draws}{ number samples required }
}
%\details{
%}
\value{
a matrix of \code{rows} rows delivering specified samples
}
%\references{ ~put references to the literature/web site here ~ }
\author{ John Wallace, 17 Feb 1997 S-news , mods by Chasalow }
\examples{
 n <- c(100,20,10)
 p <- matrix(c(.3,.1,.5,.1,.1,.2,.6,.8,.3),3)
 rmultinomial(n,p)
}

\keyword{ models }

