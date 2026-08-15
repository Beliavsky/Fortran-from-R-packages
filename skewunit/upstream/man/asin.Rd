\name{asin}
\alias{dasin}
\alias{pasin}
\alias{rasin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The ArcSin distribution.
}
\description{
Density, distribution function and random generation for the ArcSin distribution.
}
\usage{
dasin(x, log=FALSE)
pasin(q, lower.tail=TRUE, log.p=FALSE)
rasin(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{vector of quantiles.
}
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
\item{log, log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
\item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X>x]}.
}
}
\details{
The ArcSin distribution has density
\deqn{
f(x)=\frac{1}{\pi \sqrt{x(1-x)}}, \quad x \in (0,1),
}
and cumulative distribution function
\deqn{
F(x)=\frac{2}{\pi}Arcsin(\sqrt{x}), \quad x \in (0,1).
}
}
\value{
dasin gives the density, pasin gives the distribution function, and rasin generates random deviates.
The length of the result is determined by n for rasin, and is the maximum of the lengths of the numerical arguments for the other functions.
The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
}
\references{
Arnold, B.C. and Groeneveld, R.A. (1980). Some Properties of the Arcsine Distribution. Journal
of the Ammerican Statistical Association, 75, 173-175.
}
\author{
Diego Gallardo
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
dasin(0.5)
pasin(0.5)
rasin(5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
