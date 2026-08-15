\name{skewunit}
\alias{dskewunit}
\alias{pskewunit}
\alias{rskewunit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Family of Skew Distributions with Bounded Support
}
\description{
Density, distribution function and random generation for a family of skew distributions with bounded support.
}
\usage{
dskewunit(x, lambda = 0, delta = 1, delta2 = 1, family1 = "asin", family2 = "asin", 
          log = FALSE)
pskewunit(x, lambda = 0, delta = 1, delta2 = 1, family1 = "asin", family2 = "asin", 
          log.p = FALSE)
rskewunit(n, lambda = 0, delta = 1, delta2 = 1, family1 = "asin", family2 = "asin")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles.
}
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
  \item{lambda}{
skewness parameter such as \eqn{-1\leq \lambda \leq 1}.
}
  \item{delta, delta2}{
shape parameters.
}
  \item{family1}{
first family of distributions related to \eqn{f} (asin by default). See details Section.
}
  \item{family2}{
second family of distributions related to \eqn{G} (asin by default, NULL returns the model f). See details Section.
}
  \item{log, log.p}{
logical; if TRUE, probabilities p are given as log(p).
}
}
\details{
The Family of Skew Distributions with Bounded Support is defined by its density function given by
\deqn{
f(x)=2 G(\lambda(x-0.5)+0.5), \quad x \in (0,1), \lambda \in (-1,1),
}
where \eqn{f} is symmetric around 0.5, i.e., \eqn{f(x-0.5)=f(x+0.5)}.
The available options for family1 and family2 are asin, Uquad, triang, JSB and sbeta.
When family2 is NULL, f corresponds to the model specified in family1.
%
}
\value{
dskewunit gives the density, pskewunit gives the distribution function and rskewunit generates random deviates.
The length of the result is determined by n for rnorm, and is the maximum of the lengths of the numerical arguments for the other functions.
The numerical arguments other than n are recycled to the length of the result. Only the first elements of the logical arguments are used.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Diego Gallardo, Emilio Gomez-Deniz, Osvaldo Venegas and Hector W. Gomez
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
dskewunit(c(0.2,0.8), lambda = 0.5, family1 = "asin", family2 = "asin")
pskewunit(c(0.2,0.8), lambda = 0.5, family1 = "asin", family2 = "asin")
rskewunit(100, lambda = -0.4, delta = 1, family1 = "triang", family2 = "JSB")
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
