\name{TTT}
\alias{TTT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TTT function
}
\description{
There are several behaviors that the failure rate function of a random variable T can take. In this context, the graph of total test time (TTT curve) proposed by Aarset (1987) may be used for obtaining empirical behavior of the function failure rate.
}
\usage{
TTT(x, lwd = 2, lty = 2, col = "black", grid = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Data vector;
}
  \item{lwd}{
  Thickness of the TTT curve. The argument \code{lwd} must be a nonnegative real number;
}
  \item{lty}{
  The argument \code{lty} modifies the style of the diagonal line chart TTT. Possible values are: 0 [blank], 1 [solid (default)], 2 [dashed], three [dotted], 4 [dotdash], 5 [longdash], 6 [twodash];
}
  \item{col}{
Color used in the TTT curve;
}
  \item{grid}{
If \code{grid = FALSE} graphic appears without the grid;
}
  \item{\dots}{
    Other arguments passed by the user and available for the function \code{plot}. More details in \code{\link{par}}.
}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Aarset, M. V. (1987). How to identify bathtub hazard rate. IEEE Transactions Reliability, 36, 106-108.
}
\author{
Pedro Rafael Diniz Marinho (pedro.rafael.marinho@gmail.com);

Marcelo Bourguignon (m.p.bourguignon@gmail.com).
}
\note{

The graphic TTT may have various forms. Aarset (1987) showed that if the curve approaches a straight diagonal function constant failure rate is adequate. When the curve is convex or concave the failure rate function is monotonically increasing or decescente respectively is adequate. If the failure rate function is convex and concave, the failure rate function in format \eqn{U} is adequate, otherwise the failure rate function unimodal is more appropriate.

The TTT curve is constructed by values \eqn{r/n} and \eqn{G(r/n)}, wherein
\deqn{
G(r/n) =  \frac{[\sum_{i=1}^{r} T_{i:n} + (n-r)T_{r:n}]}{\sum_{i=1}^{n}T_{i:n}}, r = 1, \ldots, n, T_{1:n} = 1, \ldots, n.
}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(carbone)
TTT(carbone, col = "red", lwd = 2.5, grid = TRUE, lty = 2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
