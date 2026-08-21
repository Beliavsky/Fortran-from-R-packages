\name{descriptive}
\alias{descriptive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
descriptive - Calculation of descriptive statistics
}
\description{
The function \code{descriptive} calculates the main descriptive statistics of a vector of data.
}
\usage{
descriptive(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data vector.
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
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988). The New S Language. Wadsworth & Brooks/Cole.

}
\author{
Pedro Rafael Diniz Marinho \email{pedro.rafael.marinho@gmail.com}

Marcelo Bourguignon \email{m.p.bourguignon@gmail.com}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(carbone)
descriptive(carbone)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
