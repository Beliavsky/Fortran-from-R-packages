\name{perm-internal}
\alias{ksample.exact.mc}
\alias{ksample.pclt}
\alias{trend.exact.mc}
\alias{trend.pclt}
\alias{twosample.exact.ce}
\alias{twosample.exact.mc}
\alias{twosample.exact.network}
\alias{twosample.pclt}
\alias{getcnt}

\title{Internal Functions for perm package: Not to be called by User}
\description{
These functions are where the algorithms are done. There is much room for improvement in the speed of the 
exact functions. 
}
\usage{
ksample.exact.mc(scores, group, nmc = 10^4 - 1, seed = 1234321, 
  digits = 12, p.conf.level = 0.99, setSEED = TRUE)
ksample.pclt(scores, group)
trend.exact.mc(scores, group, alternative = "two.sided", nmc = 10^3 - 1, 
  seed = 1234321, digits = 12, p.conf.level = 0.99, setSEED = TRUE)
trend.pclt(scores, group)
twosample.exact.ce(scores, group, cm = NULL, digits = 12)
twosample.exact.mc(scores, group, alternative = "two.sided", nmc = 10^4 - 1, 
  seed = 1234321, digits = 12, p.conf.level = 0.99, setSEED = TRUE)
twosample.pclt(scores, group)

twosample.exact.network(scores, group, digits = 12)
getcnt(nodehk, cnt.edge, edgesize)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{vector of response scores}
  \item{group}{covariate vector}
  \item{alternative}{one of 'less', 'greater', 'two.sided' or 'two.sidedAbs'}
  \item{nmc}{number of Monte Carlo replications}
  \item{seed}{random number seed}
  \item{digits}{digits for rounding of test statistic, equal to that many digits are called tied}
  \item{p.conf.level}{confidence level for p-value, used with mc methods}
  \item{setSEED}{logical, set to FALSE when performing simulations on mc methods}
  \item{cm}{for speed you can input the matrix created from chooseMatrix (see \code{\link{chooseMatrix}}). If NULL it is created.} 
  \item{nodehk}{ nodes for which indeces of arcs are needed}
  \item{cnt.edge}{vector of first index for each node}
  \item{edgesize}{vector of number of arcs for each node}
}
\details{
Network algorithm is very basic, only works for two group tests. The function \code{getcnt} 
(called by twosample.exact.network)
gets a vector of indeces representing arcs for set of nodes. 
}
\value{
The function \code{getcnt} returns 
a vector of indeces representing arcs for set of nodes
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
%\examples{ }
\keyword{internal}
