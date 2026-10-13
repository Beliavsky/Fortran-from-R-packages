\name{permControl}
\alias{permControl}
\title{Auxiliary for controling permutation tests}
\description{
A function to create a list of arguments for \code{\link{permTS}}, \code{\link{permKS}} or 
\code{\link{permTREND}}.
}
\usage{
permControl(cm=NULL,nmc=10^3-1,seed=1234321,digits=12,
    p.conf.level=.99,setSEED=TRUE,tsmethod="central")
}
\arguments{
\item{cm}{a choose(n,m) by n matrix, used if method='exact.ce', ignored otherwise}
\item{nmc}{number of Monte Carlo replications, used if method='exact.mc', ignored otherwise}
\item{seed}{value used in \code{set.seed} if method='exact.mc', ignored otherwise} 
\item{setSEED}{logical, set to FALSE when performing simulations that use method='exact.mc'}
\item{p.conf.level}{confidence level for p value estimate, used if method='exact.mc', ignored otherwise }
\item{digits}{number of digits to use in \code{\link{signif}} for precision of test statistics}
\item{tsmethod}{method for calculating two-sided p-values, character, either 'central' or 'abs' (see details)}
}
\details{

When \code{cm}=NULL the resulting matrix is created by \code{\link{chooseMatrix}}, it may be optionally provided here 
only so that \code{chooseMatrix} does not need to be repeatedly called in simulations.  Also when doing simulations with
method='exact.mc', use setSEED=FALSE so that the seed is not reset to the same value each time you call the permutation test
function.

See \code{\link{calcPvalsMC}} for description of how p.conf.level is used. 

The two-sided method is given by \code{tsmethod}. The default 'central' two-sided method is 
just p=min(1, 2*min(pless,pgreater)), where pless and pgreater are the one-sided 
p-values. The name 'central' follows the convention of the \code{exact2x2} and 
\code{exactci} packages; so that for example, a two-sample permutation test on a binary
response with tsmethod='central' will match the Central Fisher's Exact test (see 
Fay, 2010). The option tsmethod='abs' defines another method for defining the two-sided p-value. We define it for complete enumeration, 
but the algorithms may differ. Let Tj be the vector of length N of all possible values of the test statistic under each of N possible permutations.
The p-value with tsmethod='abs' is defined as 1/N times the number of times \code{abs(Tj-mean(Tj)) >= abs(T0-mean(Tj))}, where T0 is the observed value of the 
test statistic. This option matches the default two-sided method for the \code{coin} package. 
}
\value{
An list with the arguments as components.
}
\references{ 
Fay, M.P. (2010). Confidence intervals that match Fisher's exact or 
    Blaker's exact tests. Biostatistics 11(2):373-373.
}
%\author{Michael Fay}
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{misc}
