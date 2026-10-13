\name{perm}
\alias{permTS}
\alias{permTS.default}
\alias{permTS.formula}
%\alias{twosample.exact.network}
%\alias{twosample.exact.mc}
%\alias{twosample.exact.ce}
%\alias{twosample.pclt}
\alias{permKS}
\alias{permKS.default}
\alias{permKS.formula}
\alias{permTREND}
\alias{permTREND.default}
\alias{permTREND.formula}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exact or Asymptotic 2-sample, k-sample, and trend permutation tests}
\description{
These functions perform either: two-sample permutation tests (\code{permTS}),
k-sample permutation tests (\code{permKS}), or trend permutation tests (\code{permTREND}).
The test function can be transformed to a linear function of the scores times the covariate, where the covariate 
may be either a factor or character vector with two (\code{permTS}) or more (\code{permKS}) levels or a  
numeric vector (\code{permTREND}).   By using suitable scores one can create
for example, the permutation t-test (general scores), the Wilcoxon rank sum test (rank scores),
the logrank test (need to use other functions to create these scores). It performs either 
exact (network algorithm, complete enumeration, or Monte Carlo) asymptotic calculations (using permutational 
central limit theorem).   
}
\usage{
permTS(x, \dots)

\method{permTS}{default}(x, y, alternative = c("two.sided", "less", "greater"), 
    exact = NULL, method = NULL, methodRule = methodRuleTS1, 
    control=permControl(), \dots)

\method{permTS}{formula}(formula, data, subset, na.action, \dots)


permKS(x,\dots)


\method{permKS}{default}(x, g, exact = NULL, method = NULL, 
   methodRule = methodRuleKS1, control=permControl(), ...)

\method{permKS}{formula}(formula,data,subset, na.action,\dots)

permTREND(x,\dots)

\method{permTREND}{default}(x, y, alternative = c("two.sided", "less", "greater"), 
   exact = NULL, method = NULL, methodRule = methodRuleTREND1, control=permControl(),...)

\method{permTREND}{formula}(formula,data,subset,na.action,\dots)
}
\arguments{
\item{x}{numeric vector of respose scores for the first group}
\item{y}{numeric vector of either response scores for the second group (for permTS) 
or trend scores for each observation (for permTREND)}
\item{g}{a factor or character vector denoting group membership}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided"
(default), "greater","less" (see details)}
\item{exact}{a logical value, TRUE denotes exact test, ignored if method is not NULL}
\item{method}{a character value, one of 'pclt','exact.network','exact.ce','exact.mc'. If NULL method chosen by methodRule}
\item{methodRule}{a function used to choose the method (see details)}
\item{control}{a list with arguments that control the algortihms, see \code{\link{permControl}}}
\item{formula}{a formula of the form lhs~rhs where lhs is a numeric variable giving the response 
scores and rhs a factor with two levels giving the corresponding groups. }
\item{data}{an optional matrix or data frame containing the variables in the formula}
\item{subset}{an optional vector specifying a subset of observations to be used.}
\item{na.action}{a function which indicates what should happen when the data contain NAs. Defaults to 
getOption("na.action"). }
\item{\dots}{ further arguments to be passed to or from methods. }
}
\details{

There are 4 different methods for deciding how to determine the p-value by defining which test statistics are extreme.
For \code{alternative} there are 3 choices, "two.sided", "less" or "greater", but within alternative="two.sided" there are 2 methods defined by the 
\code{tsmethod} given within \code{control}, see \code{\link{permControl}}.   If Ti is a vector of test statistics, and T0 is the observed test statistic, then 
alternative="less" gives p.lte=Pr[Ti<=T0], alternative="greater" gives p.gte=Pr[Ti>=T0], 
alternative="two.sided" with tsmethod="central" (default) gives p.twosided=max(1, 2*min(p.lte,p.gte)), and alternative="two.sided"
with tsmethod="abs" gives p.twosidedAbs=Pr[abs(Ti - mean(Ti) ) >=abs(T0-mean(Ti))]. For \code{permTS} the test statistic 
is equivalent to the mean of one group minus the mean of the other group. For \code{permTREND} the test 
statistic is equivalent to the correlation between the response (x) and the trend scores (y).
For \code{permKS} only a twosided pvalue based on Pr[Ti>=T0] is allowed, where the test statistic, Ti, is the 
weighted sum of the square of the mean within group, where the weights are the sample size for each group.  This will
give for example, the usual Kruskal-Wallis test when the ranks are used on the responses.

Many standard statistical tests may be put into the form of the permutation test (see Graubard and Korn, 1987). 
There is a choice of four different methods to calculate the p-values (the last two are only available for 
\code{permTS}): 
\enumerate{
\item pclt: using permutational central limit theorem (see e.g., Sen, 1985).
\item exact.mc:exact using Monte Carlo.
\item exact.network: exact method using a network algorithm (see e.g., Agresti, Mehta, and Patel, 1990). Currently the network 
method does 
not implement many of the time saving suggestions such as clubbing.
\item exact.ce: exact using complete enumeration. This is good for very small sample sizes and when doing simulations, since the cm need only 
be calculated once for the simulation.
}
The \code{exact.network} and \code{exact.ce} may give errors related to running out of memory when the sample size is not small and will depend on the system you are using (e.g., about 15 in each group for \code{exact.network} or 
14 in each group for \code{exact.ce}).

These associated functions for the above methods (e.g.,  \code{twosample.pclt}, \code{twosample.exact.network}, etc),  
are internal and are not to be called directly. 

The \code{methodRule} is a function which takes the first two objects of the default implementation, and returns the 
method. This function can be used to appropriately choose the method based on the size of the data.
For explanation of the default method rules see \code{\link{methodRuleTS1}}, \code{\link{methodRuleKS1}}, or 
\code{\link{methodRuleTREND1}}. 

For more details see Fay and Shaw (2010, Section 5).
}
\value{
An object of class \code{htest} or for 'exact.mc' of class \code{mchtest}, 
a list with the following elements:
  \item{p.value}{p value associated with alternative}
  \item{alternative}{description of alternative hypothesis}
  \item{p.values}{a vector giving lower, upper, and two-sided p-values as well as 
p.equal which is the proportion equal to the observed test statistic}
  \item{method}{a character vector describing the test}
\item{estimate}{an estimate of the test statistic}
\item{statistic}{statistic used for asymptotics, either Z statistics or chi square statistic, output if method="pclt"}
\item{parameter}{degrees of freedom for chi square statistic, output if 'statistic' is the chi square statistic}
\item{data.name}{character vector describing the response and group variables}
  \item{p.conf.int}{a confidence interval on the p-value if method='exact.mc' 
    (see \code{\link{calcPvalsMC}})}
  \item{nmc}{number of Monte Carlo replications if method='exact.mc', NULL otherwise}
}
\references{ 

Agresti, A, Mehta, CR, Patel, NR (1990). JASA 85: 453-458.

    Fay, MP and Shaw, PA (2010). Exact and Asymptotic Weighted Logrank Tests for Interval Censored Data: The interval R package.
    Journal of Statistical Software. \doi{10.18637/jss.v036.i02}. 36 (2):1-34.

Graubard, BI, and Korn, EL (1987). Biometrics 43: 471-476.

Sen, PK (1985) `Permutational central limit theorems' in Encyclopedia of Statistics, Vol 6.
}
\author{Michael Fay}
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\examples{
## Example from StatExact manual
dBP<-c(94,108,110,90,80,94,85,90,90,90,108,94,78,105,88)
treatment<-c(rep("treated",4),rep("control",11))
permTS(dBP~treatment,alternative="less",method="pclt")
result<-permTS(dBP[treatment=="treated"],dBP[treatment=="control"],alternative="greater")
result
result$p.values
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric}% __ONLY ONE__ keyword per line
