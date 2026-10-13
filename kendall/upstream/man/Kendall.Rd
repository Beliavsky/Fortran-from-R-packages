\name{Kendall}
\alias{Kendall}
\title{ Kendall rank correlation }
\description{
Computes the Kendall rank correlation and its p-value 
on a two-sided test of H0: x and y are independent.
If there are no ties, the test is exact and in this case
it should agree with the base function cor(x,y,method="kendall")
and cor.test(x,y,method="kendall").
When there are ties, the normal approximation given in Kendall is
used as discussed below. 
In the case of ties,  both Kendall and cor produce the
same result but cor.test produces a p-value which is not as
accurate
}
\usage{
Kendall(x, y)
}
\arguments{ 
  \item{x}{ first variable, a vector }
  \item{y}{ second variable, a vector the same length as x }
}

\details{
In many applications x and y may be ranks or even ordered categorical 
variables. In our function x and y should be numeric vectors or factors.  
Any observations correspondings to NA in either x or y are removed.

Kendall's rank correlation measures the strength of monotonic association
between the vectors x and y. In the case of no ties in the x and y variables,
Kendall's rank correlation coefficient, tau, 
may be expressed as
\eqn{\tau = S/D}{tau=S/D}
where
\deqn{S=\sum_{i<j} (sign(x[j]-x[i])*sign(y[j]-y[i]))}
and \eqn{D=n(n-1)/2}{D=n(n-1)/2}.
S is called the score and D, the denominator, is the maximum possible value
of S. 
When there are ties, the formula for D is more complicated (Kendall, 1974, Ch. 3) and
this general forumla for ties in both reankings is implemented in our function.

The p-value of tau under the null hypothesis of no association is computed
by in the case of no ties using an exact algorithm given by Best and Gipps (1974).

When ties are present, a normal approximation with continuity correction
is used by taking S as
normally distributed with mean zero and variance var(S), where var(S)
is given by Kendall (1976, eqn 4.4, p.55).
Unless ties are very extensive and/or the data is very short, this 
approximation is adequate.
If extensive ties are present then the bootstrap provides an
expedient solution (Davis and Hinkley, 1997).
Alternatively an exact method based on exhaustive enumeration is also available 
(Valz and Thompson, 1994) but this is not implemented in this package.

An advantage of the Kendall rank correlation over the Spearman
rank correlation is that the score function S nearly normally 
distributed for small n and the distribution of S is easier to work with.

It may also be noted that usual Pearson correlation is fairly robust
and it usually agrees well in terms of statistical significance with
results obtained using Kendall's rank correlation.

An error is returned if length(x) is less than 3.
  
}

\value{
A list with class Kendall is returned with the following components:
  \item{tau}{ Kendall's tau statistic}
  \item{sl }{two-sided p-value}
  \item{S}{Kendall Score}
  \item{D }{denominator, tau=S/D}
  \item{varS }{variance of S}
}

\note{ 
Generic functions print.Kendall and summary.Kendall are provided.

If you want to use the output from Kendall,
save the result as in out<-Kendall(x,y) and then select from the list out the
value(s) needed. 
}

\references{
Best, D.J. and Gipps, P.G. (1974),
Algorithm AS 71: The Upper Tail Probabilities of Kendall's Tau
Applied Statistics, Vol. 23, No. 1. (1974), pp. 98-100.

Davison, A.C. and Hinkley, D.V. (1997) 
Bootstrap Methods and Their Application. Cambridge University Press.

Kendall, M.G. (1976). Rank Correlation Methods. 4th Ed. Griffin.

Hill, I.D. (1973),
Algorithm AS 66: The Normal Integral
Applied Statistics, Vol. 22, No. 3. (1973), pp. 424-427.

Valz, P. (1990).
Developments in Rank Correlation Procedures
with Applications to Trend Assessment in Water Resources Research,
Ph.D. Thesis, Department of Statistical and Actuarial Sciences, The
University of Western Ontario.

Valz, P.D. and Thompson, M.E. (1994),
Exact inference for Kendall's S and Spearman's rho. 
Journal of Computational and Graphical Statistics, 3, 459--472.
}
\author{ A.I. McLeod, aim@uwo.ca }

\seealso{ 
\code{\link{cor}}, 
\code{\link{print.Kendall}}, 
\code{\link{summary.Kendall}},
\code{\link{MannKendall}},
\code{\link{SeasonalMannKendall}}
   }

\examples{

#First Example
#From Kendall (1976, p.42-43, Example 3.4)
A<-c(2.5,2.5,2.5,2.5,5,6.5,6.5,10,10,10,10,10,14,14,14,16,17)
B<-c(1,1,1,1,2,1,1,2,1,1,1,1,1,1,2,2,2)
summary(Kendall(A,B))
#Kendall obtains S=34, D=sqrt(116*60), tau=0.41

#Second Example
#From Kendall (1976, p.55, Example 4.3)
x<-c(1.5,1.5,3,4,6,6,6,8,9.5,9.5,11,12)
y<-c(2.5,2.5,7,4.5,1,4.5,6,11.5,11.5,8.5,8.5,10)
summary(Kendall(x,y))
#Kendall obtains S=34 and Var(S)=203.30

}
\keyword{nonparametric}
