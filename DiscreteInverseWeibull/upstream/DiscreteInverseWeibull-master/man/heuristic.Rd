\name{heuristic}
\alias{heuristic}
\title{
Heuristic method of estimation
}
\description{
Heuristic method for the estimation of parameters of the discrete inverse Weibull
}
\usage{
heuristic(x, beta1=1, z = 0.1, r = 0.1, Leps = 0.01)
}
\arguments{
  \item{x}{
a vector of sample values
}
  \item{beta1}{
launch value of the \eqn{\beta} parameter
}
  \item{z}{
initial value of width
}
  \item{r}{
initial value of rate
}
  \item{Leps}{
tolerance error for the likelihood function
}
}
\value{
a list containig the two estimates of \eqn{q} and \eqn{\beta}
}
\details{
For a detailed description of the method, have a look at the reference
}
\seealso{
\code{\link{estdiweibull}}
}
\examples{
n<-50
q<-0.25
beta<-1.5
x<-rdiweibull(n, q, beta)
# estimates using the heuristic algorithm
par0<-heuristic(x)
par0
# change the default values of some working parameters...
par1<-heuristic(x, beta1=2)
par1
par2<-heuristic(x, z=0.5)
par2
par3<-heuristic(x, r=0.2)
par3
par4<-heuristic(x, Leps=0.1)
par4
# ...there should be just light differences among the estimates...
# ... and among the corresponding values of the loglikelihood functions
loglikediw(x, par0[1], par0[2])
loglikediw(x, par1[1], par1[2])
loglikediw(x, par2[1], par2[2])
loglikediw(x, par3[1], par3[2])
loglikediw(x, par4[1], par4[2])
}
\references{
Jazi M.A., Lai C.-D., Alamatsaz M.H. (2010) A discrete inverse Weibull distribution and estimation of its parameters, Statistical Methodology, 7: 121-132

Drapella A. (1993) Complementary Weibull distribution: unknown or just forgotten, Quality Reliability Engineering International 9: 383-385
}
\keyword{distribution, htest}
