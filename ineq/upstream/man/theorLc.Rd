\name{theorLc}
\alias{theorLc}
\alias{Lc.dagum}
\alias{Lc.singh}
\alias{Lc.pareto}
\alias{Lc.lognorm}
\alias{Lc.exp}
\title{Theoretical Lorenz Curves}
\description{
 Theoretical Lorenz curves of income distributions}
\usage{
theorLc(type=c("Singh-Maddala","Dagum","lognorm","Pareto","exponential"), parameter=0)

Lc.dagum(p, parameter=c(2,2))
Lc.singh(p, parameter=c(2,2))
Lc.pareto(p, parameter=2)
Lc.lognorm(p, parameter=1)
Lc.exp(p)
}
\arguments{
 \item{type}{character string giving the income distribution.
   Must be one of the strings in the default argument (the first character
   is sufficient). Defaults to "Singh-Maddala".}
 \item{parameter}{vector containing parameter(s) of the distributions.}
 \item{p}{vector with elements from [0,1].}
}
\details{\code{Lc.dagum}, \code{Lc.singh}, \code{Lc.pareto}, \code{Lc.lognorm},
         \code{Lc.exp} are theoretical Lorenz curves of income distributions.
         They are functions of class \code{"theorLc"} with plot- and a lines-
         method, so that they can be added into an existing Lorenz curve plot.

         \code{theorLc} returns a function of class \code{"theorLc"}, that is a
         one of the above theoretical Lorenz curves with fixed parameters.

         \code{Lc.dagum} is the Lorenz curve of the Dagum distribution (2
         parameters), \code{Lc.singh} the one of the Singh-Maddala
         distribution (2 parameters), \code{Lc.pareto} the one of the Pareto
         distribution (1 parameter), \code{Lc.lognorm} the one of the Lognormal
         distribution (1 parameter) and \code{Lc.exp} the Lorenz curve of the
         exponential distribution (no parameter).
}
\value{A function of class \code{"theorLc"} or its value at \code{p}
respectively.}
\references{
C Dagum: Income Distribution Models, 1983, in: Johnson / Kotz (Eds): Encyclopedia of Statistical Sciences Vol.4, 27-34.

J B McDonald: Some generalized functions for the size distribution of income, 1984, Econometrica 52, 647-664.}

\seealso{\code{\link{Lc}}, \code{\link{plot.Lc}}, \code{\link{plot.theorLc}}}

\examples{
## Load and attach income (and metadata) set from Ilocos, Philippines
data(Ilocos)
attach(Ilocos)
## extract income for the province "Pangasinan"
income.p <- income[province=="Pangasinan"]

## plot empirical Lorenz curve and add theoretical Lorenz curve of
## a lognormal distribution with an estimate of the standard
## deviation parameter
Lc.p <- Lc(income.p)
plot(Lc.p)
lines(Lc.lognorm, parameter=sd(log(income.p)), col=4)



# vector of percentages
p <- (1:10)*0.1
# compute values of theoretic Lorenz curve of a Dagum-distribution
Lc.dagum(p, parameter=c(3.4,2.6))
# or
mydagum <- theorLc(type="Dagum", parameter=c(3.4,2.6))
mydagum(p)
}
\keyword{misc}

