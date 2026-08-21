\name{CompoundIndicators}
\Rdversion{2.1}
\alias{CompoundIndicators}
\alias{arpt.cgb2}
\alias{arpr.cgb2}
\alias{rmpg.cgb2}
\alias{qsr.cgb2}
\alias{main.cgb2}
\title{Indicators of Poverty and Social Exclusion under the Compound Distribution based on the GB2
}
\description{
Functions to calculate four primary social welfare indicators under the compound GB2 distribution, i.e. the at-risk-of-poverty threshold,
 the at-risk-of-poverty rate, the relative median at-risk-of-poverty gap, and the income quintile share ratio.
}
\usage{
arpt.cgb2(prop, shape1, scale, shape2, shape3, pl0, pl, decomp="r")
arpr.cgb2(prop, shape1, shape2, shape3, pl0, pl, decomp="r")
rmpg.cgb2(arpr, shape1, shape2, shape3, pl0, pl, decomp="r")
qsr.cgb2(shape1, shape2, shape3, pl0, pl, decomp="r")
main.cgb2(prop, shape1, scale, shape2, shape3, pl0, pl, decomp="r")
}
\arguments{
  \item{prop}{numeric; proportion (in general is set to 0.6).}
  \item{arpr}{numeric; the value of the at-risk-of-poverty rate.}
  \item{shape1,scale,shape2,shape3}{numeric; positive parameters of the GB2 distribution.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{pl}{numeric; a vector of mixture probabilities. Sums to one. If \eqn{pl=pl0} we obtain the GB2 distribution.}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution.
   By default, \code{decomp} = "r" - right tail decomposition.}
}
\details{
The four indicators are described in details in the case of the GB2. The difference here is that we need to give an initial vector of proportions, fitted proportions and define for which decomposition (left or right) the indicators should be calculated.
}
\value{
\code{arpt.cgb2} gives the ARPT, \code{arpr.cgb2} the ARPR, \code{rmpg.cgb2} the RMPG, \code{qsr.cgb2} gives the QSR and \code{main.cgb2} calculates
the median, the mean, the ARPR, the RMPG and the QSR under the compound GB2.
}
\author{
Monique Graf
}
\references{
 Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
}
\seealso{
\code{arpr.gb2} for details on the welfare indicators under the GB2.
}
\examples{
# GB2 parameters
a <- 3.9
b <- 18873
p <- 0.97
q <- 1.03

# Proportions defining the component densities
p0 <- rep(1/3,3)

# Mixture probabilities
pl <- c(0.39,0.26,0.35)

# for the right discretization
arpt <- arpt.cgb2(0.6, a, b, p, q, p0, pl)
arpr <- arpr.cgb2(0.6, a, p, q, p0, pl)
rmpg <- rmpg.cgb2(arpr, a, p, q, p0, pl)
qsr <- qsr.cgb2(a, p, q, p0, pl)

# for the left discretization
arptleft <- arpt.cgb2(0.6, a, b, p, q, p0, pl, "l")
}
\keyword{distribution}


