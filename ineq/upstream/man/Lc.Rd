\name{Lc}
\alias{Lc}
\title{Lorenz Curve}
\description{
 Computes the (empirical) ordinary and generalized Lorenz curve of a vector x}
\usage{
Lc(x, n = rep(1,length(x)), plot = FALSE)
}
\arguments{
 \item{x}{a vector containing non-negative elements.}
 \item{n}{a vector of frequencies, must be same length as \code{x}.}
 \item{plot}{logical. If TRUE the empirical Lorenz curve will be plotted.}
}
\details{\code{Lc(x)} computes the empirical ordinary Lorenz curve of \code{x}
  as well as the generalized Lorenz curve (= ordinary Lorenz curve *
  mean(x)). The result can be interpreted like this: \code{p}*100 percent
  have \code{L(p)}*100 percent of \code{x}.

  If \code{n} is changed to anything but the default \code{x} is
  interpreted as a vector of class means and \code{n} as a vector of
  class frequencies: in this case \code{Lc} will compute the minimal
  Lorenz curve (= no inequality within each group). A maximal curve can be
  computed with \code{\link{Lc.mehran}}.
}
\value{
  A list of class \code{"Lc"} with the following components:
\item{p}{vector of percentages}
\item{L}{vector with values of the ordinary Lorenz curve}
\item{L.general}{vector with values of the generalized Lorenz curve}
}
\references{B C Arnold: Majorization and the Lorenz Order: A Brief Introduction, 1987, Springer,

F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef.}

\seealso{\code{\link{plot.Lc}}, \code{\link{Lc.mehran}},
\code{\link{plot.theorLc}}}

\examples{
## Load and attach income (and metadata) set from Ilocos, Philippines
data(Ilocos)
attach(Ilocos)
## extract and rescale income for the provinces "Pangasinan" und "La Union"
income.p <- income[province=="Pangasinan"]/10000
income.u <- income[province=="La Union"]/10000

## compute the Lorenz curves
Lc.p <- Lc(income.p)
Lc.u <- Lc(income.u)
## it can be seen the the inequality in La Union is higher than in
## Pangasinan because the respective Lorenz curve takes smaller values.
plot(Lc.p)
lines(Lc.u, col=2)
## the picture becomes even clearer with generalized Lorenz curves
plot(Lc.p, general=TRUE)
lines(Lc.u, general=TRUE, col=2)
## inequality measures emphasize these results, e.g. Atkinson's measure
ineq(income.p, type="Atkinson")
ineq(income.u, type="Atkinson")
## or Theil's entropy measure
ineq(income.p, type="Theil", parameter=0)
ineq(income.u, type="Theil", parameter=0)



# income distribution of the USA in 1968 (in 10 classes)
# x vector of class means, n vector of class frequencies
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
n <- c(482, 825, 722, 690, 661, 760, 745, 2140, 1911, 1024)

# compute minimal Lorenz curve (= no inequality in each group)
Lc.min <- Lc(x, n=n)
# compute maximal Lorenz curve (limits of Mehran)
Lc.max <- Lc.mehran(x,n)
# plot both Lorenz curves in one plot
plot(Lc.min)
lines(Lc.max, col=4)

# add the theoretic Lorenz curve of a Lognormal-distribution with variance 0.78
lines(Lc.lognorm, parameter=0.78)
# add the theoretic Lorenz curve of a Dagum-distribution
lines(Lc.dagum, parameter=c(3.4,2.6))
}
\keyword{misc}

