\name{plot.Lc}
\alias{plot.Lc}
\alias{lines.Lc}
\title{Plot Lorenz Curve}
\description{
 plotting method for objects of class \code{"Lc"} (Lorenz curves)}
\usage{
\method{plot}{Lc}(x, general=FALSE, lwd=2, xlab="p", ylab="L(p)",
    main="Lorenz curve", las=1, ...) }
\arguments{
 \item{x}{an object of class \code{"Lc"}}
 \item{general}{logical. If TRUE the generalized Lorenz curve will be plotted}
 \item{lwd, xlab, ylab, main, las, ...}{high-level \code{\link{plot}} function parameters.}
}
\references{B C Arnold: Majorization and the Lorenz Order: A Brief Introduction, 1987, Springer,

F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef.}

\seealso{\code{\link{Lc}}, \code{\link{Lc.mehran}},
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
## plot both Lorenz curves
plot(Lc.p)
lines(Lc.u, col=2)
}
\keyword{hplot}

