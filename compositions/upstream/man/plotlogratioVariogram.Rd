\name{plotlogratioVariogram}
\alias{plot.logratioVariogram}
\title{Empirical variograms for compositions}
\description{
  Plots a logratioVariogram.
}
\usage{
\method{plot}{logratioVariogram}(x,...,type="l",lrvg=NULL,
                fcols=2:length(lrvg),oma=c(4, 4, 4, 4),gap=0,ylim=NULL)
          }
\arguments{
  \item{x}{The logratioVariogram created by
    \code{\link{logratioVariogram}}}
  \item{\dots}{further parameters for \code{plot.default}}
  \item{type}{as in \code{plot.default}}
  \item{lrvg}{a model function for a logratiovariogram or a list of
    several, to be added to the plot.}
  \item{fcols}{the colors for the different lrvg variograms}
  \item{oma}{The outer margin of the paneled plot}
  \item{gap}{The distance of the plot panals used to determin \code{mar}}
  \item{ylim}{The limits of the Y-axis. If zero it is automatically
    computed.}
}
\value{
  Nothing. 
}
\details{
  see \code{\link{logratioVariogram}}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vgram2lrvgram}}, \code{\link{CompLinModCoReg}}
}

\examples{
\dontrun{
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
lrv <- logratioVariogram(comp,X,maxdist=1,nbins=10)
fff <- CompLinModCoReg(~nugget()+sph(0.5)+R1*exp(0.7),comp)
fit <- vgmFit(lrv,fff)
fit
fff(1:3)
plot(lrv,lrvg=vgram2lrvgram(fit$vg))
}
}
\keyword{multivariate}
