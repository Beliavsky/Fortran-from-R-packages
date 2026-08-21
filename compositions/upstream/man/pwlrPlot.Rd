\name{pwlrPlot}
\alias{pwlrPlot}
\alias{pwlrplot}
\title{Plots of pairwise logratio against a covariable. }
\description{
  Creates a matrix of plots, with each pairwise logratio against a covariable. The covariable
  can be numeric or factor, and play the role of X or Y axis.
}
\usage{
pwlrPlot(x,y,...,add.line=FALSE,line.col=2,add.robust=FALSE,rob.col=4)
}
\arguments{
  \item{x}{a vector, a column of a data.frame, or an acomp representing the first set
          of things to be displayed. Either \code{x} or \code{y} must be an 
          \code{\link{acomp}} object, and the other must be a covariable. 
          Both factors and continuous covariables allowed here.}
  \item{y}{a vector, a column of a data.frame, or an acomp representing the first set
          of things to be displayed. Either \code{x} or \code{y} must be an 
          \code{\link{acomp}} object, and the other must be a covariable. 
          Factors to be used here with caution.}
  \item{\dots}{furter parameters to the panel function}
  \item{add.line}{logical, to control the addition of a regression line
         in each panel. Ignored if covariable is a factor.}
  \item{line.col}{in case the regression line is added, which color should be
        used? Defaults to red. }
  \item{add.robust}{logical, to control the addition of a robust regression line
         in each panel. Ignored if covariable is a factor. This is nowadays 
         based on \code{\link[robustbase]{lmrob}}, but this can change in the future.}
  \item{rob.col}{in case the robust regression line is added, which color should be
        used? Defaults to blue. }
  }
\details{ This function generates a matrix of plots of all possible pairwise
logratios of the \code{\link{acomp}} argument, plotted against a covariable. The 
covariable can be a factor or a numeric vector, or a column of a matrix or data.frame.
Covariable and composition can both be represented in X or Y axis:
a factor on X axis generates a \code{\link[graphics]{boxplot}}; a factor on Y axis generates a
\code{\link[graphics]{spineplot}}; if the covariable is numeric, a default scatterplot is generated.
All dot arguments are passed to these plotting functions. In any of these cases, the diagram 
shows the logratio of the component in the row divided by the component in
the column. In the case of a numeric covariable, both classical and
robust regression lines can be added.
}
\references{
  Boogaart, K.G. v.d. , R. Tolosana (2008) Mixing Compositions and Other
  scales, Proceedings of CodaWork 08.
  
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork05/}

  \url{https://ima.udg.edu/Activitats/CoDaWork08/}

}
\author{Raimon
Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{plot.aplus}}, \code{\link{pairwisePlot}}, \code{\link[graphics]{boxplot}}, \code{\link[graphics]{spineplot}}, \code{\link[graphics]{plot.default}}
}
\examples{

data(Hydrochem)
xc = acomp(Hydrochem[,c("Ca","Mg","Na","K")])
fk = Hydrochem$River
pH = -log10(Hydrochem$H)
## x=acomp, y=factor 
pwlrPlot(xc, fk, border=2:5)
## x=factor, y=acomp
pwlrPlot(fk,xc, col=2:5)
## x=acomp, y=numeric, with colors by river
pwlrPlot(xc, pH, col=as.integer(fk)+1)
## x=numeric, y=acomp, with line
pwlrPlot(pH, xc, add.robust=TRUE)

}
\keyword{hplot}
