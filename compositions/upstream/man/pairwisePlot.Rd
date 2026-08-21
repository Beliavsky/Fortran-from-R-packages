\name{pairwiseplot}
\alias{pairwisePlot}
\alias{pairwisePlot.default}
\title{Creates a paneled plot like pairs for two different datasets. }
\description{
  Creates a plot for each element of two lists or each column of each
  dataset against each of the second.
}
\usage{
pairwisePlot(X,Y,...)
\method{pairwisePlot}{default}(X,Y=X,...,
                xlab=deparse(substitute(X)),ylab=deparse(substitute(Y)),
                nm=c(length(Y),length(X)),panel=plot,
                add.line=FALSE, line.col=2,add.robust=FALSE,rob.col=4)}
\arguments{
  \item{X}{a list, a data.frame, or a matrix representing the first set
          of things to be displayed.}
  \item{Y}{a list, a data.frame, or a matrix representing the second set
          of things to be displayed.}
  \item{\dots}{furter parameters to the panel function}
  \item{xlab}{The sequence of labels for the elements of
          X. Alternatively the labels can be given as colnames or names
          of X. This option takes precedence if specified.}
  \item{ylab}{The sequence of labels for the elements of
          Y. Alternatively the labels can be given as colnames or names
          of Y. This option takes precedence if specified.}
  \item{nm}{the parameter to be used in the call
           \code{par(mfrow=nm)}. If NULL no parameter is setted and a
       sequence of plots can be generated. }
  \item{panel}{The panel function to plot the individual panels.
    If the panel function admits a formula interface, it is called
    as \code{panel(y~x, xlab=xlab,ylab=ylab,...)}, otherwise
    as \code{panel(x, y,xlab=xlab,ylab=ylab,...)}. Thus the
    panel function must be capable of taking these arguments. It must
    also set up its own plot. There is no negotiation on coordinate
    system.     
   }
  \item{add.line}{logical, to control the addition of a regression line
         in each panel}
  \item{line.col}{in case the regression line is added, which color should be
        used? defaults to red. }
  \item{add.robust}{logical, to control the addition of a robust regression line
         in each panel. Ignored if covariable is a factor. This is nowadays 
         based on \code{\link[robustbase]{lmrob}}, but this can change in the future.}
  \item{rob.col}{in case the robust regression line is added, which color should be
        used? Defaults to blue. }   
  }
\details{ This is a light-weight convenience function to plot several
aspects of one dataset against several aspects of another dataset.  It
is far more straight-forward than e.g. the \code{pairs} function and
does not do any internal computation rather than organizing the names.
Of course, the rows of the two data sets must be the same.

The current implementation may display a warning about the function
\code{panel} dispatching methods for generic \code{plot}. It can be
ignored without harm.

Optionally, classical and/or robust regression lines can be drawn, though
only for non-factor covariables.

It may be convenient to use \code{\link[graphics]{par}} capabilities to fit the device
characteristics to the plot, in particular arguments \code{mar} and \code{oma}.
}
\references{
  Boogaart, K.G. v.d. , R. Tolosana (2008) Mixing Compositions and Other
  scales, Proceedings of CodaWork 08.

  
  \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr
  
  \url{https://ima.udg.edu/Activitats/CoDaWork05/}

  \url{https://ima.udg.edu/Activitats/CoDaWork08/}

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{plot.aplus}}, \code{\link{balance}}, \code{\link{pwlrPlot}}
}
\examples{

X <- rnorm(100)
Y <- rnorm.acomp(100,acomp(c(A=1,B=1,C=1)),0.1*diag(3))+acomp(t(outer(c(0.2,0.3,0.4),X,"^")))

pairs(cbind(ilr(Y),X),panel=function(x,y,...) {points(x,y,...);abline(lm(y~x))})
pairs(cbind(balance(Y,~A/B/C),X),
         panel=function(x,y,...) {points(x,y,...);abline(lm(y~x))})
pairwisePlot(balance(Y,~A/B/C),X)
pairwisePlot(X,balance(Y,~A/B/C),
         panel=function(x,y,...) {plot(x,y,...);abline(lm(y~x))})
pairwisePlot(X,balance01(Y,~A/B/C))

# A function to extract a portion representation of subcompsitions
# with two elements:
subComps <- function(X,...,all=list(...)) {
  X <- oneOrDataset(X)
  nams <- sapply(all,function(x) paste(x[[2]],x[[3]],sep=","))
  val  <- sapply(all,function(x){ 
             a = X[,match(as.character(x[[2]]),colnames(X)) ]
             b = X[,match(as.character(x[[2]]),colnames(X)) ]
             c = X[,match(as.character(x[[3]]),colnames(X)) ]
             return(a/(b+c))
           })
  colnames(val)<-nams
  val
}

pairwisePlot(X,subComps(Y,A~B,A~C,B~C))

## using Hydrochemical data set as illustration of mixed possibilities
data(Hydrochem)
xc = acomp(Hydrochem[,c("Ca","Mg","Na","K")])
fk = Hydrochem$River
pH = -log10(Hydrochem$H)
covars = data.frame(pH, River=fk)
pairwisePlot(clr(xc), pH)
pairwisePlot(clr(xc), pH, col=fk)
pairwisePlot(pH, ilr(xc), add.line=TRUE)
pairwisePlot(covars, ilr(xc), add.line=TRUE, line.col="magenta")
pairwisePlot(clr(xc), covars, add.robust=TRUE)

}
\keyword{hplot}
