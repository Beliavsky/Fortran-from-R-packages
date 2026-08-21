\name{replot}
\alias{replot}
\alias{replotable}
\alias{noreplot}
\title{Modify parameters of compositional plots.}
\description{
  Display only a subset of the plots.
}
\usage{
replot(...,dev=dev.cur(),plot=TRUE,envir=NULL,add=FALSE)
replotable(expr,add=FALSE)
noreplot(expr,dev=dev.cur())
}
\arguments{
  \item{expr}{A (unquoted) expression that does the
    plotting. \code{replotable} will make the generated plot replotable
    and noreplot will do the inverse and avoid that the plots overwrites
    the current database entry.}
  \item{...}{Plot parameters to be modified. E.g. onlyPanel}
  \item{dev}{The device that currently contains the plot. It will be
    plotted in the current device.}
  \item{plot}{ logical or call or list of calls.
    If plot is TRUE, the new version of the plot is plotted in the
    current environment (and typically stores itself here).\cr
    If plot is FALSE the modified plot is simply stored, rather
    than actually plotted (in its own old plotting environment).\cr
    If the parameter is something else, it is stored to the internal plot
    database for the given device (but not plotted or evaluated).
  }
  \item{envir}{a new enviroment to be assigned to the plot. Rarely
    needed.}
  \item{add}{either a logical to indicating that the plot adds something
    to the plot. Or a number / name of the added thing to be modified.}
}
\value{
  \code{replot} returns an invisible copy of the modified
  call. \code{replotable} and \code{noreplot} return the result of
  expression. 
}
\details{
  Some of the plot routines of compositions internally store their
  call as a mean for replaying the plot when information is added or
  parameters are modified. The stored call can be modified by this
  function, which pretty much works like a simplified version of
  \code{\link{update}}.

  \code{replot} allows to redo the plot typically in a different
  device or with different parameters. The function provides this
  functionallity at a totally different level than
  \code{\link{dev.copy}} and allows for the modification of high level
  parameters on the fly. 

  Plots can be stored in the internal database by calling \code{replot}
  with a parameter \code{plot} set to the call of that plot. Plotting
  functions without this functionallity can be filtered through
  replotable(). However in this case all parameter names should be given
  explicitly.   

  There are actually three levels of possible replay: The \code{dev.copy}
  level on which graphic actions are replayed. The \code{gsi.pairs}
  function level that organizes panels plots and uses an internal replotting facility
  to allow modification of the parameter, e.g. addings lines .... And than
  there is the high level of the actual function call generating the
  plot. 


} 
\note{
  The function works by revaluating the call in its environment. Thus
  the plot will change!!! if the data has changed.

  The function always handles the latest plot from the package. If
  another plot ignorant of the replot system has meanwhile be used it will be
  ignored. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{plot.acomp}}, \code{\link{plot.aplus}},
  \code{\link{boxplot.acomp}}
}
\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals5))
straight(acomp(c(1,1,1,1,1)),acomp(c(1,2,3,4,5)))
replot(onlyPanel=c(2,3))
oldPlot <- replot(plot=FALSE) # get the plotting call
replotable(plot(x=1:10)) # To make a graphic replottable
replot(col=1:10)
replot(plot=oldPlot)    # Restore the old plot (without plotting)
replot(onlyPanel=NULL)    # View the whole plot again
replot(pch=20)  # Acctually plot it 
replot(col=20)  # since the actual plot is gsi.pairs not a plot.acomp
 
\dontrun{
# The following line in  a plotting function stores the plot for replotting.
replot(plot=match.call()) # Store current call as plot 
replot()                  # simply plot once again
replot(dev=otherdev)      # redo a plot from an other device here.
replot(onlyPanel=c(3,4))  # modify the plot (and replot it)
replot(onlyPanel=c(3,4),dev=7,plot=FALSE) # modify a stored plot
}
}
\keyword{aplot}
