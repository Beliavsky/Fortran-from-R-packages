\name{plot.potd}
\alias{plot.potd}
\title{Plot Fitted POT Model}
\description{
  The plot method plot.potd provides seven different plots 
  for assessing fitted POT model. The user selects the plot
  type from a menu. See the examples below. 
}
\usage{
\method{plot}{potd}(x, \dots)
}
\arguments{
  \item{x}{an object returned by the function \code{pot}}
  \item{\dots}{other graphics parameters}
}
\details{
  Plot 1 displays the exceedance process of the chosen threshold. 
  Plots 2-4  assess the Poisson nature of the exceedance process 
  by looking at the scaled gaps between exceedances, which should
  be iid unit exponentially distributed. Plots 5-6 assess the GPD
  nature of the excesses by looking at suitably defined residuals,
  which should again be iid unit exponentially distributed. Option
  8 allows the user to call GPD plotting functions. 
}
\value{
  If plot 1 or 2 from the GPD plots is selected as the final plot
  (i.e. option 8 is selected, followed by option 1 or 2), a list
  object containing details of the plot is returned invisibly.
  This object should be used as the first argument of \code{gpd.q}
  or \code{gpd.sfall} to add quantile estimates or expected shortfall
  estimates to the plot.
}
\seealso{\code{\link{gpd}}, \code{\link{pot}},
  \code{\link{plot.gpd}}}

\examples{
data(danish)
out <- pot(danish,10) 
\dontrun{plot(out)} 
  
\dontrun{Make a plot selection (or 0 to exit):} 
\dontrun{1: plot: Point Process of Exceedances} 
\dontrun{2: plot: Scatterplot of Gaps} 
\dontrun{3: plot: Qplot of Gaps} 
\dontrun{4: plot: ACF of Gaps} 
\dontrun{5: plot: Scatterplot of Residuals} 
\dontrun{6: plot: Qplot of Residuals} 
\dontrun{7: plot: ACF of Residuals} 
\dontrun{8: plot: Go to GPD Plots}   
}
\keyword{hplot}
