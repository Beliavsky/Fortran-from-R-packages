\name{get_set_qqplot_param}
\alias{get_set_qqplot_param}
\title{Plotting parameters for QQ Plots}
\description{
  Plotting parameters for the method \code{plot()} of the class 
  \code{qqplot_maha}.
}
\usage{
get_set_qqplot_param(plot.pars = list(log = ""))
}
\arguments{
  \item{plot.pars}{\code{\link{list}} specifying plotting
    parameters to be set; see below under details.}
}
\value{
  \code{get_set_qqplot_param()} returns a \code{\link{list}} with 13 elements
  that is passed to \code{\link{qqplot_maha}()}, more specifically, to the
  underlying \code{plot()} method. 
  Parameter values passed to \code{get_set_qqplot_param()} via the
  \code{plot.pars} argument overwrite the defaults; for parameters not
  specified in the \code{plot.pars} argument, the default values are being
  returned.
}
\details{
  This function provides a convenient way to set plotting parameters in the 
  argument \code{plot.pars} of the 
  function \code{\link{qqplot_maha}()} (more precisely, the underlying
  \code{plot()} method), such as logarithmic plotting, colors, linetypes and more. 
  
  The input list \code{plot.pars} can contain any of the following:

  \describe{
    \item{\code{log}}{\code{\link{character}} specifying the logarithmic 
    axes. Just like for the generic \code{plot}, must be one of \code{""},
    \code{"x"}, \code{"y"} or \code{"xy"}. }
    \item{\code{xlim, ylim}}{The x- and y-limits for plotting.}
    \item{\code{xlab, ylab}}{\code{\link{character}} specifying the x- and 
    y-axis labels. Default to \code{"Theoretical quantiles"} and 
    \code{"Sample quantiles"}, respectively.}
    \item{\code{sub, main}}{\code{\link{character}} specifying title and subtitle
    of the plot; default to \code{""}, so no titles.}
    \item{\code{plot_legend, plot_test, plot_line}}{\code{\link{logical}} 
    specifying if a legend should be plotted; if the test result of 
    the GoF test should be displayed on the 3rd axis and if the plot should 
    contain a fitted line. All default to \code{TRUE}.}
    \item{\code{pch}}{specification of the plotting symbol; see 
    \code{?\link{points}()}. Defaults to \code{1}.}
    \item{\code{lty}}{3-\code{\link{vector}} containing the specification of
    the linetypes for i) the diagonal, ii) the asymptotic CI and iii)
    the bootstrap CI; see also \code{?\link{par}()}. Defaults to 
    \code{1:3}.}
    \item{\code{col}}{4-\code{\link{vector}} specifying the colors to be used 
    for i) the points in the QQ plot; ii) the diagonal; iii) the asymptotic
    CI and iv) the bootstrap CI. Defaults to \code{c("black", "red", "azure4", 
    "chocolate4")}.}
  }
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\seealso{
  \code{\link{qqplot_maha}()}}
\examples{
get_set_qqplot_param(plot.pars = list()) # obtain defaults

## See ?qqplot_maha() for more examples.
}
\keyword{distribution}