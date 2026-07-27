\name{plot.SVDNF}
\alias{plot.SVDNF}
\alias{plot.DNFOptim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DNF Filtering Distribution Plot Function}
\description{
This function plots the median of the filtering and prediction distributions estimated from the \code{DNF} along with user-selected upper and lower percentiles.
}
\usage{
\method{plot}{SVDNF}(x, lower_p = 0.05, upper_p = 0.95, tlim = 'default', type = "l",
  location = 'topright', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \code{SVDNF} object. The \code{plot} plots the median and selected percentiles from the filtering distribution.}
  \item{lower_p}{Lower percentile of the filtering distribution to plot.}
  \item{upper_p}{Upper percentile of the filtering distribution to plot.}
  \item{tlim}{The \code{tlim} argument gives the range over which to plot the filtering and prediction distributions are displayed.

For example to plot the first \code{500} steps, set \code{tlim = c(1, 500)}. By default, filtering and prediction distribution estimates for every step in the time-series are generated.
  
If \code{tlim} is set to a single number (e.g., \code{tlim = c(5)}), \code{plot} graphs the estimated probability density functions of the filtering (in magenta) and prediction (in blue) distributions at that timestep.

If the data are passed in an \code{xts} object, standard date-based subsetting can be used for subsetting (e.g., setting \code{tlim = "2005"} to get the filtering and prediction distributions for that year).
}
  \item{location}{Location keyword passed to the \code{legend} function to determine the location of the legend. The keyword should be selected from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"}, and \code{"center"}.}
  \item{type}{Line type argument passed to the \code{plot} function. See \code{help(plot.default)} for more details.}
  \item{...}{Other parameters to be passed through to function.}
}
\value{
No return value; this function generates two plots.

The first has the median of the volatility factor obtained from the prediction distribution as well as its upper and lower percentiles from \code{lower_p} and \code{upper_p}.

The second has the median of the volatility factor obtained from the filtering distribution as well as its upper and lower percentiles from \code{lower_p} and \code{upper_p}.
}
\examples{
set.seed(1)
# Generate 500 returns from the Bates model.
Bates_mod <- dynamicsSVM(model = "Bates")
Bates_sim <- modelSim(t = 500, dynamics = Bates_mod)

# Runs DNF on the data.
dnf_filter <- DNF(data = Bates_sim$returns, dynamics = Bates_mod)

# Plot whole interval (default)
plot(dnf_filter, ylim = c(0, 0.15),
ylab = "Volatility Factor", xlab = 'Time')

# Plot specific interval
tlim <- c(100,350)
plot(dnf_filter, ylim = c(0, 0.15),
  ylab = "Volatility Factor", xlab = 'Time', tlim = tlim)

# Plot specific point 
tlim <- c(100)
plot(dnf_filter, ylim = c(0, 0.15), type = 'l',
  ylab = "Volatility Factor", xlab = 'Time', tlim = tlim)

}