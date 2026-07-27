\name{plot.predict.DNFOptim}
\alias{plot.predict.DNFOptim}
\alias{plot.predict.SVDNF}
\title{Plot Predictions from \code{DNFOptim} or \code{SVDNF} Objects}
\usage{
  \method{plot}{predict.DNFOptim}(x, ...)
  
  \method{plot}{predict.SVDNF}(x, ...)

}
\arguments{
  \item{x}{an object of class \code{predict.DNFOptim} or  \code{predict.SVDNF}.}
  \item{...}{further arguments passed to or from other methods.}
}
\description{
  Plot predictions from a \code{DNFOptim} or \code{SVDNF} object, including volatility and return mean predictions with confidence intervals. 
}
\details{
  This function plots the volatility and return predictions with confidence intervals obtained from a \code{DNFOptim} object.
  
For the volatility plot, it displays the DNF's filtering distribution median volatility for all time points in the series and, after the last observation, plots the predicted mean volatility with its confidence interval.

For the returns plot, it displays the observed returns for all time points in the series and, after the last observation, plots the predicted mean return with its confidence interval.
}
\value{
No return value; this function generates two plots.

The first has the median volatility from the filtering distribution as well as the mean predicted volatility from Monte Carlo simulated paths with its confidence interval.

The second has the observed asset returns as well as the mean predicted returns from Monte Carlo simulated paths with its confidence interval.
}
\examples{
# Generating return data
Taylor_mod <- dynamicsSVM(model = "Taylor", phi = 0.9,
                           theta = -7.36, sigma = 0.363)
Taylor_sim <- modelSim(t = 100, dynamics = Taylor_mod, init_vol = -7.36)
 
# Run the DNF
DNF_Taylor <- DNF(dynamics = Taylor_mod, data = Taylor_sim$returns)
 
# Predict the next 10 time steps
predict_Taylor <- predict(DNF_Taylor, n_ahead = 10)

# Plot the predictions
plot(predict_Taylor)

}

