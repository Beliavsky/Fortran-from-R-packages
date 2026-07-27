\name{predict.DNFOptim}
\alias{predict.DNFOptim}
\alias{predict.SVDNF}
\title{
  Predict Method for \code{DNFOptim} and \code{SVDNF} Objects
}
\description{
  This function generates Monte Carlo predictions for \code{DNFOptim} objects. The function does this by sampling volatilities from the discrete nonlinear filter's last filtering distribution. Then, using these volatilities as inital values for the \code{modelSim} function, the predict method generates \code{n_sim} path and estimates the means and confidence intervals for future volatility factor and return values. 
}
\usage{
  \method{predict}{DNFOptim}(object, n_ahead = 15, new_data = NULL, n_sim = 1000, confidence = 0.95, ...)
}
\arguments{
  \item{object}{An object of class \code{"DNFOptim"}.}
  \item{n_ahead}{Number of periods ahead to predict.}
  \item{new_data}{An optional data frame of new predictor values for models with factors (e.g., the CAPM-SV model) for which we want to forecast returns and volatility values.}
  \item{n_sim}{Number of simulated paths used to estimate the future volatility factor and return means and confidence intervals.}
  \item{confidence}{Confidence level for prediction intervals. Should be between 0 and 1.}
  \item{...}{Other parameters to be passed through to function.}
}
\details{
  This function uses Monte Carlo paths simulated from the MLE dynamics obtained via a \code{DNFOptim} object to generate predictions for a specified number of periods ahead.
  It returns predicted mean volatility and return values based on simulations with confidence intervals.
}
\value{
  A list containing the following components:
  \item{volatility_pred}{A list with mean volatility values and confidence intervals. Contains the following components: 
  \code{UB_vol}: the upper bound of the confidence interval for volatility,  \code{mean_vol_pred}: the mean prediction for volatility, and 
    \code{LB_vol}: the lower bound of the confidence interval for volatility.}
  \item{ret_pred}{A list with mean return values and confidence intervals. Contains the following components:
 (1) \code{UB_ret}: the upper bound of the confidence interval for mean returns,
    (2) \code{mean_ret_pred}: the mean prediction for mean returns, and
    (3) \code{LB_ret}: the lower bound of the confidence interval for mean returns.
    }

  \item{object}{The \code{DNFOptim} object input to the \code{predict} function.}

  \item{confidence}{The specified confidence level.}
}
\seealso{
  \code{\link{DNFOptim}},
}
\examples{
set.seed(1)

# Generating return data
Taylor_mod <- dynamicsSVM(model = "Taylor", phi = 0.9,
  theta = -7.36, sigma = 0.363)
Taylor_sim <- modelSim(t = 30, dynamics = Taylor_mod, init_vol = -7.36)

# Initial values and optimization bounds
init_par <- c( 0.7, -5, 0.3)
lower <- c(0.01, -20, 0.1); upper <- c(0.99, 0, 1)

# Running DNFOptim to get MLEs
optim_test <- DNFOptim(data = Taylor_sim$returns,
  dynamics = Taylor_mod,
  par = init_par, lower = lower, upper = upper, method = "L-BFGS-B")
                
# Parameter estimates
summary(optim_test)

# Predict 5 steps ahead
preds <- predict(optim_test, n_ahead = 5)

# Plot predictions with 95 percent confidence interval
plot(preds)
}

