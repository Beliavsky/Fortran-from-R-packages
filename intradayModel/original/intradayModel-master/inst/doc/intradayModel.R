## ---- echo=FALSE, warning=FALSE-----------------------------------------------
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.retina = 2,
  out.width = "100%",
  dpi = 96 ,
  pngquant = "--speed=1"
)
knit_hooks$set(pngquant = hook_pngquant)  # brew install pngquant

## ---- message = FALSE---------------------------------------------------------
library(intradayModel)
data(volume_aapl)
volume_aapl[1:5, 1:5] # print the head of data

volume_aapl_training <- volume_aapl[, 1:104]
volume_aapl_testing <- volume_aapl[, 105:124]

## -----------------------------------------------------------------------------
model_fit <- fit_volume(volume_aapl_training)

## ---- out.width="100%"--------------------------------------------------------
analysis_result <- decompose_volume(purpose = "analysis", model_fit, volume_aapl_training)

# visualization
plots <- generate_plots(analysis_result)
plots$log_components

## ---- out.width="100%"--------------------------------------------------------
forecast_result <- forecast_volume(model_fit, volume_aapl_testing)

# visualization
plots <- generate_plots(forecast_result)
plots$original_and_forecast

## ---- warning=FALSE-----------------------------------------------------------
data(volume_fdx)
head(volume_fdx)
tail(volume_fdx)

## -----------------------------------------------------------------------------
# set fixed value
fixed_pars <- list()
fixed_pars$"x0" <- c(13.33, -0.37)

# set initial value 
init_pars <- list()
init_pars$"a_eta" <- 1

volume_fdx_training <- volume_fdx['2019-07-01/2019-11-30']
model_fit <- fit_volume(volume_fdx_training, verbose = 2, control = list(acceleration = TRUE))

## -----------------------------------------------------------------------------
analysis_result <- decompose_volume(purpose = "analysis", model_fit, volume_fdx_training)

str(analysis_result)

## -----------------------------------------------------------------------------
plots <- generate_plots(analysis_result)
plots$log_components
plots$original_and_smooth

## -----------------------------------------------------------------------------
# use training data for burn-in
forecast_result <- forecast_volume(model_fit, volume_fdx, burn_in_days = 105) 

str(forecast_result)

## -----------------------------------------------------------------------------
plots <- generate_plots(forecast_result)
plots$log_components
plots$original_and_forecast

