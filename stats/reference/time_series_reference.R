# Regenerate deterministic reference values for selected R stats time-series functions.
options(digits = 17)

x <- c(2, 1, 4, 3, 7, 5, 8, 6)
y <- c(1, 3, 2, 5, 4, 8, 6, 7)

print(as.numeric(filter(x, c(0.2, 0.3, 0.5), sides = 2)))
print(as.numeric(filter(x, c(0.2, 0.3, 0.5), sides = 1)))
print(as.numeric(filter(x, c(0.5, -0.2), method = "recursive")))
print(as.numeric(filter(x, c(0.5, -0.2), method = "recursive", init = c(10, 20))))

fit <- acf(x, lag.max = 4, plot = FALSE)
print(as.numeric(fit$lag))
print(as.numeric(fit$acf))

fit <- pacf(x, lag.max = 4, plot = FALSE)
print(as.numeric(fit$lag))
print(as.numeric(fit$acf))
fit <- acf(x, lag.max = 4, type = "covariance", plot = FALSE)
print(as.numeric(fit$acf))

fit <- ccf(x, y, lag.max = 3, plot = FALSE)
print(as.numeric(fit$lag))
print(as.numeric(fit$acf))
fit <- ccf(x, y, lag.max = 3, type = "covariance", plot = FALSE)
print(as.numeric(fit$acf))

print(ARMAacf(ar = c(0.6, -0.2), ma = 0.4, lag.max = 6))
print(ARMAacf(ma = c(0.5, -0.25), lag.max = 5))

innovations <- c(1, -0.5, 0.25, 0, 0.75, -1, 0.5, 0.2)
print(as.numeric(arima.sim(model = list(ar = c(0.6, -0.2), ma = 0.4),
                           n = length(innovations), innov = innovations,
                           n.start = 3, start.innov = c(0.1, -0.2, 0.3))))
for (difference_order in 1:2) {
  print(as.numeric(arima.sim(
    model = list(ar = c(0.6, -0.2), ma = 0.4,
                 order = c(2, difference_order, 1)),
    n = length(innovations),
    innov = innovations,
    n.start = 3,
    start.innov = c(0.1, -0.2, 0.3)
  )))
}

spectral_x <- c(2, 1, 4, 3, 7, 5, 8, 6, 5, 4, 7, 9, 6, 3, 2, 1)
fit <- spec.pgram(spectral_x, taper = 0, detrend = FALSE, demean = TRUE,
                  fast = FALSE, plot = FALSE)
print(fit$freq)
print(fit$spec)
print(c(df = fit$df, bandwidth = fit$bandwidth))

fit <- spec.pgram(spectral_x, spans = c(3, 5), taper = 0,
                  detrend = FALSE, demean = TRUE, fast = FALSE, plot = FALSE)
print(fit$freq)
print(fit$spec)
print(c(df = fit$df, bandwidth = fit$bandwidth))

spectral_y <- c(1, 3, 2, 5, 4, 8, 6, 7, 8, 5, 3, 6, 9, 7, 4, 2)
fit <- spec.pgram(cbind(spectral_x, spectral_y), spans = c(3, 5), taper = 0,
                  detrend = FALSE, demean = TRUE, fast = FALSE, plot = FALSE)
print(fit$freq)
print(fit$spec)
print(fit$coh)
print(fit$phase)
print(c(df = fit$df, bandwidth = fit$bandwidth))

fit <- spec.pgram(spectral_x, plot = FALSE)
print(fit$freq)
print(fit$spec)
print(c(df = fit$df, bandwidth = fit$bandwidth))

fit <- spec.pgram(c(spectral_x, 0.5), taper = 0, detrend = FALSE,
                  demean = TRUE, fast = TRUE, plot = FALSE)
print(fit$freq)
print(fit$spec)

fit <- ar(spectral_x, method = "yule-walker", order.max = 2,
          aic = FALSE, demean = TRUE)
print(fit$ar)
print(fit$var.pred)
print(fit$x.mean)
print(fit$aic)
print(fit$order)
print(fit$resid)
ar_spectrum <- spec.ar(fit, n.freq = 8, plot = FALSE)
print(ar_spectrum$freq)
print(ar_spectrum$spec)

fit <- ar(spectral_x, method = "yule-walker", order.max = 5,
          aic = TRUE, demean = TRUE)
print(fit$order)
print(fit$ar)
print(fit$var.pred)
print(fit$aic)

fit <- ar.burg(spectral_x, order.max = 2, aic = FALSE, demean = TRUE)
print(fit$ar)
print(fit$var.pred)
print(fit$aic)
print(fit$resid)
print(spec.ar(fit, n.freq = 8, plot = FALSE)$spec)

fit <- ar.burg(spectral_x, order.max = 2, aic = FALSE, demean = TRUE,
               var.method = 2)
print(fit$var.pred)
print(fit$aic)

fit <- ar.burg(spectral_x, order.max = 5, aic = TRUE, demean = TRUE)
print(fit$order)
print(fit$ar)
print(fit$var.pred)
print(fit$aic)

fit <- ar.ols(spectral_x, order.max = 2, aic = FALSE, demean = TRUE,
              intercept = TRUE)
print(fit$ar)
print(fit$var.pred)
print(fit$x.mean)
print(fit$x.intercept)
print(fit$aic)
print(fit$resid)
print(spec.ar(fit, n.freq = 8, plot = FALSE)$spec)

fit <- ar.ols(spectral_x, order.max = 5, aic = TRUE, demean = TRUE,
              intercept = TRUE)
print(fit$order)
print(fit$ar)
print(fit$var.pred)
print(fit$x.intercept)
print(fit$aic)

# Exact stationary Gaussian AR(2) reference. ar.mle() fixes arima0's
# historical delta at 0.01, so call arima0 directly to request delta = 0.
fit <- arima0(spectral_x, order = c(2, 0, 0), include.mean = TRUE,
              delta = 0)
print(fit$coef)
print(fit$sigma2)
print(fit$aic)
spectral_fit <- ar(spectral_x, order.max = 2, aic = FALSE)
spectral_fit$ar <- fit$coef[1:2]
spectral_fit$var.pred <- fit$sigma2
print(spec.ar(spectral_fit, n.freq = 8, plot = FALSE)$spec)

exact_aic <- numeric(6)
for (candidate_order in 0:5) {
  candidate <- arima0(spectral_x, order = c(candidate_order, 0, 0),
                      include.mean = TRUE, delta = 0)
  exact_aic[candidate_order + 1] <- candidate$aic
}

print(exact_aic - min(exact_aic))
print(which.min(exact_aic) - 1)

# Classical additive and multiplicative seasonal decomposition.
decomposition_x <- ts(c(12, 15, 14, 18, 13, 17, 16, 20,
                        15, 19, 18, 22, 17, 21, 20, 24), frequency = 4)
for (decomposition_type in c("additive", "multiplicative")) {
  decomposition <- decompose(decomposition_x, type = decomposition_type)
  print(decomposition$trend)
  print(decomposition$figure)
  print(decomposition$seasonal)
  print(decomposition$random)
}

custom_decomposition <- decompose(
  decomposition_x,
  filter = c(0.1, 0.2, 0.3, 0.4)
)
print(custom_decomposition$trend)
print(custom_decomposition$figure)
print(custom_decomposition$random)

odd_decomposition <- decompose(
  ts(c(8, 11, 10, 9, 12, 11, 10, 13, 12, 11, 14, 13), frequency = 3)
)
print(odd_decomposition$trend)
print(odd_decomposition$figure)

# STL with evolving, periodic, and robust seasonal smoothing.
stl_x <- ts(c(10, 13, 12, 16, 11, 15, 14, 18,
              13, 16, 15, 20, 14, 18, 17, 21,
              16, 19, 18, 23, 17, 21, 20, 24), frequency = 4)
stl_fits <- list(
  stl(stl_x, s.window = 7),
  stl(stl_x, s.window = "periodic"),
  stl(stl_x, s.window = 7, robust = TRUE, outer = 3)
)
for (stl_fit in stl_fits) {
  print(stl_fit$time.series)
  print(stl_fit$weights)
  print(stl_fit$win)
  print(stl_fit$deg)
  print(stl_fit$jump)
}

# Holt-Winters filtering with fixed smoothing parameters and initial states.
hw_x <- ts(c(12, 15, 14, 18, 13, 17, 16, 20,
             15, 19, 18, 22, 17, 21, 20, 24), frequency = 4)
for (seasonal_type in c("additive", "multiplicative")) {
  initial_season <- if (seasonal_type == "additive") {
    c(-1, 0.5, 1, -0.5)
  } else {
    c(0.9, 1.05, 1.1, 0.95)
  }
  fit <- HoltWinters(hw_x, alpha = 0.4, beta = 0.2, gamma = 0.3,
                     seasonal = seasonal_type, l.start = 10,
                     b.start = 0.5, s.start = initial_season)
  print(fit$fitted)
  print(fit$SSE)
  print(fit$coefficients)
  print(predict(fit, n.ahead = 6))
  print(predict(fit, n.ahead = 6, prediction.interval = TRUE,
                level = 0.9))

  auto_fit <- HoltWinters(hw_x, alpha = 0.4, beta = 0.2, gamma = 0.3,
                          seasonal = seasonal_type, start.periods = 2)
  print(auto_fit$fitted)
  print(auto_fit$SSE)
  print(auto_fit$coefficients)
  print(predict(auto_fit, n.ahead = 6))

  optimized_fit <- HoltWinters(hw_x, seasonal = seasonal_type,
                               start.periods = 2,
                               optim.control = list(maxit = 2000))
  print(c(optimized_fit$alpha, optimized_fit$beta, optimized_fit$gamma))
  print(optimized_fit$SSE)
  print(optimized_fit$coefficients)
  print(predict(optimized_fit, n.ahead = 6))
}

# Odd seasonal periods exercise the unweighted centered moving average.
odd_hw_x <- ts(c(8, 11, 10, 9, 12, 11, 10, 13, 12, 11, 14, 13),
               frequency = 3)
odd_fit <- HoltWinters(odd_hw_x, alpha = 0.35, beta = 0.15,
                       gamma = 0.25, seasonal = "additive",
                       start.periods = 2)
print(odd_fit$SSE)
print(odd_fit$coefficients)
print(predict(odd_fit, n.ahead = 4))

# Nonseasonal Holt and simple exponential-smoothing paths.
reduced_x <- ts(c(10, 12, 11, 14, 13, 15, 16, 15))
holt_fit <- HoltWinters(reduced_x, alpha = 0.4, beta = 0.2,
                        gamma = FALSE, l.start = 10, b.start = 0.5)
print(holt_fit$fitted)
print(holt_fit$SSE)
print(holt_fit$coefficients)
print(predict(holt_fit, n.ahead = 4))
print(predict(holt_fit, n.ahead = 4, prediction.interval = TRUE,
              level = 0.9))

simple_fit <- HoltWinters(reduced_x, alpha = 0.4, beta = FALSE,
                          gamma = FALSE, l.start = 10)
print(simple_fit$fitted)
print(simple_fit$SSE)
print(simple_fit$coefficients)
print(predict(simple_fit, n.ahead = 4))
print(predict(simple_fit, n.ahead = 4, prediction.interval = TRUE,
              level = 0.9))

# Seasonal filtering can also omit the trend recurrence.
season_only_x <- ts(as.numeric(hw_x[1:12]), frequency = 4)
season_only_fit <- HoltWinters(season_only_x, alpha = 0.4, beta = FALSE,
                               gamma = 0.3, seasonal = "additive",
                               l.start = 10,
                               s.start = c(-1, 0.5, 1, -0.5))
print(season_only_fit$fitted)
print(season_only_fit$SSE)
print(season_only_fit$coefficients)
print(predict(season_only_fit, n.ahead = 6))

auto_season_only_fit <- HoltWinters(hw_x, alpha = 0.4, beta = FALSE,
                                    gamma = 0.3, seasonal = "additive")
print(auto_season_only_fit$SSE)
print(auto_season_only_fit$coefficients)
print(predict(auto_season_only_fit, n.ahead = 4))

# Nonseasonal conditional-sum-of-squares ARIMA fitting and point forecasts.
arima_x <- c(2, 1, 4, 3, 7, 5, 8, 6, 5, 4, 7, 9,
             6, 3, 2, 1, 4, 6, 5, 7, 8, 6, 4, 3)
for (model_order in list(c(1, 0, 1), c(0, 0, 2), c(1, 1, 0))) {
  arima_fit <- arima(arima_x, order = model_order, method = "CSS",
                     optim.control = list(maxit = 4000, reltol = 1e-12))
  print(arima_fit$coef)
  print(arima_fit$sigma2)
  print(arima_fit$loglik)
  print(arima_fit$var.coef)
  print(arima_fit$n.cond)
  print(residuals(arima_fit))
  print(predict(arima_fit, n.ahead = 5)$pred)
}

# Partially and fully fixed coefficients exercise the explicit Fortran mask API.
for (method in c("CSS", "ML")) {
  fixed_fit <- arima(arima_x, order = c(1, 0, 1),
                     fixed = c(0.25, NA, NA), method = method,
                     optim.control = list(maxit = 5000, reltol = 1e-12))
  print(fixed_fit$coef)
  print(fixed_fit$sigma2)
  print(fixed_fit$loglik)
  print(fixed_fit$var.coef)

  all_fixed_fit <- arima(arima_x, order = c(1, 0, 1),
                         fixed = c(0.25, -0.1, 4.5), method = method)
  print(all_fixed_fit$coef)
  print(all_fixed_fit$sigma2)
  print(all_fixed_fit$loglik)
  print(all_fixed_fit$var.coef)
}

# Stationary and integrated regression models exercise xreg and newxreg.
regressor <- sin(0.4 * seq_along(arima_x))
future_regressor <- sin(0.4 * (length(arima_x) + seq_len(5)))
for (method in c("CSS", "ML")) {
  regression_fit <- arima(arima_x, order = c(1, 0, 0), xreg = regressor,
                          method = method,
                          optim.control = list(maxit = 5000, reltol = 1e-12))
  regression_forecast <- predict(regression_fit, n.ahead = 5,
                                 newxreg = future_regressor)
  print(regression_fit$coef)
  print(regression_fit$sigma2)
  print(regression_fit$loglik)
  print(regression_fit$var.coef)
  print(regression_forecast$pred)
  print(regression_forecast$se)
}

trend_regressor <- seq_along(arima_x)/10
future_trend <- (length(arima_x) + seq_len(5))/10
integrated_regression_fit <- arima(
  arima_x,
  order = c(0, 1, 0),
  xreg = trend_regressor,
  method = "ML",
  optim.control = list(maxit = 5000, reltol = 1e-12)
)
integrated_regression_forecast <- predict(
  integrated_regression_fit,
  n.ahead = 5,
  newxreg = future_trend
)
print(integrated_regression_fit$coef)
print(integrated_regression_fit$sigma2)
print(integrated_regression_fit$loglik)
print(integrated_regression_fit$var.coef)
print(integrated_regression_forecast$pred)
print(integrated_regression_forecast$se)

# Exact stationary Gaussian likelihood fixtures.
for (model_order in list(c(1, 0, 0), c(1, 0, 1), c(0, 0, 2))) {
  exact_fit <- arima(arima_x, order = model_order, method = "ML",
                     optim.control = list(maxit = 5000, reltol = 1e-12))
  print(exact_fit$coef)
  print(exact_fit$sigma2)
  print(exact_fit$loglik)
  print(exact_fit$var.coef)
  print(residuals(exact_fit))
  exact_forecast <- predict(exact_fit, n.ahead = 5)
  print(exact_forecast$pred)
  print(exact_forecast$se)
}

# Exact stationary likelihood with internal response gaps.
missing_arima_x <- arima_x
missing_arima_x[c(4, 11)] <- NA_real_
missing_fit <- arima(missing_arima_x, order = c(1, 0, 1), method = "ML",
                     optim.control = list(maxit = 5000, reltol = 1e-12))
missing_forecast <- predict(missing_fit, n.ahead = 5)
print(missing_fit$coef)
print(missing_fit$sigma2)
print(missing_fit$loglik)
print(residuals(missing_fit))
print(missing_fit$var.coef)
print(missing_forecast$pred)
print(missing_forecast$se)

# Exact diffuse likelihood after ordinary differencing with internal gaps.
for (model_order in list(c(0, 1, 0), c(1, 1, 0))) {
  diffuse_fit <- arima(missing_arima_x, order = model_order, method = "ML",
                       optim.control = list(maxit = 5000, reltol = 1e-12))
  diffuse_forecast <- predict(diffuse_fit, n.ahead = 5)
  print(diffuse_fit$coef)
  print(diffuse_fit$sigma2)
  print(diffuse_fit$loglik)
  print(residuals(diffuse_fit))
  print(diffuse_fit$var.coef)
  print(diffuse_forecast$pred)
  print(diffuse_forecast$se)
}

# Multiplicative seasonal CSS models. The Fortran predictor uses CSS residuals;
# R's seasonal-MA predictions use its state-space residual initialization and
# can consequently differ while the fitted CSS coefficients and objective agree.
seasonal_i <- seq_len(80)
seasonal_innovations <- sin(1.7 * seasonal_i) + 0.3 * cos(0.23 * seasonal_i)
seasonal_x <- rep(10, length(seasonal_i))
for (i in seasonal_i) {
  value <- 10 + seasonal_innovations[i]
  if (i > 1) value <- value + 0.3 * (seasonal_x[i - 1] - 10)
  if (i > 4) value <- value + 0.45 * (seasonal_x[i - 4] - 10)
  if (i > 5) value <- value - 0.135 * (seasonal_x[i - 5] - 10)
  seasonal_x[i] <- value
}
seasonal_models <- list(
  list(order = c(1, 0, 0), seasonal = c(1, 0, 0)),
  list(order = c(0, 0, 1), seasonal = c(0, 0, 1)),
  list(order = c(0, 0, 0), seasonal = c(0, 1, 0))
)
for (model in seasonal_models) {
  seasonal_fit <- arima(
    seasonal_x,
    order = model$order,
    seasonal = list(order = model$seasonal, period = 4),
    method = "CSS",
    optim.control = list(maxit = 5000, reltol = 1e-12)
  )
  print(seasonal_fit$coef)
  print(seasonal_fit$sigma2)
  print(seasonal_fit$loglik)
  print(seasonal_fit$n.cond)
  print(predict(seasonal_fit, n.ahead = 8)$pred)
}

fixed_seasonal_fit <- arima(
  seasonal_x,
  order = c(1, 0, 0),
  seasonal = list(order = c(1, 0, 0), period = 4),
  fixed = c(NA, 0.5, NA),
  method = "CSS",
  optim.control = list(maxit = 5000, reltol = 1e-12)
)
print(fixed_seasonal_fit$coef)
print(fixed_seasonal_fit$sigma2)
print(fixed_seasonal_fit$loglik)
print(fixed_seasonal_fit$var.coef)

# Exact seasonal and integrated likelihood paths.
exact_seasonal_models <- list(
  list(order = c(1, 0, 0), seasonal = c(1, 0, 0)),
  list(order = c(0, 0, 1), seasonal = c(0, 0, 1)),
  list(order = c(0, 0, 0), seasonal = c(0, 1, 0)),
  list(order = c(0, 1, 0), seasonal = c(0, 1, 0))
)
for (model in exact_seasonal_models) {
  exact_seasonal_fit <- arima(
    seasonal_x,
    order = model$order,
    seasonal = list(order = model$seasonal, period = 4),
    method = "ML",
    optim.control = list(maxit = 5000, reltol = 1e-12)
  )
  print(exact_seasonal_fit$coef)
  print(exact_seasonal_fit$sigma2)
  print(exact_seasonal_fit$loglik)
  print(residuals(exact_seasonal_fit))
  exact_seasonal_forecast <- predict(exact_seasonal_fit, n.ahead = 8)
  print(exact_seasonal_forecast$pred)
  print(exact_seasonal_forecast$se)
}

# Exact diffuse likelihood after seasonal differencing with internal gaps.
seasonal_missing <- seasonal_x
seasonal_missing[c(13, 34)] <- NA_real_
diffuse_seasonal_fit <- arima(
  seasonal_missing,
  order = c(0, 0, 0),
  seasonal = list(order = c(0, 1, 0), period = 4),
  method = "ML"
)
diffuse_seasonal_forecast <- predict(diffuse_seasonal_fit, n.ahead = 8)
print(diffuse_seasonal_fit$sigma2)
print(diffuse_seasonal_fit$loglik)
print(residuals(diffuse_seasonal_fit))
print(diffuse_seasonal_forecast$pred)
print(diffuse_seasonal_forecast$se)

diffuse_combined_fit <- arima(
  seasonal_missing,
  order = c(0, 1, 0),
  seasonal = list(order = c(0, 1, 0), period = 4),
  method = "ML"
)
diffuse_combined_forecast <- predict(diffuse_combined_fit, n.ahead = 8)
print(diffuse_combined_fit$sigma2)
print(diffuse_combined_fit$loglik)
print(diffuse_combined_forecast$pred)
print(diffuse_combined_forecast$se)
