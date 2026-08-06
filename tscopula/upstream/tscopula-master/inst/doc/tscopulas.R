## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(tscopula)

## ----setup-tscopula-----------------------------------------------------------
ar1 <- armacopula(list(ar = 0.7))
ar1

## ----simulate-copula, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
set.seed(13)
data1 <- sim(ar1, 1000)
ts.plot(data1)

## ----fit-copula---------------------------------------------------------------
ar1spec <- armacopula(list(ar = 0))
ar1fit <- fit(ar1spec, data1)
ar1fit

## ----arma11example, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
arma11 <- armacopula(list(ar = 0.95, ma = -0.85))
data2 <- sim(arma11, 1000)
ts.plot(data2)

arma11spec <- armacopula(list(ar = 0.1, ma = 0.1))
arma11fit <- fit(arma11spec, data2, tsoptions = list(hessian = TRUE))
arma11fit

## ----arma11plots, fig.show='hold'---------------------------------------------
coef(arma11fit)
res <- resid(arma11fit)
acf(res)
acf(abs(res))
plot(arma11fit)
plot(arma11fit, plottype = "kendall")
mu_t <- resid(arma11fit, trace = TRUE)
ts.plot(mu_t)

## ----kalman-------------------------------------------------------------------
head(kfilter(arma11fit@tscopula, data2))

## ----dvinecopula-setup--------------------------------------------------------
copmod <- dvinecopula(
  family = c("Clayton","Frank", "Gaussian"),
  pars = list(1.2, 2, 0.15),
  rotation = c(180,0, 0)
)
copmod

## ----dvinecopula-sim, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
set.seed(29)
data1 <- sim(copmod, n = 2000)
hist(data1)
ts.plot(data1)

## ----dvinecopula-est----------------------------------------------------------
copspec <- dvinecopula(
  family = c("Clayton","Frank", "Gaussian"),
  pars = list(0.5, 1, 0),
  rotation = c(180, 0, 0)
)
copfit <- fit(copspec, data1, 
              tsoptions = list(hessian = TRUE),
              control = list(maxit = 2000))
copfit
coef(copfit)
coef(copmod)

## ----dvinecopula-fit, fig.show='hold', dev.args =list(pointsize=9)------------
plot(copfit)
plot(copfit, plottype = "kendall")

## ----genlagplot, fig.show='hold', fig.width = 6, fig.height = 6, dev.args =list(pointsize=9)----
plot(copfit, plottype = "glag")

## ----dvinecopula2-setup-------------------------------------------------------
copmod <- dvinecopula2(family = "joe",
                       kpacf = "kpacf_arma",
                       pars = list(ar = 0.9, ma = -0.85),
                       maxlag = 20)
copmod

## ----dvinecopula2-sim, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
set.seed(13)
data1 <- sim(copmod, n = 2000)
hist(data1)
ts.plot(data1)

## ----dvinecopula2-est---------------------------------------------------------
copspec_Gauss <- dvinecopula2(family = "gauss",
                              pars = list(ar = 0, ma = 0),
                              maxlag = 20)
fitGauss <- fit(copspec_Gauss, data1)
fitGauss

copspec_Joe <- dvinecopula2(family = "joe",
                            pars = list(ar = 0, ma = 0),
                            maxlag = 20)
fitJoe <- fit(copspec_Joe, data1)
fitJoe

AIC(fitGauss, fitJoe)

## ----dvinecopula2-plot, fig.show='hold', dev.args =list(pointsize=9)----------
plot(fitJoe)
plot(fitJoe, plottype = "kendall")

