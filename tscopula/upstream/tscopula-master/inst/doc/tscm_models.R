## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(tscopula)

## ----standard-arma, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
set.seed(13)
data1 <- 0.5 + 2*arima.sim(list(ar =0.95, ma =-0.85), 1000)
ts.plot(data1)

## ----standard-arma-fit--------------------------------------------------------
copspec <- armacopula(pars = list(ar =0.01, ma =0.01))
margspec <- margin("norm")
fullspec <- tscm(copspec, margspec)
modfit <- fit(fullspec, data1, method = "full")
modfit

## ----standardarma-plot, fig.show='hold', dev.args =list(pointsize=9)----------
plot(modfit, plottype = "residual")
plot(modfit, plottype = "kendall")
plot(modfit, plottype = "margin")

## ----tscm-dvine---------------------------------------------------------------
copmod <- dvinecopula2(family = "joe",
                       kpacf = "kpacf_arma",
                       pars = list(ar = 0.9, ma = -0.8),
                       maxlag = 20)
vcopmod <- vtscopula(copmod,
                     Vtransform = V2p(delta = 0.5, kappa = 2))
margmod <- margin("slaplace",
                  pars = c(mu = 1, scale = 2, gamma = 0.7))
tscmmod <- tscm(vcopmod, margmod)
tscmmod

## ----tscm-dvine-sim, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
set.seed(13)
data2 <- sim(tscmmod, n= 2000)
hist(data2)
ts.plot(data2)

## ----margin-fit---------------------------------------------------------------
margfit <- fit(margmod, data2)

## ----stepwise-fit-------------------------------------------------------------
tscmfit_step <- fit(tscmmod, data2)
tscmfit_step
coef(tscmfit_step)
coef(tscmmod)

## ----full-optimization--------------------------------------------------------
tscmfit_full <- fit(tscmfit_step, data2, method = "full")
tscmfit_full

## ----compare-models-----------------------------------------------------------
AIC(margfit, tscmfit_step, tscmfit_full)

## ----plots1, fig.show='hold', dev.args =list(pointsize=9)---------------------
plot(tscmfit_full)
plot(tscmfit_full, plottype = "kendall")

## ----plot2, fig.show='hold', dev.args =list(pointsize=9)----------------------
plot(tscmfit_full, plottype = "margin")

## ----plots3, fig.show='hold', dev.args =list(pointsize=9)---------------------
plot(tscmfit_full, plottype = "vtransform")
plot(tscmfit_full, plottype = "volprofile")

## ----plot4, fig.show='hold', dev.args =list(pointsize=9)----------------------
plot(tscmfit_full, plottype = "volproxy")

