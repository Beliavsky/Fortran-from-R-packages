## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(tscopula)

## ----arma11setup, fig.width = 6, fig.height = 3, fig.show='hold', dev.args =list(pointsize=9)----
vtarma11 <- vtscopula(armacopula(list(ar = 0.95, ma = -0.85)),
  Vtransform = Vlinear(delta = 0.6))
vtarma11
set.seed(19)
data <- sim(vtarma11, 2000)
ts.plot(data)

## ----arma11fit, fig.show='hold', dev.args =list(pointsize=9)------------------
vtarma11spec <- vtscopula(armacopula(list(ar = 0, ma = 0)), Vtransform = Vlinear(delta = 0.6))
vtarma11fit <- fit(vtarma11spec, data)
vtarma11fit

plot(vtarma11fit)
plot(vtarma11fit, plottype = "vtransform")
plot(vtarma11fit, plottype = "kendall" )

## ----profileplot, fig.show='hold', dev.args =list(pointsize=9)----------------
profilefulcrum(data, tscopula = vtarma11spec, locations = seq(from = 0, to = 1, length = 11))
abline(v = 0.6)

## ----vtdvinecopula2-----------------------------------------------------------
copmod <- dvinecopula2(family = "joe",
                       kpacf = "kpacf_arma",
                       pars = list(ar = 0.9, ma = -0.85),
                       maxlag = 20)
vcopmod <- vtscopula(copmod,
  Vtransform = V2p(delta = 0.6, kappa = 0.8)
)
vcopmod

## ----showvtdvinedata, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
set.seed(13)
data2 <- sim(vcopmod, n = 2000)
hist(data2)
ts.plot(data2)

## ----fitvtdvinecopula2--------------------------------------------------------
copspec_Joe <- dvinecopula2(family = "joe",
                            pars = list(ar = 0, ma = 0),
                            maxlag = 30)
vcopspec <- vtscopula(copspec_Joe, V2p(delta = 0.6))
vcopfit <- fit(vcopspec, data2, 
               tsoptions = list(hessian = TRUE),
               control = list(maxit = 1000))
vcopfit
coef(vcopfit)
coef(vcopmod)

## ----finalplots, fig.show='hold', dev.args =list(pointsize=9)-----------------
plot(vcopfit, plottype = "vtransform")
plot(vcopfit, plottype = "kendall")
plot(vcopfit, plottype = "residual")

