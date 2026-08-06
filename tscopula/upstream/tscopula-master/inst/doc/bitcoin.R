## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(tscopula)
tsoptions <- list(hessian = TRUE, method = "Nelder-Mead", avoidzero= FALSE)


## ----define-data, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
data(bitcoin)
X <- (diff(log(bitcoin))[-1]) * 100 # log-returns (as percentages)
length(X)
plot(X, type = "h")

## ----arma-copula-profile, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
U <- strank(X)

copmod_Gauss <- armacopula(pars = list(ar = 0.95, ma = -0.85))
pts <- seq(from=0.01, to=0.99, length=51) + 0.005
profilefulcrum(U, copmod_Gauss, locations = pts)
abline(v = 0.45)

## ----arma-copula-fit----------------------------------------------------------
mod_Gauss <- vtscopula(copmod_Gauss, Vlinear(0.45))
fit_Gauss <- fit(mod_Gauss, U, tsoptions = tsoptions)
fit_Gauss

## ----dvinecopula2-------------------------------------------------------------
copmod_Frank <- dvinecopula2(family = "frank",
                             pars = list(ar = 0.95, ma = -0.85),
                             maxlag = 30)
mod_Frank <- vtscopula(copmod_Frank, Vlinear(0.45))
fit_Frank <- fit(mod_Frank, U, tsoptions = tsoptions)

AIC(fit_Gauss, fit_Frank)

## ----dvinecopula2-plots, fig.show='hold', fig.width = 6, fig.height = 3, dev.args =list(pointsize=9)----
plot(fit_Frank, plottype = "residual")
plot(fit_Frank, plottype = "kendall")

## ----margmods-----------------------------------------------------------------
marg_st <- fit(margin("st"), X)
marg_sst <- fit(margin("sst"), X)
marg_lp <- fit(margin("laplace", 
                      pars = c(mu = 0.2, scale = 2.7)), X)
marg_slp <- fit(margin("slaplace", 
                pars = c(mu = 0.2, scale = 2.7, gamma = 0.9)), X)
marg_dw <- fit(margin("doubleweibull", 
                      pars = c(mu = 0.2, shape = 0.8, scale = 2.7)), X)
marg_sdw <- fit(margin("sdoubleweibull", 
                pars = c(mu = 0.2, shape = 0.8, scale = 2.7, gamma = 0.9)), X)
AIC(marg_st, marg_sst, marg_slp, marg_lp, marg_dw, marg_sdw)

## ----fullmod-fit,fig.show='hold', dev.args =list(pointsize=9)-----------------
fullmod <- tscm(fit_Frank, margin = marg_dw)
fullmod <- fit(fullmod, as.numeric(X), 
               method = "full", tsoptions = tsoptions)
fullmod
AIC(marg_dw, fullmod)

## ----fullmod-plots1,fig.show='hold', dev.args =list(pointsize=9)--------------
plot(fullmod, plottype = "residual")
plot(fullmod, plottype = "kendall")
plot(fullmod, plottype = "margin")
plot(fullmod, plottype = "vtransform")
plot(fullmod, plottype = "volprofile")
plot(fullmod, plottype = "volproxy")

## ----fullmod-plots2, fig.show='hold', fig.width = 6, fig.height = 6, dev.args =list(pointsize=9)----
plot(fullmod, plottype = "glag")

