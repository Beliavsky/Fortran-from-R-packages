## ----echo=FALSE, message=FALSE------------------------------------------------
library(knitr)
library(poweRlaw)
options(replace.assign = FALSE)

opts_chunk$set(fig.path = "knitr_figure_poweRlaw/graphicsa-",
               cache.path = "knitr_cache_poweRlaw_a/",
               fig.align = "center",
               dev = "pdf", fig.width = 5, fig.height = 5,
               fig.show = "hold", cache = FALSE, par = TRUE,
               out.width = "0.4\\textwidth")
knit_hooks$set(crop = hook_pdfcrop)

knit_hooks$set(par = function(before, options, envir) {
  if (before && options$fig.show != "none") {
    par(mar = c(3, 4, 2, 1), cex.lab = .95, cex.axis = .9,
        mgp = c(3, .7, 0), tcl = -.01, las = 1)
  }}, crop = hook_pdfcrop)

set.seed(1)
palette(c(rgb(170, 93, 152, maxColorValue = 255),
          rgb(103, 143, 57, maxColorValue = 255),
          rgb(196, 95, 46, maxColorValue = 255),
          rgb(79, 134, 165, maxColorValue = 255),
          rgb(205, 71, 103, maxColorValue = 255),
          rgb(203, 77, 202, maxColorValue = 255),
          rgb(115, 113, 206, maxColorValue = 255)))

## ----installation, eval=FALSE-------------------------------------------------
#  install.packages("poweRlaw")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("csgillespie/poweRlaw")

## -----------------------------------------------------------------------------
library("poweRlaw")

## ----eval=FALSE---------------------------------------------------------------
#  help(package = "poweRlaw")

## ----results='hide', eval=FALSE-----------------------------------------------
#  vignette(package = "poweRlaw")

## ----results='hide', eval=FALSE-----------------------------------------------
#  browseVignettes("poweRlaw")

## ----tidy=FALSE, eval=FALSE---------------------------------------------------
#  ?displ

## ----eval=FALSE---------------------------------------------------------------
#  example(displ)

## ----eval=FALSE---------------------------------------------------------------
#  demo(package = "poweRlaw")
#  data(package = "poweRlaw")

## ----results='hide', eval=FALSE-----------------------------------------------
#  citation("poweRlaw")

## ----echo=FALSE---------------------------------------------------------------
data(bootstrap_moby)
bs = bootstrap_moby
data(bootstrap_p_moby)
bs_p = bootstrap_p_moby

## ----example_word-------------------------------------------------------------
data("moby")

## ----fitting------------------------------------------------------------------
m_m = displ$new(moby)

## ----tidy=FALSE---------------------------------------------------------------
m_m$getXmin()
m_m$getPars()

## -----------------------------------------------------------------------------
m_m$setXmin(5)
m_m$setPars(2)

## -----------------------------------------------------------------------------
(est = estimate_pars(m_m))

## ----m_m, echo=1--------------------------------------------------------------
(est = estimate_xmin(m_m))
m_m$setXmin(est)

## ----echo=FALSE, fig.width=8, fig.height=8, out.width="0.8\\textwidth"--------
par(mfrow = c(2, 2))
plot(m_m, xlab = "x", ylab = "CDF",
     pch = 21, bg = 1, cex = 0.6,
     panel.first = grid())
lines(m_m, col = 2, lwd = 2)
hist(bs$bootstraps[, 2], xlab = expression(x[min]), ylim = c(0, 1600),
     xlim = c(0, 30), main = NULL, breaks = "fd")
grid()
hist(bs$bootstraps[, 3], xlab = expression(alpha),
     ylim = c(0, 500), xlim = c(1.8, 2.1), main = NULL, breaks = "fd")
grid()
plot(jitter(bs$bootstraps[, 2], factor = 1.2), bs$bootstraps[, 3],
     xlab = expression(x[min]), ylab = expression(alpha),
     xlim = c(0, 30), ylim = c(1.8, 2.1), cex = 0.35,
     pch = 21, bg = 1, panel.first = grid())

## ----m_m, echo=2, eval=FALSE, results='hide'----------------------------------
#  (est = estimate_xmin(m_m))
#  m_m$setXmin(est)

## ----fig.keep='none'----------------------------------------------------------
## Plot the data (from xmin)
plot(m_m)
## Add in the fitted distribution
lines(m_m, col = 2)

## ----fig.keep='none'----------------------------------------------------------
dd = plot(m_m)
head(dd, 3)

## ----uncertainty, eval=FALSE--------------------------------------------------
#  bs = bootstrap(m_m, no_of_sims = 1000, threads = 1)

## -----------------------------------------------------------------------------
parallel::detectCores()

## ----fig.keep='none'----------------------------------------------------------
hist(bs$bootstraps[, 2], breaks = "fd")
hist(bs$bootstraps[, 3], breaks = "fd")

## ----fig.keep='none'----------------------------------------------------------
plot(jitter(bs$bootstraps[, 2], factor = 1.2), bs$bootstraps[, 3])

## ----do_we_have_a_power,echo=FALSE, fig.width=8, fig.height=4, out.width="\\textwidth"----
par(mfrow = c(1, 3))
hist(bs_p$bootstraps[, 2], xlab = expression(x[min]), ylim = c(0, 1600),
     xlim = c(0, 45), main = NULL, breaks = "fd")
grid()
hist(bs_p$bootstraps[, 3], xlab = expression(alpha),
     ylim = c(0, 500), xlim = c(1.80, 2.05), main = NULL, breaks = "fd")
grid()

plot(jitter(bs_p$bootstraps[, 2], factor = 1.2), bs_p$bootstraps[, 3],
     xlab = expression(x[xmin]), ylab = expression(alpha),
     xlim = c(0, 40), ylim = c(1.8, 2.05), cex = 0.35,
     pch = 21, bg = 1,
     panel.first = grid())

## ----eval=FALSE, tidy=FALSE---------------------------------------------------
#  ## This may take a while
#  ## Use the mle to estimate the parameters
#  bs_p = bootstrap_p(m_m, no_of_sims = 1000, threads = 2)

## ----distribution_objects-----------------------------------------------------
m_m = displ$new(moby)

## ----echo=FALSE---------------------------------------------------------------
m_m$setXmin(est)

## ----echo=FALSE, results='hide', message=FALSE, warning=FALSE, error=FALSE----
# Downloaded from Clausetts webiste - no longer there
blackouts = c(570, 210.882, 190, 46, 17, 360, 74, 19, 460, 65, 18.351, 25,
25, 63.5, 1, 9, 50, 114.5, 350, 25, 50, 25, 242.91, 55, 164.5,
877, 43, 1140, 464, 90, 2100, 385, 95.63, 166, 71, 100, 234,
300, 258, 130, 246, 114, 120, 24.506, 36.073, 10, 160, 600, 12,
203, 50.462, 40, 59, 15, 1.646, 500, 60, 133, 62, 173, 81, 20,
112, 600, 24, 37, 238, 50, 50, 147, 32, 40.911, 30.5, 14.273,
160, 230, 92, 75, 130, 124, 120, 11, 235, 50, 94.285, 240, 870,
70, 50, 50, 18, 51, 51, 145, 557.354, 219, 191, 2.9, 163, 257.718,
1660, 1600, 1300, 80, 500, 10, 290, 375, 95, 725, 128, 148, 100,
2, 48, 18, 5.3, 32, 250, 45, 38.5, 126, 284, 70, 400, 207.2,
39.5, 363.476, 113.2, 1500, 15, 7500, 8, 56, 88, 60, 29, 75,
80, 7.5, 82.5, 272, 272, 122, 145, 660, 50, 92, 60, 173, 106.85,
25, 146, 158, 1500, 40, 100, 300, 1.8, 300, 70, 70, 29, 18.819,
875, 100, 50, 1500, 650, 58, 142, 350, 71, 312, 25, 35, 315,
500, 404, 246, 43.696, 71, 65, 29.9, 30, 20, 899, 10.3, 490,
115, 2085, 206, 400, 26.334, 598, 160, 91, 33, 53, 300, 60, 55,
60, 66.005, 11.529, 56, 4.15, 40, 320.831, 30.001, 200) * 1000

## ----loading_data,eval=FALSE--------------------------------------------------
#  blackouts = read.table("blackouts.txt")$V1

## -----------------------------------------------------------------------------
m_bl = conpl$new(blackouts)

## ----echo=FALSE---------------------------------------------------------------
est = estimate_xmin(m_bl)
m_bl$setXmin(est)
plot(m_bl, panel.first = grid())
lines(m_bl, col = 2)

## ----clean-up, include=FALSE--------------------------------------------------
# R compiles all vignettes in the same session, which can be bad
rm(list = ls(all = TRUE))

## ----echo=FALSE, eval=FALSE---------------------------------------------------
#  ##Used to generate the figures for github
#  ppi = 50
#  png("../graphics/figure1.png", width = 6 * ppi, height = 4 * ppi, res = ppi)
#  setnicepar(mfrow = c(1, 2))
#  plot(m_m, xlab = "x", ylab = "CDF")
#  lines(m_m, col = 2, lty = 2)
#  grid()
#  plot(m_bl, xlab = "x", ylab = "CDF")
#  lines(m_bl, col = 2, lty = 2)
#  grid()
#  sink = dev.off()
#  
#  png("../graphics/figure2.png", width = 6 * ppi, height = 4 * ppi, res = ppi)
#  setnicepar(mfrow = c(1, 2))
#  hist(bs$bootstraps[, 2], xlab = expression(x[min]), ylim = c(0, 2000),
#       xlim = c(0, 30), main = NULL, breaks = "fd")
#  grid()
#  hist(bs$bootstraps[, 3], xlab = expression(alpha),
#       ylim = c(0, 500), xlim = c(1.8, 2.1), main = NULL,
#       breaks = "fd")
#  grid()
#  sink = dev.off()

