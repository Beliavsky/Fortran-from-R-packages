## ----setup, message = FALSE---------------------------------------------------
library(nvmix)
library(RColorBrewer)
doPDF <- FALSE

## ----fig.align = "center", fig.width = 7, fig.height = 7, fig.show = "hold"----
set.seed(1) # for reproducibility
qmix  <- function(u, df) 1/qgamma(1-u, shape = df/2, rate = df/2)
df    <- 3.5 
n     <- 20
level <- seq(from = 0.9, to = 0.995, length.out = n)
VaR_true <- VaR_nvmix(level, qmix = "inverse.gamma", df = df)
VaR_est  <- VaR_nvmix(level, qmix = qmix, df = df)
stopifnot(all.equal(VaR_true, qt(level, df = df)))
## Prepare plot
pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
cols <- pal(2) # colors
if(doPDF) pdf(file = (file <- "fig_VaR_nvmix_comparison.pdf"),
              width = 7, height = 7)
plot(NA, xlim = range(level), ylim = range(VaR_true, VaR_est), 
     xlab = expression(alpha), ylab = expression(VaR[alpha]))
lines(level, VaR_true, col = cols[1], lty = 2, type = 'b')
lines(level, VaR_est,  col = cols[2], lty = 3, lwd = 2)
legend('topleft', c("True VaR", "Estimated VaR"), col = cols, lty = c(2,3),
       pch = c(1, NA))
if(doPDF) dev.off()

## ----fig.align = "center", fig.width = 7, fig.height = 7, fig.show = "hold"----
ES_true <- ES_nvmix(level, qmix = "inverse.gamma", df = df)
ES_est  <- ES_nvmix(level, qmix = qmix, df = df)
## Prepare plot
if(doPDF) pdf(file = (file <- "fig_ES_nvmix_comparison.pdf"),
              width = 7, height = 7)
plot(NA, xlim = range(level), ylim = range(ES_true, ES_est), 
     xlab = expression(alpha), ylab = expression(ES[alpha]))
lines(level, ES_true, col = cols[1], lty = 2, type = 'b')
lines(level, ES_est,  col = cols[2], lty = 3, lwd = 2)
legend('topleft', c("True ES", "Estimated ES"), col = cols, lty = c(2,3),
       pch = c(1, NA))
if(doPDF) dev.off()

