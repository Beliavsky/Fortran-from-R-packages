## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)

## ----setup--------------------------------------------------------------------
library(QuantBondCurves)

## -----------------------------------------------------------------------------
coupon.dates(maturity = "2029-10-01", analysis.date = "2027-08-01", freq = 1, convention = "F")

## -----------------------------------------------------------------------------
coupon.dates(maturity = 2, analysis.date = "2023-03-01", freq = 1)

## -----------------------------------------------------------------------------
coupon.dates(maturity = "2025-02-28", analysis.date = "2023-09-29", 
             asset.type = "Fixed Income", freq = 4,
             trade.date = "2023-09-29", coupon.schedule = "LF")

## -----------------------------------------------------------------------------
coupons(maturity = 4.08, analysis.date = "2021-02-28", coupon.rate = 0.03,
        asset.type = "IBR", daycount = "ACT/360", trade.date = "2020-02-29",
        coupon.schedule = "LF")

## -----------------------------------------------------------------------------
discount.factors(dates = c("2020-09-10", "2020-12-10", "2021-03-10"),
                 rates = c(0.01, 0.015, 0.017), analysis.date = "2010-09-01",
                 rate.type = 1, freq = 12)

## -----------------------------------------------------------------------------
valuation.bonds(maturity = "2026-06-01", coupon.rate = 0.06, rates = 0.08,
                principal = 1000, analysis.date= "2023-06-01")

## -----------------------------------------------------------------------------
valuation.bonds(maturity = "2026-01-05", coupon.rate = c(0.04,0.043,0.05,0.052),
                rates = c(0.06,0.061,0.063,0.067), principal = 1000,
                analysis.date = "2025-01-05", asset.type = "IBR", freq = 4)

## -----------------------------------------------------------------------------
bond.price2rate(maturity = "2023-01-03", analysis.date = "2021-01-03",
                price = 0.98, coupon.rate = 0.04, principal = 1,
                asset.type = "IBR", freq = 4, daycount = "ACT/365")

## -----------------------------------------------------------------------------
sens.bonds(input = "price", price = 1.02, maturity = "2023-01-03",
           analysis.date = "2020-01-03", coupon.rate = 0.04,
           principal = 1, asset.type = "FixedIncome", freq = 1, 
           rate.type = 1, daycount = "ACT/365", dirty = 1)

## -----------------------------------------------------------------------------
sens.bonds(input = "rate", price = c(0.04,0.05), maturity = "2025-02-05",
           analysis.date = "2024-06-12", coupon.rate = 0.03, 
           principal = 1, asset.type = "FixedIncome", freq = 2,
           rate.type = 0)

## -----------------------------------------------------------------------------
average.life(input = c("rate"), price = c(0.043,0.05), maturity = "2023-01-03", 
             analysis.date = "2021-01-03", coupon.rate = 0.04, principal = 1,
             asset.type = "FixedIncome", freq = 1, rate.type = 0)

## ----results='hide', message=FALSE--------------------------------------------
# The `yield.curve` input is created for the IRR's of the market assets.
yield.curve <- c(0.1233,0.1280,0.131,0.1315,0.132,0.1322,0.1325,0.1323,0.1321,0.132)
# The output terms desired are established.
nodes <- seq(0,10, by = 0.001)
# Since for TES, IRR's and coupon rates differ, `market.assets` input is required.
# Below it is constructed.
market.assets       <- matrix(NA,nrow = 10,ncol = 2)
market.assets[1,2]  <- "2020-01-03"
market.assets[2,2]  <- "2021-01-03"
market.assets[3,2]  <- "2022-01-03"
market.assets[4,2]  <- "2023-01-03"
market.assets[5,2]  <- "2024-01-03"
market.assets[6,2]  <- "2025-01-03"
market.assets[7,2]  <- "2026-01-03"
market.assets[8,2]  <- "2027-01-03"
market.assets[9,2]  <- "2028-01-03"
market.assets[10,2] <- "2029-01-03"

market.assets[1,1]  <- 0.1233
market.assets[2,1]  <- 0.1280
market.assets[3,1]  <- 0.131
market.assets[4,1]  <- 0.1315
market.assets[5,1]  <- 0.132
market.assets[6,1]  <- 0.1322
market.assets[7,1]  <- 0.1325
market.assets[8,1]  <- 0.1323
market.assets[9,1]  <- 0.1321
market.assets[10,1] <- 0.132
# Calibration
curve.calibration(yield.curve = yield.curve, market.assets = market.assets,
                  analysis.date = "2019-01-03" , asset.type = "TES", freq = 1,
                  daycount = "ACT/365", fwd = 1, nodes = nodes, approximation = "constant")

## ----echo = FALSE, message=FALSE , fig.width = 7.2, fig.height = 5, results = 'hide',fig.align='center', warning=FALSE----

yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325, 0.1320)
names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,7,10)
nodes <- seq(0,10, by = 0.001)
y1 <- curve.calibration(yield.curve = yield.curve, asset.type = "IBRSwaps",
                        analysis.date = "2019-01-03" , nodes = nodes, 
                        daycount = "ACT/365", freq = 4, fwd = 0, approximation = "constant")
x <- as.numeric(names(y1))
y2 <- curve.calibration(yield.curve = yield.curve, asset.type = "IBRSwaps",
                       analysis.date = "2019-01-03" , nodes = nodes, 
                       daycount = "ACT/365", freq = 4, fwd = 1, approximation = "constant")

df <- data.frame(x = x, y1 = y1, y2 = y2)

# Plot
ggplot(df, aes(x)) +
  geom_line(aes(y = y1, color = "Spot"), size = 1.5) +
  geom_line(aes(y = y2, color = "Forward"), size = 1.5) +
  scale_color_manual(values = c("red", "blue"), 
                     name = "", 
                     labels = c("Forward", "Spot")) +
  labs(x = "term", y = "rate", title = "Bootstrapping Curve Calibration") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.box = "vertical", 
    legend.box.background = element_rect(color = "black", size = 1.2)
  )

## ----results='hide'-----------------------------------------------------------
yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325, 0.1320)
names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,7,10)
nodes <- seq(0,10, by = 0.001)
# Calibration
curve.calibration (yield.curve = yield.curve, market.assets = NULL, 
                   analysis.date = "2019-01-03", asset.type = "IBRSwaps",
                   freq = 4, npieces = 2, fwd = 0, obj = "Rates",
                   piece.term = 3, nodes = nodes, approximation = "linear")

## ----echo = FALSE, message=FALSE , fig.width = 7.2, fig.height = 5, results = 'hide',fig.align='center', warning=FALSE----
yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325, 0.1320)
names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,7,10)
nodes <- seq(0,10, by = 0.001)
y1 <- curve.calibration (market.assets = NULL, yield.curve = yield.curve, asset.type = "IBRSwaps",
                         analysis.date = "2019-01-03" , nodes = nodes, fwd = 0,
                         freq = 4, approximation = "linear",
                         npieces = 2, piece.term = 3, obj = "Rates")
x <- as.numeric(names(y1))
y2 <- curve.calibration (market.assets = NULL, yield.curve = yield.curve, asset.type = "IBRSwaps",
                        analysis.date = "2019-01-03" , nodes = nodes, fwd = 1,
                        freq = 4, approximation = "linear",
                        npieces = 6, piece.term = c(0.44,3.75,6.76,8.79,9.60), obj = "Rates")

df <- data.frame(x = x, y1 = y1, y2 = y2)


ggplot(df, aes(x)) +
  geom_line(aes(y = y1, color = "Spot"), size = 1.5) +
  geom_line(aes(y = y2, color = "Forward"), size = 1.5) +
  scale_color_manual(values = c("red", "blue"), 
                     name = "", 
                     labels = c("Forward", "Spot")) +
  labs(x = "term", y = "rate", title = "RSS Curve Calibration") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.box = "vertical", 
    legend.box.background = element_rect(color = "black", size = 1.2)
  )

## -----------------------------------------------------------------------------
# `previous.curve` input
previous.curve <- matrix(0.04,nrow = 2,ncol = 8)
rownames(previous.curve) <- c("2014-01-01","2015-01-01")
colnames(previous.curve) <- c(0, 0.25, 0.5, 1:5)
# `serie` input
serie <- matrix(NA,nrow = 4,ncol = 6)
rownames(serie) <- c("2014-01-01","2015-01-01","2016-01-01","2017-01-01")
colnames(serie) <- c(0, 0.08333, 0.25, 0.5, 1, 2)
serie[1,1] <- 0.04
serie[1,2] <- 0.05
serie[1,3] <- 0.06
serie[1,4] <- 0.065
serie[1,5] <- 0.07
serie[1,6] <- 0.075
serie[2,1] <- 0.03
serie[2,2] <- 0.04
serie[2,3] <- 0.05
serie[2,4] <- 0.063
serie[2,5] <- 0.074
serie[2,6] <- 0.08
serie[3,1] <- 0.06
serie[3,2] <- 0.065
serie[3,3] <- 0.07
serie[3,4] <- 0.08
serie[3,5] <- 0.084
serie[3,6] <- 0.09
serie[4,1] <- 0.02
serie[4,2] <- 0.03
serie[4,3] <- 0.04
serie[4,4] <- 0.042
serie[4,5] <- 0.045
serie[4,6] <- 0.05
# `market.assets` input
market.assets <- matrix(NA,nrow = 10,ncol = 2)
market.assets[1,1]  <- 0.04
market.assets[2,1]  <- 0.05
market.assets[3,1]  <- 0.06
market.assets[4,1]  <- 0.07
market.assets[5,1]  <- 0.08
market.assets[6,1]  <- 0.09
market.assets[7,1]  <- 0.06
market.assets[8,1]  <- 0.07
market.assets[9,1]  <- 0.075
market.assets[10,1] <- 0.07
market.assets[1,2]  <- "2016-01-01"
market.assets[2,2]  <- "2016-02-01"
market.assets[3,2]  <- "2016-04-01"
market.assets[4,2]  <- "2016-07-01"
market.assets[5,2]  <- "2017-01-01"
market.assets[6,2]  <- "2017-02-01"
market.assets[7,2]  <- "2017-04-01"
market.assets[8,2]  <- "2017-07-01"
market.assets[9,2]  <- "2018-01-01"
market.assets[10,2] <- "2019-01-01"
# Calculation
curve.calculation(serie = serie, market.assets = market.assets, noSpots = 1,  
                  previous.curve = previous.curve, asset.type = "TES",
                  freq = 1, rate.type = 1, fwd = 0,
                  nodes = c(0, 0.25, 0.5, 1:5), approximation = "linear")

## ----results = 'hide'---------------------------------------------------------
# Inputs for calibration of spot curve
yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
nodes <- seq(0,10,0.001)
# Calibration
spot <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
                           analysis.date = "2019-01-03" , asset.type = "IBRSwaps",
                           freq = 4, rate.type = 0, fwd = 0, npieces = NULL, 
                           nodes = nodes, approximation = "linear")
# Spot to Forward
dates <- names(spot)
spot2forward(dates, spot, approximation = "linear")

## ----echo = FALSE, message=FALSE , fig.width = 7.2, fig.height = 5, results = 'hide', fig.align='center',warning=FALSE----
yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
nodes <- seq(0,10,0.001)
# Calibration
y1 <- curve.calibration (market.assets = NULL, yield.curve = yield.curve, asset.type = "IBRSwaps",
                         analysis.date = "2019-01-03" , nodes = nodes,
                         fwd = 0, freq = 4, approximation = "linear", 
                         npieces = NULL, rate.type = 0, obj = "Price")
# Spot to Forward
dates <- as.numeric(names(y1))
y2 <- spot2forward(dates, y1, approximation = "linear")
x  <- as.numeric(names(y1))

df <- data.frame(x = x, y1 = y1, y2 = y2)

# Create the plot with ggplot2
ggplot(df, aes(x)) +
  geom_line(aes(y = y1, color = "Spot"), size = 1.5) +
  geom_line(aes(y = y2, color = "Forward"), size = 1.5) +
  scale_color_manual(values = c("red", "blue"), 
                     name = "", 
                     labels = c("Forward", "Spot")) +
  labs(x = "term", y = "rate", title = "Spot to Forward") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.box = "vertical", 
    legend.box.background = element_rect(color = "black", size = 1.2)
  )

## ----results = 'hide'---------------------------------------------------------
# Inputs for calibration of forward curve
yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
nodes <- seq(0,10,0.001)
# Calibration
fwd <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
                          analysis.date = "2019-01-03", asset.type = "LIBORSwaps",
                          freq = 4, rate.type = 0, daycount = "ACT/365", 
                          npieces = NULL, fwd = 1, nodes = nodes, 
                          approximation = "linear")
# Forward to Spot
dates <- names(fwd)
fwd2spot(dates, fwd, approximation = "linear")

## ----echo = FALSE, message=FALSE , fig.width = 7.2, fig.height = 5, results = 'hide', fig.align='center',warning=FALSE----
yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
nodes <- seq(0,10,0.001)
# Calibration
y2 <- curve.calibration (market.assets = NULL, yield.curve = yield.curve, asset.type = "LIBORSwaps",
                          analysis.date = "2019-01-03" , nodes = nodes,fwd = 1,
                          daycount = "ACT/365", npieces = NULL,rate.type = 0,
                          approximation = "linear", freq = 4 , obj = "Rates")
# Forward to Spot
dates <- as.numeric(names(y2))
y1 <- fwd2spot(dates, y2, approximation = "linear")
x  <- as.numeric(names(y1))

df <- data.frame(x = x, y1 = y1, y2 = y2)

# Create the plot with ggplot2
ggplot(df, aes(x)) +
  geom_line(aes(y = y1, color = "Spot"), size = 1.5) +
  geom_line(aes(y = y2, color = "Forward"), size = 1.5) +
  scale_color_manual(values = c("red", "blue"), 
                     name = "", 
                     labels = c("Forward", "Spot")) +
  labs(x = "term", y = "rate", title = "Forward to Spot") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.box = "vertical", 
    legend.box.background = element_rect(color = "black", size = 1.2)
  )

## -----------------------------------------------------------------------------
ex.rate <- 4814
swaps <- rbind(c("2024-03-01", "FF", 0.07 , 0.0325, NA   , NA    , 2000 * ex.rate, 2000),
               c("2025-03-01", "VV", NA   , NA    , 0.015, 0.0175, 2000 * ex.rate, 2000),
               c("2026-03-01", "FF", 0.075, 0.03  , NA   ,  NA   , 5000000, 5000000 / ex.rate),
               c("2027-03-01", "VV", NA   , NA    , 0.01 , 0.015 , 5000000, 5000000 / ex.rate),
               c("2028-03-01", "FF", 0.08 ,0.035  , NA   , NA    , 3000000, 3000000 / ex.rate),
               c("2029-03-01", "VV", NA   , NA    , 0.01 , 0.0125, 3000000, 3000000 / ex.rate))
colnames(swaps) <- c("Mat"  ,"Legs", "C1" , "C2", "spread1", "spread2", "prin1", "prin2")

## ----results='hide'-----------------------------------------------------------
# Inputs for calibration of spot curve
yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
nodes <- seq(0,10,0.001)
# Calibration of local spot curve
rates <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
                           analysis.date = "2019-01-03" , asset.type = "IBRSwaps", 
                           freq = 4, rate.type = 0, fwd = 0, npieces = NULL, 
                           obj = "Price", nodes = nodes, approximation = "linear")
# Calibration of Basis Curve
nodes <- seq(0,10,0.001)
basis.curve(swaps = swaps, ex.rate = 4814, analysis.date = "2023-03-01", 
            rates = rates, rates2 = rates / 4, freq = c(2,2,2,2,1,1),
            rate.type = 1, npieces = 4, obj = "Price", Weights = NULL, 
            nsimul = 10, nodes = nodes, approximation = "linear")

## ----echo = FALSE, message=FALSE , fig.width = 7.2, fig.height = 5, results = 'hide', fig.align='center', warning = FALSE----
ex.rate <- 4814
swaps <- rbind(c("2024-03-01", "VV", NA   , NA , 0.015  , 0.015    , 2000 * ex.rate, 2000),
               c("2025-03-01", "VV", NA   , NA , 0.0175 , 0.0175, 2000 * ex.rate, 2000),
               c("2026-03-01", "VV", NA   , NA , 0.015   , 0.02   , 5000000, 5000000 / ex.rate),
               c("2027-03-01", "VV", NA   , NA , 0.0125 , 0.0225 , 5000000, 5000000 / ex.rate),
               c("2028-03-01", "VV", NA   , NA , 0.0150 , 0.025    , 3000000, 3000000 / ex.rate),
               c("2029-03-01", "VV", NA   , NA , 0.02   , 0.0275, 3000000, 3000000 / ex.rate))
colnames(swaps) <- c("Mat"  ,"Legs", "C1" , "C2", "spread1", "spread2", "prin1", "prin2")

yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325)
names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,6)
nodes <- seq(0, 10, by = 0.001)
rates <- curve.calibration (market.assets = NULL, yield.curve = yield.curve, asset.type = "IBRSwaps",
                            analysis.date = "2023-03-01" , nodes = nodes,
                            daycount = "ACT/365", fwd = 0, freq = 4,
                            approximation = "constant", npieces = NULL, nsimul = nsimul, 
                            rate.type = 0, obj = "Rates")

nodes <- seq(0, 6, by = 0.001)
y1 <- basis.curve(swaps, analysis.date = "2023-03-01", freq = c(2,2,2,2,1,1),
                  rate.type = 1, ex.rate = 4814, rates = rates, nodes = nodes,
                  rates2 = rates / 4, npieces = NULL,
                  Weights = NULL, nsimul = 1, obj = "Price", approximation = "linear")
y2 <- basis.curve(swaps, analysis.date = "2023-03-01", freq = c(2,2,2,2,1,1),
                  rate.type = 1, ex.rate = 4814, rates = rates, nodes = nodes,
                  rates2 = rates / 4, npieces = NULL,
                  Weights = NULL, nsimul = 1, obj = "Price", approximation = "constant")
x <- as.numeric(names(y1))
df <- data.frame(x = x, y1 = y1, y2= y2)
ggplot(df, aes(x)) +
  geom_line(aes(y = y1, color = "Linear"), linewidth = 1.5) +
  geom_line(aes(y = y2, color = "Constant"), linewidth = 1.5) +
  scale_color_manual(values = c("blue", "red"), name = "") +
  labs(x = "term", y = "rate", title = "Piecewise Linear Basis Curves") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.box = "vertical", 
    legend.box.background = element_rect(color = "black", size = 1.2)
  )


## -----------------------------------------------------------------------------
valuation.swaps(maturity = "2026-07-01", analysis.date = "2023-01-01",
                asset.type = "IBRSwaps", freq = 4, coupon.rate = 0.04, 
                rates = rep(0.04,14), float.rate = 500)

## -----------------------------------------------------------------------------
valuation.swaps(maturity = "2026-07-01", analysis.date = "2023-02-01",
                asset.type = "IBRSwaps", freq = 4, coupon.rate = 0.04,
                rates = rep(0.04,14), float.rate = 0.042)

## ----results = 'hide'---------------------------------------------------------
# Curve Calibration for `rates` input
 yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325)
 names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,6)
 nodes <- seq(0, 10, by = 0.001) # Our curve has nodes with three decimals.
 rates <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
                             analysis.date = "2023-03-01", asset.type = "IBRSwaps",
                             freq = 4, rate.type = 0, daycount = "ACT/365", fwd = 0, 
                             npieces = NULL, obj = "Rates", nsimul = nsimul, 
                             nodes = nodes,  approximation = "constant")

## ----results='hide'-----------------------------------------------------------
# Curve Calibration for `basis.rates` input
nodes  <- seq(0, 10, by = 0.001)
rates2 <- rates/4 # It is assumed foreign curve is proportional to local spot curve.
# Swaps input for calibration
ex.rate <- 4814
swaps <- rbind(c("2024-03-01", "FF", 0.07 , 0.0325, NA   , NA    , 2000 * ex.rate, 2000),
               c("2025-03-01", "VV", NA   , NA    , 0.015, 0.0175, 2000 * ex.rate, 2000),
               c("2026-03-01", "FF", 0.075, 0.03  , NA   ,  NA   , 5000000, 5000000 / ex.rate),
               c("2027-03-01", "VV", NA   , NA    , 0.01 , 0.015 , 5000000, 5000000 / ex.rate),
               c("2028-03-01", "FF", 0.08 ,0.035  , NA   , NA    , 3000000, 3000000 / ex.rate),
               c("2029-03-01", "VV", NA   , NA    , 0.01 , 0.0125, 3000000, 3000000 / ex.rate))
colnames(swaps) <- c("Mat"  ,"Legs", "C1" , "C2", "spread1", "spread2", "prin1", "prin2")
# Calibration
basis.rates <- basis.curve(swaps, ex.rate = 4814, analysis.date = "2023-03-01",
                           rates = rates, rates2 = rates2, freq = c(2,2,2,2,1,1), 
                           rate.type = 1, npieces = NULL, obj = "Price", 
                           Weights = NULL, nodes = nodes, approximation = "linear")

## -----------------------------------------------------------------------------
# Valuation
valuation.swaps (maturity = "2024-03-01", analysis.date = "2023-03-01", asset.type = "CCS",
                 freq = 2, coupon.rate = NA, rates = rates, float.rate = NULL, spread = 0.015,
                 principal = 2000 * ex.rate, Legs = "VV", ex.rate = ex.rate, 
                 basis.rates = basis.rates, coupon.rate2 = NA, rates2 = rates2, 
                 float.rate2 = NULL, spread2 = 0.0175, principal2 = 2000, rate.type = 0, 
                 daycount = "ACT/365", loc = "BOG")

