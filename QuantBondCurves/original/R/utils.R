curve.rates <- function() {
  yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325)
  names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,6)
  nodes <- seq(0, 10, by = 0.001)
  yield.curve <- curve.calibration (market.assets = NULL,
                                    yield.curve = yield.curve,
                                    asset.type = "IBRSwaps",
                                    analysis.date = "2023-03-01" ,
                                    nodes = nodes,
                                    daycount = "ACT/365",
                                    fwd = 0,
                                    freq = 4,
                                    approximation = "constant",
                                    npieces = NULL,
                                    rate.type = 0,
                                    obj = "Rates")
return(yield.curve)
}

swaps.input <- function() {
  ex.rate <- 4814
  swaps <- rbind(c("2024-03-01", "FF", 0.07 , 0.0325, NA   , NA    , 2000 * ex.rate, 2000),
                 c("2025-03-01", "VV", NA   , NA    , 0.015, 0.0175, 2000 * ex.rate, 2000),
                 c("2026-03-01", "FF", 0.075, 0.03  , NA   ,  NA   , 5000000, 5000000 / ex.rate),
                 c("2027-03-01", "VV", NA   , NA    , 0.01 , 0.015 , 5000000, 5000000 / ex.rate),
                 c("2028-03-01", "FF", 0.08 ,0.035  , NA   , NA    , 3000000, 3000000 / ex.rate),
                 c("2029-03-01", "VV", NA   , NA    , 0.01 , 0.0125, 3000000, 3000000 / ex.rate))
  colnames(swaps) <- c("Mat"  ,"Legs", "C1" , "C2", "spread1", "spread2", "prin1", "prin2")
  return(swaps)
}

curve.basis <- function() {
  yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325)
  names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,6)
  nodes <- seq(0, 10, by = 0.001)

  rates <- curve.calibration (market.assets = NULL,
                              yield.curve = yield.curve,
                              asset.type = "IBRSwaps",
                              analysis.date = "2023-03-01" ,
                              nodes = nodes,
                              daycount = "ACT/365",
                              fwd = 0,
                              freq = 4,
                              approximation = "constant",
                              npieces = NULL,
                              rate.type = 0,
                              obj = "Rates")
  rates2 <- rates/4

  basis.rates <- basis.curve(swaps.input(), analysis.date = "2023-03-01", freq = c(2,2,2,2,1,1),
                             rate.type = 1, ex.rate = 4814, rates = rates, nodes = nodes,
                             rates2 = rates2, npieces = NULL,
                             Weights = NULL, daycount = "ACT/365", obj = "Price", approximation = "linear")
  return(basis.rates)
}


market.assets.input <- function() {

market.assets       <- matrix(NA,nrow = 10,ncol = 2)
market.assets[1,1]  <- 0.103
market.assets[2,1]  <- 0.1044
market.assets[3,1]  <- 0.1083
market.assets[4,1]  <- 0.101
market.assets[5,1]  <- 0.112
market.assets[6,1]  <- 0.113
market.assets[7,1]  <- 0.115
market.assets[8,1]  <- 0.116
market.assets[9,1]  <- 0.115
market.assets[10,1] <- 0.13
market.assets[1,2]  <- "2019-01-03"
market.assets[2,2]  <- "2019-02-03"
market.assets[3,2]  <- "2019-04-03"
market.assets[4,2]  <- "2019-07-03"
market.assets[5,2]  <- "2020-01-03"
market.assets[6,2]  <- "2021-01-03"
market.assets[7,2]  <- "2022-01-03"
market.assets[8,2]  <- "2024-07-03"
market.assets[9,2]  <- "2026-01-03"
market.assets[10,2] <-"2029-01-03"

return(market.assets)
}

market.assets.input2 <- function() {
market.assets      <- matrix(NA,nrow = 9,ncol = 2)
market.assets[1,1] <- 0.04
market.assets[2,1] <- 0.06
market.assets[3,1] <- 0.07
market.assets[4,1] <- 0.08
market.assets[5,1] <- 0.09
market.assets[6,1] <- 0.10
market.assets[7,1] <- 0.11
market.assets[8,1] <- 0.12
market.assets[9,1] <- 0.115
market.assets[1,2] <- "2020-01-03"
market.assets[2,2] <- "2021-01-03"
market.assets[3,2] <- "2022-01-03"
market.assets[4,2] <- "2023-01-03"
market.assets[5,2] <- "2024-01-03"
market.assets[6,2] <- "2025-01-03"
market.assets[7,2] <- "2026-01-03"
market.assets[8,2] <- "2027-01-03"
market.assets[9,2] <- "2028-01-03"

return(market.assets)
}

market.assets.input3 <- function() {
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

market.assets[1,1]  <- 0.103
market.assets[2,1]  <- 0.1044
market.assets[3,1]  <- 0.1083
market.assets[4,1]  <- 0.101
market.assets[5,1]  <- 0.112
market.assets[6,1]  <- 0.113
market.assets[7,1]  <- 0.115
market.assets[8,1]  <- 0.116
market.assets[9,1]  <- 0.115
market.assets[10,1] <- 0.13
return(market.assets)
}

serie.input <- function() {
serie<- matrix(NA,nrow = 4,ncol = 6)
rownames(serie) <- c("2014-01-01","2015-01-01","2016-01-01","2017-01-01")
colnames(serie) <- c(0, 0.08333, 0.25, 0.5, 1, 2)
serie[1,1]<- 0.04
serie[1,2]<- 0.05
serie[1,3]<- 0.06
serie[1,4]<- 0.065
serie[1,5]<- 0.07
serie[1,6]<- 0.075
serie[2,1]<- 0.03
serie[2,2]<- 0.04
serie[2,3]<- 0.05
serie[2,4]<- 0.063
serie[2,5]<- 0.074
serie[2,6]<- 0.08
serie[3,1]<- 0.06
serie[3,2]<- 0.065
serie[3,3]<- 0.07
serie[3,4]<- 0.08
serie[3,5]<- 0.084
serie[3,6]<- 0.09
serie[4,1]<- 0.02
serie[4,2]<- 0.03
serie[4,3]<- 0.04
serie[4,4]<- 0.042
serie[4,5]<- 0.045
serie[4,6]<- 0.05

return(serie)
}

market.assets.input4 <- function() {
market.assets       <- matrix(NA,nrow = 10,ncol = 2)
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


market.assets[1,2] <- "2016-01-01"
market.assets[2,2] <- "2016-02-01"
market.assets[3,2] <- "2016-04-01"
market.assets[4,2] <- "2016-07-01"
market.assets[5,2] <- "2017-01-01"
market.assets[6,2] <- "2017-02-01"
market.assets[7,2] <- "2017-04-01"
market.assets[8,2] <- "2017-07-01"
market.assets[9,2] <- "2018-01-01"
market.assets[10,2]<- "2019-01-01"
return(market.assets)
}

spot.input <- function() {
yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
nodes <- seq(0,10,1)

spot <- curve.calibration (market.assets = NULL,
                           yield.curve = yield.curve,
                           asset.type = "IBRSwaps",
                           analysis.date = "2019-01-03" ,
                           nodes = nodes,
                           daycount = "ACT/365",
                           fwd = 0,
                           freq = 4,
                           approximation = "constant",
                           npieces = NULL,
                           rate.type = 0,
                           obj = "Rates")
return(spot)
}

fwd.input <- function() {
  yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
  names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
  nodes <- seq(0,10,2)

  fwd <- curve.calibration (market.assets = NULL,
                            yield.curve = yield.curve,
                            asset.type = "IBRSwaps",
                            analysis.date = "2019-01-03" ,
                            nodes = nodes,
                            daycount = "ACT/365",
                            fwd = 1,
                            freq = 4,
                            approximation = "linear",
                            npieces = NULL,
                            rate.type = 0,
                            obj = "Rates")


  return(fwd)
}
