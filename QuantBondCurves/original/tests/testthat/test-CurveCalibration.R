# 13. Curve Calibration Spot --------------------------------------------------

test_that("output has to have always the same length as nodes input vector", {
  yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
  nodes <- c(0.117,0.15,0.25,0.35,2,3,6,8)
  expect_length(curve.calibration(market.assets = market.assets.input2(),
                                  yield.curve = yield.curve,
                                  asset.type = "TES",
                                  analysis.date = "2019-01-03" ,
                                  nodes = nodes,
                                  daycount = "ACT/365",
                                  rate.type = 1,
                                  fwd = 0,
                                  freq = 1,
                                  approximation = "constant") , length(nodes))
})


test_that("output has to have always the same length as nodes input vector", {
  yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
  nodes <- c(0.117)
  expect_length(curve.calibration(market.assets = market.assets.input2(),
                                  yield.curve = yield.curve,
                                  asset.type = "TES",
                                  analysis.date = "2019-01-03" ,
                                  nodes = nodes,
                                  fwd = 0,
                                  daycount = "ACT/365",
                                  rate.type = 1,
                                  freq = 1,
                                  approximation = "linear"), length(nodes))
})

test_that("if node date is after last maturty of bonds, then it is assumed that the rate associated
          is the same as last final node", {
            yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
            nodes <- c(0.1,1,2,3,10,100)
            expect_equal(as.numeric(curve.calibration(market.assets = market.assets.input2(),
                                                      yield.curve = yield.curve,
                                                      asset.type = "LIBORSwaps",
                                                      analysis.date = "2019-01-03" ,
                                                      nodes = nodes,
                                                      daycount = "ACT/365",
                                                      rate.type = 1,
                                                      freq = 4,
                                                      fwd = 0,
                                                      approximation = "linear")[5]),
                         as.numeric(curve.calibration(market.assets = market.assets.input2(),
                                                      yield.curve = yield.curve,
                                                      asset.type = "LIBORSwaps",
                                                      analysis.date = "2019-01-03" ,
                                                      nodes = nodes,
                                                      daycount = "ACT/365",
                                                      rate.type = 1,
                                                      freq = 4,
                                                      approximation = "linear")[6]))
          })

test_that("constant option throws a constant curve", {
  yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
  nodes <- seq(0.255, 0.5, by = 0.025)
  expect_setequal(curve.calibration(market.assets = market.assets.input2(),
                                    yield.curve = yield.curve,
                                    asset.type = "IBRSwaps",
                                    analysis.date = "2019-01-03" ,
                                    nodes = nodes,
                                    daycount = "ACT/365",
                                    rate.type = 1,
                                    freq = 4,
                                    approximation = "constant"
                                    ), rep(yield.curve[1],length(nodes)))
})


test_that("linear option throws a linear curve", {
  yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
  nodes <- seq(7.25, 8, by = 0.25)
  market.assets <- market.assets.input2()
  for (i in seq_along(nodes)) {
    expect_gte(curve.calibration(market.assets = market.assets,
                                 yield.curve = yield.curve,
                                 asset.type = "IBRSwaps",
                                 analysis.date = "2015-01-03" ,
                                 nodes = nodes,
                                 daycount = "ACT/365",
                                 rate.type = 0,
                                 freq = 4,
                                 approximation = "linear"
                                 )[i], yield.curve[3])
  }
})

test_that("linear option throws a linear curve", {
  yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
  nodes <- seq(7.25, 8, by = 0.25)
  market.assets <- market.assets.input2()
  for (i in seq_along(nodes)) {
    expect_lte(curve.calibration(market.assets = market.assets,
                                 yield.curve = yield.curve,
                                 asset.type = "IBRSwaps",
                                 analysis.date = "2015-01-03" ,
                                 nodes = nodes,
                                 daycount = "ACT/365",
                                 rate.type = 0,
                                 freq = 4,
                                 approximation = "linear"
                                 )[i], yield.curve[4])
  }
})

test_that("continously compounded rates have a shorter output rate than anually compounded", {

  yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
  nodes <- c(0.117,0.15,0.25,0.35,2,3,6,8)
  market.assets <- market.assets.input2()
  for (i in 1:length(nodes)) {
    expect_lte(curve.calibration(market.assets = market.assets,
                                 yield.curve = yield.curve,
                                 asset.type = "IBRSwaps",
                                 analysis.date = "2018-01-03" ,
                                 nodes = nodes,
                                 daycount = "ACT/365",
                                 rate.type = 0,
                                 freq = 4,
                                 approximation = "constant")[i],
               curve.calibration(market.assets = market.assets,
                                 yield.curve = yield.curve,
                                 asset.type = "IBRSwaps",
                                 analysis.date = "2018-01-03" ,
                                 nodes = nodes,
                                 daycount = "ACT/365",
                                 rate.type = 1,
                                 freq = 4,
                                 approximation = "constant")[i])
  }
})

test_that("asset.type isnt that relevant for methodology of calibration, can be determined with
          the other parameters", {

  yield.curve <- c(0.103,0.1044,0.1083,0.101,0.112,0.113,0.115,0.116,0.115,0.13)
  names(yield.curve) <- c("1","2","3","4","5","6","7","8","9","10")
  nodes <- c(0.117,0.15,0.25,0.35,2,3,6,8)
  for (i in 1:length(nodes)) {
    expect_lte(curve.calibration(market.assets = NULL,
                                 yield.curve = yield.curve,
                                 asset.type = "IBRSwaps",
                                 analysis.date = "2018-01-03" ,
                                 nodes = nodes,
                                 daycount = "ACT/365",
                                 rate.type = 0,
                                 freq = 4,
                                 approximation = "linear")[i],
               curve.calibration(market.assets = NULL,
                                 yield.curve = yield.curve,
                                 asset.type = "FixedIncome",
                                 analysis.date = "2018-01-03" ,
                                 nodes = nodes,
                                 daycount = "ACT/365",
                                 rate.type = 0,
                                 freq = 4,
                                 noSpots = 4,
                                 approximation = "linear")[i])
  }
})



test_that("freq can be a vector, for bonds", {
  yield.curve <- c(0.04,0.06,0.07,0.08,0.09,0.10, 0.11,0.12,0.125)
  nodes <- c(0.117,0.15,0.25,0.35,2,3,6,8)
  expect_length(curve.calibration(market.assets = market.assets.input2(),
                                  yield.curve = yield.curve,
                                  asset.type = "FixedIncome",
                                  analysis.date = "2019-01-03",
                                  nodes = nodes,
                                  daycount = "ACT/365",
                                  noSpots = 1,
                                  rate.type = 1,
                                  fwd = 0,
                                  freq = c(1,2,4,2,1,1,2,4),
                                  approximation = "constant"), length(nodes))
})

# 13 b. Curve Calibration Fwd -----------------------------------------------
test_that("output has to have always the same length as nodes input vector", {
  yield.curve <- c(0.103,0.1044,0.1083,0.101,0.112,0.113,0.115,0.116,0.115,0.13)
  names(yield.curve) <- c("1","2","3","4","5","6","7","8","9","10")
  nodes <- seq(0,10, by = 0.001)
  expect_length(curve.calibration(market.assets = market.assets.input3(),
                                  yield.curve = yield.curve,
                                  asset.type = "TES",
                                  analysis.date = "2019-01-03" ,
                                  nodes = nodes,
                                  daycount = "ACT/365",
                                  freq = 1,
                                  fwd = 1,
                                  approximation = "linear"), length(nodes))
})

test_that("output has to have always the same length as nodes input vector", {
  yield.curve <- c(0.103,0.1044,0.1083,0.101,0.112,0.113,0.115,0.116,0.115,0.13)
  names(yield.curve) <- c("1","2","3","4","5","6","7","8","9","10")
  nodes <- c(0.117)
  expect_length(curve.calibration (market.assets = market.assets.input3(),
                                   yield.curve = yield.curve,
                                   asset.type = "TES",
                                   analysis.date = "2019-01-03" ,
                                   nodes = nodes,
                                   daycount = "ACT/365",
                                   freq = 4,
                                   fwd = 1), length(nodes))
})

test_that("if node date is after last maturty of bonds, then it is assumed that the rate associated
          is the same as last final node", {
            nodes <- c(0.1,1,2,5,7,8,10,100)
            yield.curve <- cbind(0.02,0.023,0.025,0.03,0.035,0.036,0.039,0.42,0.04)
            colnames(yield.curve) <- cbind("0","0.08","0.25","0.5","5","6","7","8","9")

            market.assets <- market.assets.input()[-10,]
            expect_equal(as.numeric(curve.calibration (market.assets = market.assets,
                                                       yield.curve = yield.curve,
                                                       asset.type = "LIBORSwaps",
                                                       analysis.date = "2019-01-03" ,
                                                       nodes = nodes,
                                                       daycount = "ACT/365",
                                                       fwd = 1,
                                                       freq = 4)[7]),
                         as.numeric(curve.calibration(market.assets = market.assets,
                                                      yield.curve = yield.curve,
                                                      asset.type = "LIBORSwaps",
                                                      analysis.date = "2019-01-03" ,
                                                      nodes = nodes,
                                                      daycount = "ACT/365",
                                                      freq = 4,
                                                      fwd = 1)[8]))
          })

#Piecewise linear spot vs forward
test_that("plot of curves is correct", {
  skip("Test is for plotting curves")
  skip_on_cran()
  yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325, 0.1320)
  names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,7,10)

  for (i in 1:9) {
    if (i %in% c(1:4) ) {
      xf <- curve.calibration (market.assets = NULL,
                               yield.curve = yield.curve,
                               asset.type = "IBRSwaps",
                               analysis.date = "2019-01-03" ,
                               nodes = nodes,
                               daycount = "ACT/365",
                               fwd = 1,
                               freq = 4,
                               approximation = "linear",
                               npieces = i,
                               nsimul = nsimul,
                               rate.type = 0,
                               obj = "Rates")
      y <- as.numeric(names(xf))
      plot(y,xf,type = "p", col = "purple", xlab = "Tasa", ylab = "t", ylim = c(0,0.20))
      legend("bottomright", legend = c("Forward"), col = c("purple"),lty=1:1,cex = 0.8)

    } else {
      x <- curve.calibration (market.assets = market.assets,
                              yield.curve = yield.curve,
                              asset.type = "IBRSwaps",
                              analysis.date = "2019-01-03" ,
                              nodes = nodes,
                              daycount = "ACT/365",
                              fwd = 0,
                              freq = 4,
                              approximation = "linear",
                              npieces = i - 4,
                              nsimul = nsimul,
                              rate.type = 0,
                              obj = "Rates")

      y <- as.numeric(names(x))

      xf <- curve.calibration (market.assets = market.assets,
                               yield.curve = yield.curve,
                               asset.type = "IBRSwaps",
                               analysis.date = "2019-01-03" ,
                               nodes = nodes,
                               daycount = "ACT/365",
                               fwd = 1,
                               freq = 4,
                               approximation = "linear",
                               npieces = i,
                               nsimul = nsimul,
                               rate.type = 0,
                               obj = "Rates")

      plot(y,x,type = "p", col = "blue", xlab = "t", ylab = "Tasa", ylim = c(0,0.20))
      points(y,xf, col = "purple")
      legend("bottomright", legend = c("Spot","Forward"), col = c("blue","purple"),lty=1:2,cex = 0.8)

    }
  }

  #Bootstrapping Method

  xf <- curve.calibration (market.assets = market.assets,
                           yield.curve = yield.curve,
                           asset.type = "IBRSwaps",
                           analysis.date = "2019-01-03",
                           nodes = nodes,
                           daycount = "ACT/365",
                           fwd = 1,
                           freq = 4,
                           approximation = "linear",
                           npieces = NULL,
                           rate.type = 0,
                           obj = "Price")
  y <- as.numeric(names(xf))

  x <- curve.calibration (market.assets = market.assets,
                          yield.curve = yield.curve,
                          asset.type = "IBRSwaps",
                          analysis.date = "2019-01-03" ,
                          nodes = nodes,
                          daycount = "ACT/365",
                          fwd = 0,
                          freq = 4,
                          approximation = "linear",
                          npieces = NULL,
                          rate.type = 0,
                          obj = "Rates")

  plot(y,x,type = "p", col = "blue", xlab = "t", ylab = "Tasa", ylim = c(0,0.20))
  points(y,xf, col = "purple")
  legend("topleft", legend = c("Spot","Forward"), col = c("blue","purple"),lty=1:2,cex = 0.8)
})

#Piecewise constant spot vs forward

test_that("plot of constant curves is correct", {
  skip("Test is for plotting curves")
  skip_on_cran()

  for (i in 1:9) {

    if (i %in% c(1:4) ) {
      xf <- curve.calibration (market.assets = market.assets,
                               yield.curve = yield.curve,
                               asset.type = "IBRSwaps",
                               analysis.date = "2019-01-03" ,
                               nodes = nodes,
                               daycount = "ACT/365",
                               fwd = 1,
                               freq = 4,
                               approximation = "constant",
                               npieces = i,
                               nsimul = nsimul,
                               rate.type = 0,
                               obj = "Rates")
      y <- as.numeric(names(xf))
      plot(y,xf,type = "p", col = "purple", xlab = "Tasa", ylab = "t",ylim = c(0,0.20))
      legend("bottomright", legend = c("Forward"), col = c("purple"),lty=1:1,cex = 0.8)

    } else {
      x <- curve.calibration (market.assets = market.assets,
                              yield.curve = yield.curve,
                              asset.type = "IBRSwaps",
                              analysis.date = "2019-01-03" ,
                              nodes = nodes,
                              daycount = "ACT/365",
                              fwd = 0,
                              freq = 4,
                              approximation = "constant",
                              npieces = i - 4,
                              nsimul = nsimul,
                              rate.type = 0,
                              obj = "Rates")

      y <- as.numeric(names(x))

      xf <- curve.calibration (market.assets = market.assets,
                               yield.curve = yield.curve,
                               asset.type = "IBRSwaps",
                               analysis.date = "2019-01-03" ,
                               nodes = nodes,
                               daycount = "ACT/365",
                               fwd = 1,
                               freq = 4,
                               approximation = "constant",
                               npieces = i,
                               nsimul = nsimul,
                               rate.type = 0,
                               obj = "Rates")

      plot(y,x,type = "p", col = "blue", xlab = "t", ylab = "Tasa", ylim = c(0,0.20))
      points(y,xf, col = "purple")
      legend("bottomright", legend = c("Spot","Forward"), col = c("blue","purple"),lty=1:2,cex = 0.8)
    }
  }

  #Bootstrapping Method

  xf <- curve.calibration (market.assets = market.assets,
                           yield.curve = yield.curve,
                           asset.type = "IBRSwaps",
                           analysis.date = "2019-01-03" ,
                           nodes = nodes,
                           daycount = "ACT/365",
                           fwd = 1,
                           freq = 4,
                           approximation = "constant",
                           npieces = NULL,
                           rate.type = 0,
                           obj = "Rates")
  y <- as.numeric(names(xf))

  x <- curve.calibration (market.assets = market.assets,
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

  plot(y,x,type = "p", col = "blue", xlab = "t", ylab = "Tasa", ylim = c(0,0.15))
  points(y,xf, col = "purple")
  legend("bottomright", legend = c("Spot","Forward"), col = c("blue","purple"),lty=1:2,cex = 0.8)
})

# 14. Curve Calculation ---------------------------------------------------
test_that("historical rates that begin before the previous curve, must be eliminated", {
  previous.curve <- matrix(0.04,nrow = 2,ncol = 28)
  rownames(previous.curve) <- c("2014-01-01","2015-01-01")
  colnames(previous.curve) <- c(0, 0.25, 0.5, 1:25)
  expect_equal(as.numeric(curve.calculation(previous.curve = previous.curve,
                                            market.assets = market.assets.input4(),
                                            serie = serie.input(),
                                            asset.type = "TES",
                                            rate.type = 1,
                                            freq = 1,
                                            fwd = 0,
                                            daycount = NULL,
                                            approximation = "linear",
                                            nodes = c(0, 0.25, 0.5, 1:25))[1,]), rep(0.04,28) )
})

test_that("historical rates that begin before the previous curve, must be eliminated", {
  previous.curve <- matrix(0.04,nrow = 2,ncol = 28)
  rownames(previous.curve) <- c("2014-01-01","2015-01-01")
  colnames(previous.curve) <- c(0, 0.25, 0.5, 1:25)
  expect_equal(as.numeric(curve.calculation(previous.curve = previous.curve,
                                            market.assets = market.assets.input4(),
                                            serie = serie.input(),
                                            asset.type = "TES",
                                            noSpots = 1,
                                            rate.type = 1,
                                            freq = 1,
                                            fwd = 0,
                                            daycount = NULL,
                                            approximation = "linear",
                                            nodes = c(0, 0.25, 0.5, 1:25))[2,]), rep(0.04,28) )
})


test_that("historical rates that begin before the previous curve, must be eliminated,
          when constructing an instantaneous forward curve", {
            previous.curve <- matrix(0.04,nrow = 2,ncol = 28)
            rownames(previous.curve) <- c("2014-01-01","2015-01-01")
            colnames(previous.curve) <- c(0, 0.25, 0.5, 1:25)
            expect_equal(as.numeric(curve.calculation(previous.curve = previous.curve,
                                                      market.assets = market.assets.input4(),
                                                      serie = serie.input(),
                                                      asset.type = "IBRSwaps",
                                                      noSpots = 4,
                                                      rate.type = 0,
                                                      freq = 4,
                                                      fwd = 1,
                                                      daycount = NULL,
                                                      approximation = "constant",
                                                      nodes = c(0, 0.25, 0.5, 1:25))[1,]), rep(0.04,28) )
          })
