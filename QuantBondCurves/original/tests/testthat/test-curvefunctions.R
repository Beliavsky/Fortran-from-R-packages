# 15. Spot2fwd ------------------------------------------------------------
test_that("an error has to emerge is length of dates and spot doesn't match.", {
  dates <- names(spot.input())
  spot2 <- spot.input()[-2]
  expect_error(spot2forward(dates, spot2, approximation = "linear"),"must have the same length")
})


test_that("forward rate is always greater or equal than spot rate", {
  dates <- names(spot.input())
  spot  <- spot.input()
  for (i in seq_along(spot)) {
    expect_gte(spot2forward(dates, spot, approximation = "linear")[i],
               spot[i])
  }
})

test_that("if input length is 1 and the maturity of that input
            is different from zero, then output has to have length 2", {
              spot <- spot.input()[2]
              dates <- names(spot.input())
              dates <- dates[2]
              expect_length(spot2forward(dates, spot, approximation = "linear"),2)
            })

test_that("if spot rate is selected to be constant, then a warning message appears
          and output is the same as spot rates.", {
            dates <- names(spot.input())
            expect_warning(spot2forward(dates, spot.input(), approximation = "constant"),"Zero spot constant")
          })


#Plot spot 2 forward linear
test_that("forward rate", {
  skip("Test is for plotting the forward curve")
  skip_on_cran()
  nodes <- seq(0, 10, by = 0.001)
  yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
  names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
  x <- curve.calibration (market.assets = NULL,
                          yield.curve = yield.curve,
                          asset.type = "IBRSwaps",
                          analysis.date = "2019-01-03" ,
                          nodes = nodes,
                          daycount = "ACT/365",
                          fwd = 0,
                          freq = 4,
                          approximation = "linear",
                          npieces = 7,
                          rate.type = 0,
                          obj = "Price")

  y <- names(x)

  xf <- spot2forward(y, x, approximation = "linear")

  plot(y,x,type = "p", col = "blue", xlab = "t", ylab = "Tasa", ylim = c(0,0.20))
  points(y,xf, col = "purple")
  legend("topright", legend = c("Spot","Forward"), col = c("blue","purple"),lty=1:2,cex = 0.8)
})
# 16. fwd2spot ------------------------------------------------------------
test_that("an error has to emerge is length of dates and forward doesn't match.", {
  dates <- names(fwd.input())
  fwd2 <- fwd.input()[-2]
  expect_error(fwd2spot(dates, fwd2, approximation = "linear"),"must have the same length")
})

test_that("if input length is 1 and the maturity of that input
            is different from zero, then output has to have length 2", {
              dates <- names(fwd.input())
              fwd2 <- fwd.input()[2]
              dates2 <- dates[2]
              expect_length(fwd2spot(dates2, fwd2, approximation = "linear"),2)
            })

test_that("spot curve is always lower or equal than forward rate", {
  skip_on_cran()
  fwd <- fwd.input()
  dates <- names(fwd)
  for (i in seq_along(fwd)) {
    expect_lte(fwd2spot(dates, fwd, approximation = "linear")[i],
               fwd[i])
  }
})


test_that("spot curve is always lower or equal than forward rate", {
  skip_on_cran()
  dates <- names(fwd.input())
  fwd <- fwd.input()
  x <- fwd2spot(dates, fwd, approximation = "constant")
  for (i in seq_along(fwd)) {
    expect_lte(x[i],fwd[i])
  }
})

test_that("forward piecewise constant has a larger area under the curve,
          than linear piecewise constant", {
            fwd <- fwd.input()
            dates <- names(fwd)
            for (i in seq_along(fwd)) {
              expect_gte(fwd2spot(dates, fwd, approximation = "constant")[i],
                         fwd2spot(dates, fwd, approximation = "linear")[i])
            }
          })

#Plot forward 2 spot linear
test_that("spot linear plot is correct", {
  skip("Test is for plotting the spot linear curve")
  skip_on_cran()
  nodes <- seq(0, 10, by = 0.001)
  x <- curve.calibration (market.assets = NULL,
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

  y <- names(x)

  xf <- fwd2spot(y, x, approximation = "linear")

  plot(y,x,type = "p", col = "purple", xlab = "t", ylab = "Tasa", ylim = c(0,0.20))
  points(y,xf, col = "blue")
  legend("topright", legend = c("Forward","Spot"), col = c("purple","blue"),lty=1:2,cex = 0.8)
})


#Plot forward 2 spot constant
test_that("spot constant plot is correct", {
  skip("Test is for plotting the spot constant curve")
  skip_on_cran()
  x <- curve.calibration(market.assets = NULL,
                         yield.curve = yield.curve,
                         asset.type = "IBRSwaps",
                         analysis.date = "2019-01-03",
                         nodes = nodes,
                         daycount = "ACT/365",
                         fwd = 1,
                         freq = 4,
                         approximation = "constant",
                         npieces = NULL,
                         rate.type = 0,
                         obj = "Rates")

  y <- names(x)

  xf <- fwd2spot(y, x, approximation = "constant")


  plot(y,x,type = "p", col = "purple", xlab = "t", ylab = "Tasa", ylim = c(0,0.20))
  points(y,xf, col = "blue")
  legend("topright", legend = c("Forward","Spot"), col = c("purple","blue"),lty=1:2,cex = 0.8)
})
