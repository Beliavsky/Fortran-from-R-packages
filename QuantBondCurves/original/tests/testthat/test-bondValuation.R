# 6. Valuation.bonds ---------------------------------------------------------

test_that("if asset is a TES and coupon.rate has a length longer than 1,
          then an error has to emerge", {
            expect_error(valuation.bonds(maturity = "2026-01-05",
                                         coupon.rate = c(0.04,0.05),
                                         rates = 0.06,
                                         analysis.date = "2023-01-03",
                                         asset.type = "TES"), "unique fixed rate")
          })

test_that("if asset is FixedIncome and coupon.rate has a length longer than 1,
          then an error has to emerge", {
            expect_error(valuation.bonds(maturity = "2026-01-05",
                                         coupon.rate = c(0.02,0.03),
                                         rates = 0.04,
                                         analysis.date = "2023-01-03",
                                         asset.type = "FixedIncome"), "unique fixed rate")
          })

test_that("if length of coupon payments is longer than length of rates and method iz Zero,
          then an error has to emerge", {
            expect_error(valuation.bonds(maturity = "2026-01-05",
                                         coupon.rate = 0.02,
                                         rates = c(0.04, 0.05),
                                         analysis.date = "2023-01-03",
                                         asset.type = "FixedIncome",
                                         freq = 1), "Missing rates for all")
          })

test_that("if length of coupon payments is longer than length of rates and method iz Zero,
          then an error has to emerge", {
            expect_error(valuation.bonds(maturity = "2026-01-05",
                                         coupon.rate = 0.02,
                                         rates = c(0.04, 0.05),
                                         analysis.date = "2023-01-03",
                                         asset.type = "FixedIncome",
                                         freq = 1), "Missing rates for all")
          })

test_that("if there a spread for TES asset, then a warning has to emerge", {
  expect_warning(valuation.bonds(maturity = "2026-01-05",
                                 coupon.rate = 0.02,
                                 rates = 0.04,
                                 analysis.date = "2023-01-03",
                                 asset.type = "TES",
                                 freq = 1,
                                 spread = 0.1), "Spread should be")
})

test_that("bond valuation is equal either inserting a unique coupon.rate
          or a vector of the same coupon.rate for every coupon date.", {
            expect_equal(valuation.bonds(maturity = "2026-01-05",
                                         coupon.rate = 0.04,
                                         rates = 0.04,
                                         analysis.date = "2023-01-03",
                                         asset.type = "IBR",
                                         freq = 4,
                                         spread = 0.1),
                         valuation.bonds(maturity = "2026-01-05",
                                         coupon.rate = rep(0.04,13),
                                         rates = 0.04,
                                         analysis.date = "2023-01-03",
                                         asset.type = "IBR",
                                         freq = 4,
                                         spread = 0.1))
          })

test_that("bond valuation is lower when using continous compounded.rates vs using
          anually compounded", {
            expect_gte(valuation.bonds(maturity = "2026-01-05",
                                       coupon.rate = 0.04,
                                       rates = 0.04,
                                       analysis.date = "2023-01-03",
                                       asset.type = "IBR",
                                       freq = 4,
                                       spread = 0.05,
                                       rate.type = 1),
                       valuation.bonds(maturity = "2026-01-05",
                                       coupon.rate = rep(0.04,13),
                                       rates = 0.04,
                                       analysis.date = "2023-01-03",
                                       asset.type = "IBR",
                                       freq = 4,
                                       spread = 0.05,
                                       rate.type = 0))
          })


test_that("bond valuation is equal either using a unique rate or a vector with the same unique rate", {
  expect_equal(valuation.bonds(maturity = "2026-01-05",
                               coupon.rate = 0.04,
                               rates = rep(0.04,13),
                               analysis.date = "2023-01-03",
                               asset.type = "IBR",
                               freq = 4,
                               spread = 0.01,
                               rate.type = 1),
               valuation.bonds(maturity = "2026-01-05",
                               coupon.rate = rep(0.04,13),
                               rates = 0.04,
                               analysis.date = "2023-01-03",
                               asset.type = "IBR",
                               freq = 4,
                               spread = 0.01,
                               rate.type = 1))
})

test_that("bond valuation is lower or equal when calculating clean price vs dirty price", {
  expect_lte(valuation.bonds(maturity = "2026-01-05",
                             coupon.rate = 0.04,
                             rates = rep(0.06,13),
                             analysis.date = "2023-01-03",
                             asset.type = "IBR",
                             freq = 4,
                             spread = 0.01,
                             rate.type = 0,
                             dirty = 0),
             valuation.bonds(maturity = "2026-01-05",
                             coupon.rate = rep(0.04,13),
                             rates = 0.06,
                             analysis.date = "2023-01-03",
                             asset.type = "IBR",
                             freq = 4,
                             spread = 0.01,
                             rate.type = 0,
                             dirty = 1))
})

test_that("including spread on coupon.rate is equivalent to calculate it apart each one,
          like this function does", {
            expect_lte(valuation.bonds(maturity = "2026-01-05",
                                       coupon.rate = 0.04,
                                       rates = rep(0.04,13),
                                       analysis.date = "2023-01-03",
                                       asset.type = "IBR",
                                       freq = 4,
                                       spread = 0),
                       valuation.bonds(maturity = "2026-01-05",
                                       coupon.rate = rep(0.02,13),
                                       rates = 0.04,
                                       analysis.date = "2023-01-03",
                                       asset.type = "IBR",
                                       freq = 4,
                                       spread = 0.02,
                                       rate.type = 1,
                                       dirty = 1))
          })

test_that("if coupon.rates are equal to rates, then bond is valued par", {
  expect_equal(valuation.bonds(maturity = "2024-01-05",
                               coupon.rate = 0.04,
                               rates = rep(0.04,4),
                               analysis.date = "2021-01-05",
                               asset.type = "TES",
                               freq = 1,
                               spread = 0,
                               principal = 1), 1)
})

test_that("Output is exactly as Excels calculations", {
  expect_equal(valuation.bonds(maturity = "2026-01-05",
                               coupon.rate = 0.15,
                               rates = 0.147,
                               analysis.date = "2016-01-05",
                               asset.type = "TES",
                               freq = 1,
                               principal = 1000), 1015.02, tolerance = 0.001)
})

test_that("Output is exactly as Excels calculations", {
  expect_equal(valuation.bonds(maturity = "2024-02-29",
                               coupon.rate = 0.1256,
                               rates = seq(0.05,0.14,by = 0.005),
                               analysis.date = "2019-07-14",
                               asset.type = "FixedIncome",
                               freq = 4,
                               principal = 567,
                               daycount = "ACT/360",
                               rate.type = 0),  565.70268)
})


test_that("LF and LL can be included in bonds valuation", {
  expect_lt(valuation.bonds(maturity = "2024-02-29",
                            coupon.rate = 0.1256,
                            rates = seq(0.05,0.14,by = 0.005),
                            analysis.date = "2019-07-14",
                            asset.type = "FixedIncome",
                            freq = 4,
                            principal = 567,
                            daycount = "ACT/360",
                            rate.type = 0,
                            trade.date = "2019-07-14",
                            coupon.schedule = "LL"),
            valuation.bonds(maturity = "2024-02-29",
                            coupon.rate = 0.1256,
                            rates = seq(0.05,0.14,by = 0.005),
                            analysis.date = "2019-07-14",
                            asset.type = "FixedIncome",
                            freq = 4,
                            principal = 567,
                            daycount = "ACT/360",
                            rate.type = 0,
                            trade.date = "2019-07-14",
                            coupon.schedule = "LF"))
})

test_that("LF and LL can be included in bonds valuation and maturity as numeric", {
  expect_lt(valuation.bonds(maturity = 4.58,
                            coupon.rate = 0.1256,
                            rates = seq(0.05,0.14,by = 0.005),
                            analysis.date = "2019-07-14",
                            asset.type = "FixedIncome",
                            freq = 4,
                            principal = 567,
                            daycount = "ACT/360",
                            rate.type = 0,
                            trade.date = "2019-07-14",
                            coupon.schedule = "LL"),
            valuation.bonds(maturity = 4.58,
                            coupon.rate = 0.1256,
                            rates = seq(0.05,0.14,by = 0.005),
                            analysis.date = "2019-07-14",
                            asset.type = "FixedIncome",
                            freq = 4,
                            principal = 567,
                            daycount = "ACT/360",
                            rate.type = 0,
                            trade.date = "2019-07-14",
                            coupon.schedule = "LF"))
})

# 7. Swap Valuation -------------------------------------------------------

test_that("if format of analysis.date is not adequate, an error has to emerge", {
  expect_error(valuation.swaps(maturity = "2026-06-01",
                               coupon.rate = 0.04,
                               rates = rep(0.04,4),
                               analysis.date = "2023,01,02"), "is not valid as Date")
})

test_that("if length of cash flows is longer than length of rates, then an error has to emerge", {
  expect_error(valuation.swaps(maturity = "2026-06-01",
                               coupon.rate = 0.04,
                               rates = rep(0.04,2),
                               analysis.date = "2023-01-02"), "Missing rates")
})

test_that("if analysis.date is inside coupon dates, then fixed leg is valued par,
          and float.rate input is useless", {
            expect_equal(valuation.swaps(maturity = "2026-07-01",
                                         coupon.rate = 0.04,
                                         rates = rep(0.04,14),
                                         analysis.date = "2023-01-01",
                                         freq = 4,
                                         asset.type = "IBRSwaps",
                                         float.rate = 500), 0, tolerance = 0.005)
          })

test_that("if analysis.date is outside coupon dates, then fixed leg is not valued par,
          and float.rate input has effect", {
            expect_gte(valuation.swaps(maturity = "2026-07-01",
                                       coupon.rate = 0.04,
                                       rates = rep(0.04,14),
                                       analysis.date = "2023-01-02",
                                       freq = 4,
                                       asset.type = "IBRSwaps",
                                       float.rate = 100), 0)
          })

test_that("including spread on floating note,
is the same as substracting the same spread on coupon.rate,
          and leaving spread on zero", {
            expect_equal(valuation.swaps(maturity = "2026-07-01",
                                         coupon.rate = 0.04,
                                         rates = rep(0.05,14),
                                         analysis.date = "2023-01-02",
                                         freq = 4,
                                         asset.type = "IBRSwaps",
                                         float.rate = 0.03,
                                         spread = 0),
                         valuation.swaps(maturity = "2026-07-01",
                                         coupon.rate = 0.06,
                                         rates = rep(0.05,14),
                                         analysis.date = "2023-01-02",
                                         freq = 4,
                                         asset.type = "IBRSwaps",
                                         float.rate = 0.03,
                                         spread = 0.02))
          })

test_that("value of a swap with par fixed leg, is closer to par when analysis date is on coupon dates", {
  expect_lt(abs(valuation.swaps(maturity = "2026-07-01",
                                coupon.rate = 0.04,
                                rates = rep(0.04,13),
                                analysis.date = "2023-04-01",
                                freq = 4,
                                asset.type = "IBRSwaps",
                                float.rate = 0.04)),
            abs(valuation.swaps(maturity = "2026-07-01",
                                coupon.rate = 0.04,
                                rates = rep(0.04,13),
                                analysis.date = "2023-04-02",
                                freq = 4,
                                asset.type = "IBRSwaps",
                                float.rate = 0.04)))
})

test_that("output is equivalent to Excel calculations", {
  expect_equal(valuation.swaps(maturity = "2028-02-29",
                               coupon.rate = 0.06,
                               rates = seq(0.06,0.09,by = 0.0025),
                               analysis.date = "2025-01-02",
                               freq = 4,
                               asset.type = "IBRSwaps",
                               float.rate = 0.03,
                               spread = 0,
                               daycount = "ACT/365",
                               principal = 985,
                               rate.type = 0),73.007317)
})

test_that("output is equivalent to Excel calculations", {
  ex.rate <- 4814
  nodes <- seq(0, 10, by = 0.001)
  expect_equal(valuation.swaps (maturity = "2027-03-01", principal = 50000000000,
                                coupon.rate = NULL, float.rate = NULL, spread = 0.01,
                                rate.type = 0, asset.type = "CCS", rates = curve.rates(),
                                float.rate2 = NULL, analysis.date = "2023-03-01",
                                daycount = "ACT/365", freq = 2, loc = "BOG",
                                convention = "F", Legs = "VV", coupon.rate2 = NULL,
                                spread2 = 0.015, rates2 = curve.rates()/4, ex.rate = 4814,
                                basis.rates = curve.basis(), principal2 =  50000000000/ ex.rate), 67420756)
})
test_that("output is equivalent to Excel calculations", {
  ex.rate <- 4814
  nodes <- seq(0, 10, by = 0.001)
  expect_equal(valuation.swaps (maturity = "2024-03-01", principal = 2000 * ex.rate,
                                coupon.rate = 0.07, float.rate = NULL, rates = curve.rates(),
                                rate.type = 0, asset.type = "CCS", analysis.date = "2023-03-01",
                                float.rate2 = NULL,
                                daycount = "ACT/365", freq = 2, loc = "BOG", spread = 0,
                                convention = "F", Legs = "FF", coupon.rate2 = 0.0325,
                                spread2 = 0, rates2 = curve.rates()/4, ex.rate = 4814,
                                basis.rates = curve.basis(), principal2 = 2000),-34590.32, tolerance = 0.0001)
})

test_that("output is equivalent to Excel calculations", {
  ex.rate <- 4814
  nodes <- seq(0, 10, by = 0.001)
  expect_equal(valuation.swaps (maturity = "2029-03-01", principal = 30000000000 ,
                                coupon.rate = NA, float.rate = NULL, rates = curve.rates(),
                                rate.type = 0, asset.type = "CCS", analysis.date = "2023-03-01",
                                float.rate2 = NULL,
                                daycount = "ACT/365", freq = 1, loc = "BOG", spread = 0.01,
                                convention = "F", Legs = "VV", coupon.rate2 = 0.0325,
                                spread2 = 0.0125, rates2 = curve.rates()/4, ex.rate = 4814,
                                basis.rates = curve.basis(), principal2 = 30000000000/ ex.rate),62985905, tolerance = 0.0001)
})


# 8. price.dirty2clean -------------------------------------------------------

test_that("if analysis.date is after maturity, then an error has to emerge", {
  expect_error(price.dirty2clean(analysis.date = "2025-01-03",
                                 maturity = "2024-01-03",
                                 price = 1,
                                 coupon.rate = 0.04), "after maturity")
})

test_that("if dirty 2 clean, then clean has to be lower or equal than initial price input", {
  expect_lte(price.dirty2clean(analysis.date = "2023-01-02",
                               maturity = "2026-01-03",
                               price = 1,
                               coupon.rate = 0.04,
                               freq = 4,
                               dirty = 1), 1)
})

test_that("if clean 2 dirty, then dirty has to be bigger or equal than initial price input", {
  expect_gte(price.dirty2clean(analysis.date = "2023-01-02",
                               maturity = "2026-01-03",
                               price = 1,
                               coupon.rate = 0.04,
                               freq = 4,
                               dirty = 0), 1)
})

test_that("if analysis date in coupon dates, then dirty price is equal to clean.", {
  expect_equal(price.dirty2clean(analysis.date = "2023-01-03",
                                 maturity = "2026-01-03",
                                 price = 1,
                                 coupon.rate = 0.04,
                                 freq = 4,
                                 dirty = 0), 1)
})


# 9. bond.price2rate ---------------------------------------------------------


test_that("if TES and price par, then rate has to be coupon.rate ", {
  expect_equal(bond.price2rate(coupon.rate = 0.04,
                               analysis.date = "2021-01-03",
                               maturity = "2023-01-03",
                               freq = 1,
                               asset.type = "TES",
                               daycount = "ACT/365",
                               price = 1), 0.04, tolerance = 0.0001)
})

test_that("if price almost par, then rate has to be close to coupon.rate ", {
  expect_equal(bond.price2rate(coupon.rate = 0.04,
                               analysis.date = "2021-01-03",
                               maturity = "2026-01-03",
                               freq = 4,
                               asset.type = "IBR",
                               daycount = "ACT/365",
                               price = 1), 0.04, tolerance = 0.1)
})

test_that("if price is higher, then rate has to be lower", {
  expect_lte(bond.price2rate(coupon.rate = 0.04,
                             analysis.date = "2021-01-03",
                             maturity = "2023-01-03",
                             freq = 4,
                             asset.type = "IBR",
                             daycount = "ACT/365",
                             price = 1.03),
             bond.price2rate(coupon.rate = 0.04,
                             analysis.date = "2021-01-03",
                             maturity = "2023-01-03",
                             freq = 4,
                             asset.type = "IBR",
                             daycount = "ACT/365",
                             price = 0.98)
  )
})

test_that("bond price 2 rate is the opposite of valuation of bonds with the same discount rate
          output of bond price 2 rate", {
            expect_equal(valuation.bonds(
              rates         = bond.price2rate(coupon.rate = 0.04,
                                              analysis.date = "2019-04-12",
                                              maturity = "2023-01-03",
                                              freq = 4,
                                              asset.type = "IBR",
                                              daycount = "ACT/365",
                                              price = 0.95),
              maturity      = "2023-01-03",
              principal     = 1,
              coupon.rate   = 0.04,
              rate.type     = 1,
              asset.type    = "IBR",
              daycount      = "ACT/365",
              analysis.date = "2019-04-12",
              freq          = 4,
              spread        = 0,
              dirty         = 1), 0.95, tolerance = 0.005)
          })

test_that("rate is equivalent to Excel Calculations", {
  expect_equal(bond.price2rate(coupon.rate = 0.1256,
                               analysis.date = "2019-07-14",
                               maturity = "2024-02-29",
                               freq = 4,
                               asset.type = "IBR",
                               daycount = "ACT/360",
                               price = 755.67626,
                               principal = 567,
                               rate.type = 0), 0.05, tolerance = 0.005)
})



# 10. Bonds Sensibility ------------------------------------------------------

test_that("sensibility must be the same for the same bond, either by giving as input a rate or bond price", {
  expect_equal(sens.bonds(input = c("price"),
                          maturity = "2023-01-03",
                          analysis.date = "2019-01-05",
                          coupon.rate = 0.04,
                          dirty = 1,
                          price = 0.98,
                          asset.type = "IBR",
                          daycount = "ACT/365",
                          principal = 1,
                          rate.type = 1,
                          freq = 4),
               sens.bonds(input = "rate",
                          maturity = "2023-01-03",
                          analysis.date = "2019-01-05",
                          coupon.rate = 0.04,
                          dirty = 1,
                          price = bond.price2rate(coupon.rate = 0.04,
                                                  analysis.date = "2019-01-05",
                                                  maturity = "2023-01-03",
                                                  freq = 4,
                                                  asset.type = "IBR",
                                                  daycount = "ACT/365",
                                                  price = 0.98),
                          asset.type = "IBR",
                          daycount = "ACT/365",
                          principal = 1,
                          rate.type = 1,
                          freq = 4), tolerance = 0.005)
})


test_that("the more coupons a bond has, the larger effective sensibility is", {
  expect_gt(sens.bonds(input = c("price"),
                       maturity = "2023-01-03",
                       analysis.date = "2015-01-03",
                       coupon.rate = 0.04,
                       dirty = 1,
                       price = 1.02,
                       asset.type = "FixedIncome",
                       daycount = "ACT/365",
                       principal = 1,
                       rate.type = 1,
                       freq = 1),
            sens.bonds(input = c("price"),
                       maturity = "2023-01-03",
                       analysis.date = "2020-01-03",
                       coupon.rate = 0.04,
                       dirty = 1,
                       price = 1.02,
                       asset.type = "FixedIncome",
                       daycount = "ACT/365",
                       principal = 1,
                       rate.type = 1,
                       freq = 1))
})

test_that("the larger the frequence of coupons, the smaller sensibility is", {
  expect_lt(sens.bonds(input = c("price"),
                       maturity = "2023-01-03",
                       analysis.date = "2020-01-03",
                       coupon.rate = 0.04,
                       dirty = 1,
                       price = 1.02,
                       asset.type = "FixedIncome",
                       daycount = "ACT/365",
                       principal = 1,
                       rate.type = 1,
                       freq = 4),
            sens.bonds(input = c("price"),
                       maturity = "2023-01-03",
                       analysis.date = "2020-01-03",
                       coupon.rate = 0.04,
                       dirty = 1,
                       price = 1.02,
                       asset.type = "FixedIncome",
                       daycount = "ACT/365",
                       principal = 1,
                       rate.type = 1,
                       freq = 1))
})

test_that("sensibility can be calculated when input of price is a rate vector", {
  expect_equal(sens.bonds(input = c("rate"),
                          maturity = "2023-01-03",
                          analysis.date = "2015-01-03",
                          coupon.rate = 0.04,
                          dirty = 1,
                          price = rep(0.08,8),
                          asset.type = "FixedIncome",
                          daycount = "ACT/365",
                          principal = 1,
                          rate.type = 1,
                          freq = 1),
               sens.bonds(input = c("rate"),
                          maturity = "2023-01-03",
                          analysis.date = "2015-01-03",
                          coupon.rate = 0.04,
                          dirty = 1,
                          price = 0.08,
                          asset.type = "FixedIncome",
                          daycount = "ACT/365",
                          principal = 1,
                          rate.type = 1,
                          freq = 1))
})

test_that("sensibility is equivalent to Excel Calculations", {
  expect_equal(sens.bonds(input = "rate",
                          coupon.rate = 0.1256,
                          analysis.date = "2019-07-14",
                          maturity = "2024-02-29",
                          freq = 4,
                          asset.type = "IBR",
                          daycount = "ACT/360",
                          price = seq(0.05, 0.14, by = 0.005),
                          principal = 567,
                          rate.type = 0), 3.4360538)
})


test_that("sensibility is equivalent to Excel Calculations", {
  expect_equal(sens.bonds(input = "rate",
                          coupon.rate = 0.1256,
                          analysis.date = "2019-07-14",
                          maturity = "2024-02-29",
                          freq = 4,
                          asset.type = "IBR",
                          daycount = "ACT/360",
                          price = rep(0.05,19),
                          principal = 567,
                          rate.type = 0), 3.68073711, tolerance = 0.00001)
})

test_that("sensibility is equivalent to Excel Calculations", {
  expect_equal(sens.bonds(input = "price",
                          coupon.rate = 0.1256,
                          analysis.date = "2019-07-14",
                          maturity = "2024-02-29",
                          freq = 4,
                          asset.type = "IBR",
                          daycount = "ACT/360",
                          price = 755.67626,
                          principal = 567,
                          rate.type = 0), 3.68073711, tolerance = 0.00001)
})

# 11. Average Weighted Life -----------------------------------------------


test_that("AVWLife for a zero coupon bond has to be equal to maturity in years.", {
  expect_equal(average.life(input = c("price"),
                            maturity = "2024-01-05",
                            analysis.date = "2022-01-05",
                            coupon.rate = 0.04,
                            dirty = 1,
                            price = 1,
                            asset.type = "FixedIncome",
                            daycount = "ACT/365",
                            principal = 1,
                            rate.type = 0,
                            freq = 0.5),2)
})

test_that("AVWLife is equal to sensibility when rates are continously compounded.", {
  expect_equal(average.life(input = c("price"),
                            maturity = "2023-01-03",
                            analysis.date = "2021-01-03",
                            coupon.rate = 0.04,
                            dirty = 1,
                            price = 1,
                            asset.type = "FixedIncome",
                            daycount = "ACT/365",
                            principal = 1,
                            rate.type = 0,
                            freq = 1),
               sens.bonds(input = c("price"),
                          maturity = "2023-01-03",
                          analysis.date = "2021-01-03",
                          coupon.rate = 0.04,
                          dirty = 1,
                          price = 1,
                          asset.type = "FixedIncome",
                          daycount = "ACT/365",
                          principal = 1,
                          rate.type = 0,
                          freq = 1),tolerance = 0.01)
})

test_that("sensibility is equivalent to Excel Calculations", {
  expect_equal(average.life(input = "price",
                            coupon.rate = 0.1256,
                            analysis.date = "2019-07-14",
                            maturity = "2024-02-29",
                            freq = 4,
                            asset.type = "IBR",
                            daycount = "ACT/360",
                            price = 755.67626,
                            principal = 567,
                            rate.type = 0), 3.7337, tolerance = 0.0001)
})


# 12. Basis Curve Calibration --------------------------------------------------
test_that("basis curve is equivalent to Excel calculations", {
  ex.rate <- 4814
  nodes <- seq(0, 10, by = 1)

  xf <- basis.curve(swaps.input(), analysis.date = "2023-03-01", freq = c(2,2,2,2,1,1),
                    rate.type = 1, ex.rate = 4814, rates = curve.rates(), nodes = nodes,
                    rates2 = curve.rates() / 4, npieces = NULL,
                    Weights = NULL, daycount = "ACT/365", obj = "Price", approximation = "linear")

  res <- c(0.08062099, 0.08062099, 0.03471961, 0.07871210,0.03767375, 0.08072442,
           0.03521362, 0.03496157, 0.03496157, 0.03496157, 0.03496157)
  names(res) <- seq(0, 10, by = 1)

  expect_equal(xf, res,tolerance = 0.001)
})

test_that("basis curve plot is correct", {
  skip("Test is to construct the plot of basis curve")
  skip_on_cran()
  nodes <- seq(0, 10, by = 0.001)

  xf <- basis.curve(swaps.input(), analysis.date = "2023-03-01", freq = c(2,2,2,2,1,1),
                    rate.type = 1, ex.rate = 4814, rates = curve.rates(), nodes = nodes,
                    rates2 = curve.rates() / 4, npieces = NULL,
                    Weights = NULL, nsimul = 1, obj = "Price", approximation = "linear")

  x  <- basis.curve(swaps.input(), analysis.date = "2023-03-01", freq = c(2,2,2,2,1,1),
                    rate.type = 1, ex.rate = 4814, rates = curve.rates(), nodes = nodes,
                    rates2 = curve.rates() / 4, npieces = NULL,
                    Weights = NULL, nsimul = 1, obj = "Price", approximation = "constant")

  y <- names(xf)

  plot(y,x,type = "p", col = "blue", xlab = "t", ylab = "Basis", ylim = c(0,0.20))
  points(y,xf, col = "purple")
  legend("topright", legend = c("Constant","Linear"), col = c("blue","purple"),lty=1:2,cex = 0.8)
})
