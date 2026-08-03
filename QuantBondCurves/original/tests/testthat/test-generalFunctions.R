# 1. Coupon.dates ---------------------------------------------------------

test_that("if format of maturity is not valid, an error has to emerge", {
  expect_error(coupon.dates(maturity = "ABC"), "is not valid")
})


test_that("if format of analysis.date is not valid as date an error has to emerge", {
  expect_error(coupon.dates(maturity = "2022-01-01",
                            analysis.date = "2022-01-45"), "is not valid as Date")
})

test_that("if TES, frequency of coupon dates is anual", {
  expect_length(coupon.dates(maturity = "2025-01-05",
                             analysis.date = "2023-01-03",
                             asset.type = "TES")$dates, 3)
})

test_that("if LIBOR, frequency of coupon dates is quarterly", {
  expect_length(coupon.dates(maturity = "2025-01-05",
                             analysis.date = "2023-01-03",
                             asset.type = "LIBOR")$dates, 9)
})

test_that("if LIBORSwaps, frequency of coupon dates is quarterly", {
  expect_length(coupon.dates(maturity = "2025-01-05",
                             analysis.date = "2023-01-03",
                             asset.type = "LIBORSwaps")$dates, 9)
})

test_that("if IBRSwaps, frequency of coupon dates is quarterly", {
  expect_length(coupon.dates(maturity = "2025-01-04",
                             analysis.date = "2023-01-03",
                             asset.type = "IBRSwaps")$dates, 9)
})

test_that("if asset requires the input of frequency, the frequency value is not assumed", {
  expect_length(coupon.dates(maturity = "2025-01-04",
                             analysis.date = "2023-01-03",
                             asset.type = "FixedIncome",
                             freq = 6)$dates, 13)
})

test_that("if maturity as date is given, function assumes coupon payments have already started", {
  expect_length(coupon.dates(maturity = "2025-01-04",
                             analysis.date = "2023-01-03",
                             asset.type = "FixedIncome",
                             freq = 1)$dates, 3)
})

test_that("if maturity as numeric is given, function assumes that trade-date is on analysis.date", {
  res <- as.Date(c("2024-01-05","2025-01-05"))
  expect_equal(coupon.dates(maturity = 2,
                            analysis.date = "2023-01-05",
                            asset.type = "FixedIncome",
                            freq = 1)$dates, res)
})

test_that("if maturity as numeric is given, but trade-date is specified,
          function doesn't assume trade-date is on analysis.date +2 if LIBOR", {
            res <- as.Date(c("2023-01-05","2024-01-05"))
            expect_equal(coupon.dates(maturity = 2,
                                      analysis.date = "2023-01-03",
                                      asset.type = "LIBORSwaps",
                                      freq = 1,
                                      trade.date = "2022-01-03")$dates, res)
          })



test_that("if maturity as numeric is given, but trade-date is specified,
          function doesn't assume trade-date is on analysis.date", {
            res <- as.Date(c("2024-01-03"))
            expect_equal(coupon.dates(maturity = 2,
                                      analysis.date = "2023-01-03",
                                      asset.type = "TES",
                                      freq = 1,
                                      trade.date = "2022-01-03")$dates, res)
          })


test_that("if following business days, then effective dates is one day later after dates,
          when dates lands on a holiday", {
            expect_equal(coupon.dates(maturity = "2023-01-01",
                                      analysis.date = "2022-01-01",
                                      asset.type = "FixedIncome",
                                      freq = 0,
                                      convention = "F")$effective.dates,
                         quantdates::AddBusinessDays(date = "2023-01-01",
                                                     numDate = 1))
          })

test_that("if backward business days, then effective dates is one day before dates,
          when dates lands on a holiday", {
            expect_equal(coupon.dates(maturity = "2023-01-01",
                                      analysis.date = "2022-01-01",
                                      asset.type = "FixedIncome",
                                      freq = 0,
                                      convention = "B")$effective.dates,
                         quantdates::AddBusinessDays(date = "2023-01-01",
                                                     numDate = -1))
          })

test_that("if modified following business days, then effective dates is one day before dates,
          when dates lands on a holiday and when next business days lands on a new month", {
            expect_equal(coupon.dates(maturity = "2023-12-31",
                                      analysis.date = "2022-01-01",
                                      asset.type = "FixedIncome",
                                      freq = 0,
                                      convention = "MF")$effective.dates,
                         quantdates::AddBusinessDays(date = "2023-12-31",
                                                     numDate = -1))
          })

test_that("if modified backward business days, then effective dates is one day after dates,
          when dates lands on a holiday and when previous business days lands on a new month", {
            expect_equal(coupon.dates(maturity = "2023-01-01",
                                      analysis.date = "2022-01-01",
                                      asset.type = "FixedIncome",
                                      freq = 0,
                                      convention = "MB")$effective.dates,
                         quantdates::AddBusinessDays(date = "2023-01-01",
                                                     numDate = 0))
          })

test_that("if 1st of Jan, then following business day has to produce the same outpout as
          modified backward", {
            expect_setequal(coupon.dates(maturity = "2023-01-01",
                                         analysis.date = "2010-01-01",
                                         asset.type = "FixedIncome",
                                         freq = 1,
                                         convention = "F")$effective.dates,
                            coupon.dates(maturity = "2023-01-01",
                                         analysis.date = "2010-01-01",
                                         asset.type = "FixedIncome",
                                         freq = 1,
                                         convention = "MB")$effective.dates)
          })



test_that("short first or long first is assumed if trade-date is not specified and maturity is a date", {
  expect_setequal(coupon.dates(maturity = "2026-05-29",
                               analysis.date = "2023-02-28",
                               asset.type = "FixedIncome",
                               freq = 1, coupon.schedule = "SL")$dates,
                  coupon.dates(maturity = "2026-05-29",
                               analysis.date = "2023-02-28",
                               asset.type = "FixedIncome",
                               freq = 1, coupon.schedule = "LL")$dates)
})

test_that("Short Last works", {
  res <- as.Date(c("2020-08-29","2021-02-28","2021-08-29","2022-02-28","2022-08-29","2023-02-28","2023-08-29","2024-02-29",
                   "2024-08-29","2025-02-28","2025-08-29","2026-02-28","2026-05-29"))
  expect_setequal(coupon.dates(maturity = "2026-05-29",
                               analysis.date = "2020-02-29",
                               asset.type = "FixedIncome",
                               freq = 2, coupon.schedule = "SL",
                               trade.date = "2020-02-29")$dates, res)
})

test_that("if frequency of coupons aligns with maturity and trade.date, then coupon.schedule
          is not taken into account", {
            expect_setequal(coupon.dates(maturity = "2026-02-28",
                                         analysis.date = "2020-02-29",
                                         asset.type = "FixedIncome",
                                         freq = 2, coupon.schedule = "LL",
                                         trade.date = "2020-02-29")$dates,
                            coupon.dates(maturity = "2026-02-28",
                                         analysis.date = "2020-02-29",
                                         asset.type = "FixedIncome",
                                         freq = 2, coupon.schedule = "SF",
                                         trade.date = "2020-02-29")$dates)
          })

test_that("if frequency of coupons aligns with maturity and trade.date, then coupon.schedule
          is not taken into account", {
            expect_setequal(coupon.dates(maturity = 2,
                                         analysis.date = "2020-02-29",
                                         asset.type = "FixedIncome",
                                         freq = 2, coupon.schedule = "LF",
                                         trade.date = "2020-02-29")$dates,
                            coupon.dates(maturity = 2,
                                         analysis.date = "2020-02-29",
                                         asset.type = "FixedIncome",
                                         freq = 2, coupon.schedule = "SL",
                                         trade.date = "2020-02-29")$dates)
          })

test_that("Long First works", {
  res <- as.Date(c("2020-06-29", "2020-09-29", "2020-12-29", "2021-03-29"))
  expect_setequal(coupon.dates(maturity = 1.08,
                               analysis.date = "2020-02-29",
                               asset.type = "FixedIncome",
                               freq = 4, coupon.schedule = "LF",
                               trade.date = "2020-02-29")$dates, res)
})

test_that("Short First works", {
  res <- as.Date(c("2020-09-30","2020-12-31", "2021-03-31", "2021-06-30", "2021-09-30", "2021-12-31",
                   "2022-03-31", "2022-06-30", "2022-09-30"))
  expect_setequal(coupon.dates(maturity = 2.16,
                               analysis.date = "2020-07-31",
                               asset.type = "FixedIncome",
                               freq = 4, coupon.schedule = "SF",
                               trade.date = "2020-07-31")$dates, res)
})



# 2. Coupons --------------------------------------------------------------

test_that("If dates input is given and in a wrong format, an error has to emerge", {
  expect_warning(coupons(dates = "ABC",
                         coupon.rate = 0.04), "formats failed to parse")
})

test_that("If coupon.rate is not a unique number or doesn't has the same length as coupon dates,
          a warning message has to emerge", {
            expect_error(coupons(coupon.rate = c(0.04, 0.05),
                                 maturity = "2026-01-05",
                                 analysis.date = "2023-01-05",
                                 asset.type = "LIBOR",
                                 freq = 4), "length of dates")
          })

test_that("If freq when divided by 12 has a remainder different from zero, then an error
          has to emerge", {
            expect_error(coupons(coupon.rate = 0.04,
                                 maturity = "2026-01-05" ,
                                 asset.type = "LIBOR",
                                 freq = 7), "frequency of payments")
          })

test_that("Output is equal when using a unique coupon.rate or the same coupon rate for each coupon date", {
  expect_equal(coupons(coupon.rate = c(0.04, 0.04,0.04,0.04,0.04),
                       maturity = "2024-01-05" ,
                       analysis.date = "2023-01-03",
                       asset.type = "IBR",
                       freq = 4),
               coupons(coupon.rate = 0.04,
                       maturity = "2024-01-05" ,
                       analysis.date = "2023-01-03",
                       asset.type = "IBR",
                       freq = 4))
})

test_that("Introducing the input of coupon dates also works fine", {
  dates <- coupon.dates(maturity = "2025-01-05", analysis.date = "2023-01-05")
  expect_length(coupons(coupon.rate = 0.04,
                        maturity = "2024-01-05" ,
                        analysis.date = "2023-01-03",
                        asset.type = "IBR",
                        freq = 4,
                        dates = dates$dates),2)
})

test_that("output is equal to Excel calculations", {
  expect_equal(coupons(maturity = "2025-01-05",
                       analysis.date = "2023-01-03",
                       asset.type = "IBR",
                       freq = 4,
                       coupon.rate = 0.04),c(0.01022222, 0.01000000, 0.01011111, 0.01022222, 0.01022222, 0.01011111,
                                             0.01011111,0.01022222,1.01022222))


})

test_that("Coupon calculation is fine when bond is long first, (trade.date is introduced)", {
  expect_equal(coupons(coupon.rate = 0.04,
                       maturity = 1.08 ,
                       analysis.date = "2020-02-29",
                       asset.type = "IBR",
                       freq = 4,
                       trade.date = "2020-02-29",
                       coupon.schedule = "LF"
                       ), c(0.01344444, 0.01022222, 0.01011111, 1.01000000),tolerance = 0.0001)
})

test_that("Coupon calculation is fine when bond is short last, (trade.date is introduced)", {
  expect_equal(coupons(coupon.rate = 0.04,
                       maturity = 1.08 ,
                       analysis.date = "2020-02-29",
                       asset.type = "IBR",
                       freq = 4,
                       trade.date = "2020-02-29",
                       coupon.schedule = "SL"
                       ), c(0.01000000, 0.01022222, 0.01022222, 0.01011111, 1.00322222))
})



# 3. Discount Factors --------------------------------------------------------

test_that("If length of dates is longer than length of rates and length(rates) !=1,
          then an error has to emerge", {
  dates <- coupon.dates(maturity = "2025-01-05", analysis.date = "2015-01-05")
  expect_error(discount.factors(rates = c(0.01,0.02),
                                analysis.date = "2023-01-03",
                                dates = dates$dates), "Missing rates for all coupon")
})

test_that("Output has to have the same length as dates", {
  dates <- coupon.dates(maturity = "2025-01-05")
  expect_length(discount.factors(rates = c(0.01,0.03),
                                 analysis.date = "2023-01-03",
                                 dates = dates$dates), length(dates$dates))
})


test_that("if rates are constant, discount factor has to be each time smaller", {
  dates <- coupon.dates(maturity = "2025-01-05",
                        analysis.date = "2021-01-01",
                        asset.type = "IBR")
  for (i in 2:(length(dates$dates)-1)) {

    expect_gt(discount.factors(rates = rep(0.01,length(dates$dates)),
                               analysis.date = "2021-01-01",
                               dates = dates$dates)[i],
              discount.factors(rates = rep(0.01,length(dates$dates)),
                               analysis.date = "2015-01-01",
                               dates = dates$dates)[i+1])
  }
})



test_that("Discount factors have to be lower when using continously compounded interest rates", {
  dates <- coupon.dates(maturity = "2025-01-05", analysis.date = "2023-01-05")
  for (i in 1:(length(dates$dates))) {
    expect_lte(discount.factors(rates = c(0.01,0.04),
                                analysis.date = "2023-01-03",
                                dates = dates$dates,
                                rate.type = 0)[i],
               discount.factors(rates = c(0.01,0.04),
                                analysis.date = "2023-01-03",
                                dates = dates$dates,
                                rate.type = 1)[i])
  }
})



test_that("outpout is equal to Excel calculations", {
  dates <-  coupon.dates(maturity = "2025-01-05",analysis.date = "2023-01-03", asset.type = "IBR", freq = 4)$effective.dates
  res <- c(0.9999455,0.9962543,0.9901207,0.9815680, 0.9707170,0.9577746,0.9427650,
           0.9254729,0.9065446)
  for (i in seq_along(dates)) {
    expect_equal(discount.factors(rates = c(0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05),
                                  analysis.date = "2023-01-03",
                                  rate.type = 1,
                                  dates = dates)[i], res[i],tolerance = 0.0001)
  }
})

# 4. Discount Time --------------------------------------------------------

test_that("if format of tinitial is not adequate, an error has to emerge", {
  expect_error(discount.time(tinitial = "ABC", tfinal = "2023-01-03"), "ABC is not valid as Date")
})

test_that("if format of tfinal is not adequate, an error has to emerge", {
  expect_error(discount.time(tinitial = "2023-01-01", tfinal = "Hoy"),"Hoy is not valid as Date")
})

test_that("A year has to be a year when intial date and final are the same date,
but with final having, and additional year", {
  expect_equal(discount.time(tinitial = "2023-01-01", tfinal = "2024-01-01"),1)
})

test_that("A year has to be a year, when intial date and final are the same date,
but with final having, and additional year", {
  expect_equal(discount.time(tinitial = "2021-06-13", tfinal = "2023-06-13"),2)
})

test_that("A year has to be a year, when intial date and final are the same date,
but with final having, and additional year", {
  expect_equal(discount.time(tinitial = "2021-12-31", tfinal = "2029-12-31"),8)
})

test_that("When initial date is 29-Feb, then a year is next year on the 28-Feb", {
  expect_equal(discount.time(tinitial = "2024-02-29", tfinal = "2025-02-28"),1)
})

test_that("When initial date is 29-Feb, then next year 27-Feb has to be,
          less than a year.", {
            expect_lt(discount.time(tinitial = "2024-02-29", tfinal = "2025-02-27"),1)
          })

test_that("When initial date is 29-Feb, then next year 27-Feb has to be,
          more than a year.", {
            expect_gt(discount.time(tinitial = "2024-02-29", tfinal = "2025-03-01"),1)
          })

test_that("When initial date is 29-Feb, then four years is next leap date", {
  expect_equal(discount.time(tinitial = "2024-02-29", tfinal = "2028-02-29"),4)
})

# 5. Accrued Interest -----------------------------------------------------

test_that("if analysis date is outside time gap between the previous and upcoming cupon date,
          return is 0.", {
            expect_equal(accrued.interests(analysis.date = "2025-01-03",
                                           coupon.rate = 0.04,
                                           maturity = "2024-01-03",
                                           daycount = "ACT/365",
                                           freq = 1),0)
          })

test_that("if accrued interest is calculated on a coupon date,
          then accrued interest is equal to 0.", {
            expect_equal(accrued.interests(analysis.date = "2023-01-04",
                                           coupon.rate = 0.04,
                                           maturity = "2024-01-04",
                                           daycount = "ACT/365"),0)
          })

test_that("if accrued interest is calculated on a coupon date (with different inputs),
          then accrued interest is equal to 0.", {
            expect_equal(accrued.interests(analysis.date = "2023-01-04",
                                           coupon.rate = 0.04,
                                           maturity = "2024-01-04",
                                           daycount = "ACT/360",
                                           freq = 4,
                                           asset.type = "IBR"),0)
          })

test_that("if analysis date is very close to next coupon date,
          then accrued interest has to be very close to coupon.", {
            expect_gt(accrued.interests(analysis.date = "2023-01-03",
                                        coupon.rate = 0.04,
                                        maturity = "2024-01-04",
                                        daycount = "ACT/365",
                                        freq = 1,
                                        asset.type = "TES"),0.039)
          })

test_that("if analysis date is very close to previous coupon date,
          then accrued interest has to be very close to zero.", {
            expect_lt(accrued.interests(analysis.date = "2023-01-05",
                                        coupon.rate = 0.04,
                                        maturity = "2024-01-04",
                                        daycount = "ACT/365",
                                        freq = 1,
                                        asset.type = "TES"),0.001)
          })
