test_that("calculate optimal allocation", {

  expect_equal(
    calc_optimal_risky_asset_allocation(risky_asset_return_mean = 0.05,
                                        risky_asset_return_sd   = 0.20,
                                        safe_asset_return       = 0,
                                        risk_aversion           = 2),
    0.625)
})

test_that("calculate optimal allocation with zeros for all inputs", {

  expect_equal(
    calc_optimal_risky_asset_allocation(risky_asset_return_mean = 0,
                                        risky_asset_return_sd   = 0,
                                        safe_asset_return       = 0,
                                        risk_aversion           = 0),
    0)
})

test_that("calculate optimal allocation with zeros for some inputs", {

  expect_equal(
    calc_optimal_risky_asset_allocation(risky_asset_return_mean = 0,
                                        risky_asset_return_sd   = 0,
                                        safe_asset_return       = 0,
                                        risk_aversion             = 2),
    0)
})

test_that("calculate optimal allocation with zeros and vector inputs", {

  expect_equal(
    calc_optimal_risky_asset_allocation(risky_asset_return_mean = 0,
                                        risky_asset_return_sd   = 0,
                                        safe_asset_return       = 0,
                                        risk_aversion             =
                                          c(2, 2, 2)),
    c(0, 0, 0))
})

