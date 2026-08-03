test_that("plotting retirement ruin", {

  skip_on_cran()
  skip_on_ci()

  rr1 <- function() plot_retirement_ruin(
    portfolio_return_mean = 0.03,
    portfolio_return_sd   = 0.15,
    portfolio_value       = 1000000,
    age                   = 65,
    gompertz_mode         = 84,
    gompertz_dispersion   = 10
  ) 
  if (interactive()) print(rr1())
  vdiffr::expect_doppelganger("rr1", rr1)

  rr2 <- function() plot_retirement_ruin(
    monthly_spendings = 0.04 * 1000000 / 12,
    portfolio_return_mean = 0.03,
    portfolio_return_sd   = 0.15,
    portfolio_value       = 1000000,
    age                   = 65,
    gompertz_mode         = 84,
    gompertz_dispersion   = 10
  ) 
  if (interactive()) print(rr2())
  vdiffr::expect_doppelganger("rr2", rr2)

  rr3 <- function() plot_retirement_ruin(
    monthly_spendings = 0.07 * 1000000 / 12,
    portfolio_return_mean = 0.03,
    portfolio_return_sd   = 0.15,
    portfolio_value       = 1000000,
    age                   = 65,
    gompertz_mode         = 84,
    gompertz_dispersion   = 10
  ) 
  if (interactive()) print(rr3())
  vdiffr::expect_doppelganger("rr3", rr3)

})  

