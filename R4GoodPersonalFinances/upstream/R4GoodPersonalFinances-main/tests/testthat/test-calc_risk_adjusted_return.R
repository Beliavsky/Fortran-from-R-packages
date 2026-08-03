test_that("calculating risk adjusted return for optimal allocation", {
  
  safe_asset_return         <- 0.015
  risky_asset_return_mean   <- 0.04 + safe_asset_return
  risky_asset_return_sd     <- 0.20
  risk_aversion             <- 2
  
  optimal_allocation <-
    calc_optimal_risky_asset_allocation(
      risky_asset_return_mean = risky_asset_return_mean,
      risky_asset_return_sd   = risky_asset_return_sd,
      safe_asset_return       = safe_asset_return,
      risk_aversion           = risk_aversion)
  
  expect_equal(optimal_allocation, 0.5)

  risk_adjusted_return <- 
    calc_risk_adjusted_return(
      risky_asset_return_mean = risky_asset_return_mean,
      safe_asset_return       = safe_asset_return,
      risky_asset_allocation  = optimal_allocation,
      risky_asset_return_sd   = risky_asset_return_sd,
      risk_aversion           = risk_aversion
    )
  
  expect_equal(risk_adjusted_return, 0.025)
  
  expect_equal(
    calc_risk_adjusted_return(
      safe_asset_return       = safe_asset_return,
      risky_asset_return_mean = risky_asset_return_mean,
      risky_asset_allocation  = optimal_allocation
    ),
    risk_adjusted_return
  )
})

test_that("calculating risk adjusted return for allocation vector", {

  safe_asset_return         <- 0.015
  risky_asset_return_mean   <- 0.04 + safe_asset_return
  risky_asset_return_sd     <- 0.20
  risk_aversion             <- 2

  allocation_vector <- seq(from = 0.00, to = 1.00, by = 0.01)

  plot_ra2 <-
    function() calc_risk_adjusted_return(
      risky_asset_return_mean = risky_asset_return_mean,
      safe_asset_return       = safe_asset_return,
      risky_asset_allocation  = allocation_vector,
      risky_asset_return_sd   = risky_asset_return_sd,
      risk_aversion           = 2) |> 
      plot(x = allocation_vector)
  if (interactive()) print(plot_ra2())
    
  vdiffr::expect_doppelganger("plot_ra2", plot_ra2)

  plot_ra3 <-
    function() calc_risk_adjusted_return(
      risky_asset_return_mean = risky_asset_return_mean,
      safe_asset_return       = safe_asset_return,
      risky_asset_allocation  = allocation_vector,
      risky_asset_return_sd   = risky_asset_return_sd,
      risk_aversion           = 3) |> 
      plot(x = allocation_vector)
  if (interactive()) print(plot_ra3())
  
  vdiffr::expect_doppelganger("plot_ra3", plot_ra3)
})

test_that("both parameters cannot be NULL", {

  calc_risk_adjusted_return(
    safe_asset_return = 0.02,
    risky_asset_return_mean = 0.04,
    risky_asset_allocation = 0.5,
    risk_aversion = 2
  ) |> expect_error("cannot be NULL")
    
  calc_risk_adjusted_return(
    safe_asset_return = 0.02,
    risky_asset_return_sd = 0.20,
    risky_asset_allocation = 0.5
  ) |> expect_error("cannot be NULL")
})
