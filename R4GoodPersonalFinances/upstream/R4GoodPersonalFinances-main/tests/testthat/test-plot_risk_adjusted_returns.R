test_that("plotting risk adjusted returns for optimal allocation", {

  skip_on_cran()
  skip_on_ci()

  safe_asset_return         <- 0.015
  risky_asset_return_mean   <- 0.040 + safe_asset_return
  risky_asset_return_sd     <- 0.200
  risk_aversion             <- 2

  rar_ra2 <- 
    plot_risk_adjusted_returns(
      safe_asset_return         = safe_asset_return,
      risky_asset_return_mean   = risky_asset_return_mean,
      risky_asset_return_sd     = risky_asset_return_sd,
      risk_aversion             = risk_aversion
    )
  if (interactive()) print(rar_ra2)
  
   vdiffr::expect_doppelganger("rar_ra2", 
                                rar_ra2)
      
  rar_ra3 <- 
    plot_risk_adjusted_returns(
      safe_asset_return         = safe_asset_return,
      risky_asset_return_mean   = risky_asset_return_mean,
      risky_asset_return_sd     = risky_asset_return_sd,
      risk_aversion             = risk_aversion + 1
    )
  if (interactive()) print(rar_ra3)
  
  vdiffr::expect_doppelganger("rar_ra3", 
                              rar_ra3)
})

test_that("plotting risk adjusted returnswith current allocation", {

  skip_on_cran()
  skip_on_ci()

  safe_asset_return         <- 0.015
  risky_asset_return_mean   <- 0.040 + safe_asset_return
  risky_asset_return_sd     <- 0.200
  risk_aversion             <- 2 + 1

   rar_craa_1 <- 
    plot_risk_adjusted_returns(
      current_risky_asset_allocation = 0.6666,
      safe_asset_return              = safe_asset_return,
      risky_asset_return_mean        = risky_asset_return_mean,
      risky_asset_return_sd          = risky_asset_return_sd,
      risk_aversion                  = risk_aversion 
    )
  if (interactive()) {print(rar_craa_1)}

  vdiffr::expect_doppelganger("rar_craa_1",
                               rar_craa_1)

   rar_craa_2 <- 
    plot_risk_adjusted_returns(
      current_risky_asset_allocation = 0.3333,
      safe_asset_return              = safe_asset_return,
      risky_asset_return_mean        = risky_asset_return_mean,
      risky_asset_return_sd          = risky_asset_return_sd,
      risk_aversion                  = risk_aversion
    )
  if (interactive()) {print(rar_craa_2)}

  vdiffr::expect_doppelganger("rar_craa_2",
                               rar_craa_2)

   rar_craa_3 <- 
    plot_risk_adjusted_returns(
      current_risky_asset_allocation = 0,
      safe_asset_return              = safe_asset_return,
      risky_asset_return_mean        = risky_asset_return_mean,
      risky_asset_return_sd          = risky_asset_return_sd,
      risk_aversion                  = risk_aversion
    )
  if (interactive()) {print(rar_craa_3)}

  vdiffr::expect_doppelganger("rar_craa_3",
                               rar_craa_3)

})
