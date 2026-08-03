test_that("decreasing purchasing power due to inflation", {

  skip_on_cran()
  skip_on_ci()

  plot <- plot_purchasing_power(10, 
                                real_interest_rate = -0.02, 
                                legend_title = "Inflation rate",
                                seed = 1234)
  if (interactive()) print(plot)
  
  vdiffr::expect_doppelganger("pp_inflation", plot)
})  

test_that("multiple intrest rates", {
  
  skip_on_cran()
  skip_on_ci()

  plot <- 
    plot_purchasing_power(10, 
                          real_interest_rate = seq(from = 0.04, to = -0.04, by = -0.02),
                          seed = 1234)
  if (interactive()) print(plot)

  vdiffr::expect_doppelganger("pp_multiple", plot)
})

test_that("multiple intrest rates with overlapping labels", {
  
  skip_on_cran()
  skip_on_ci()

  plot <- 
    plot_purchasing_power(10,
                          real_interest_rate = seq(from = 0.04, to = -0.04, by = -0.01),
                          seed = 1234)
  if (interactive()) print(plot)

  vdiffr::expect_doppelganger("pp_overlapping", plot)
})

test_that("multiple intrest rates for 10 years only", {
  
  skip_on_cran()
  skip_on_ci()

  plot <- 
    plot_purchasing_power(10, 
                          years = 10,
                          real_interest_rate = seq(from = 0.04, to = -0.04, by = -0.01),
                          seed = 1234)
  if (interactive()) print(plot)

  vdiffr::expect_doppelganger("pp_10years", plot)
})

