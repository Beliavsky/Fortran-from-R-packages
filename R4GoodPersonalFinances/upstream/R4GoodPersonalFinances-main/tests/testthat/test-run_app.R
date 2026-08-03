test_that("run_app works", {
  
  output <- run_app(shinylive = TRUE)
  expect_true("shiny.appobj" %in% class(output))
})
