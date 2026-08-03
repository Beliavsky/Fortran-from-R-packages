test_that("calculating certainty equivalent return (h)", {
  
  h <- calc_certainty_equivalent_return(
    expected_return = 0.03370865,
    variance        = 0.004088745,
    risk_tolerance  = 0.25
  )
  expect_equal(
    round(h, 6),
    0.02529
  )
})

