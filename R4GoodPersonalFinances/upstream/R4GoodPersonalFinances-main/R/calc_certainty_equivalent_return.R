calc_certainty_equivalent_return <- function(
  expected_return,
  variance,
  risk_tolerance
) {

  (1 + expected_return) * exp(-variance / (2 * risk_tolerance)) - 1
}
