# simulating multiple scenarios

    Code
      dplyr::summarise(dplyr::group_by(scenarios, scenario_id), total_income = sum(
        total_income))
    Output
      # A tibble: 5 x 2
        scenario_id total_income
              <dbl>        <dbl>
      1          60      2892000
      2          65      3432000
      3          70      3972000
      4          75      4512000
      5          80      5052000

