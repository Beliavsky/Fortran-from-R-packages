# printing portfolio

    Code
      portfolio
    Message
      
      -- Portfolio -------------------------------------------------------------------
      
      -- Market assumptions --
      
      -- Expected real returns: 
    Output
      # A tibble: 2 x 3
        name                    expected_return standard_deviation
        <chr>                             <dbl>              <dbl>
      1 GlobalStocksIndexFund            0.0461               0.15
      2 InflationProtectedBonds          0.02                 0   
    Message
      
      -- Correlation matrix: 
    Output
                              GlobalStocksIndexFund InflationProtectedBonds
      GlobalStocksIndexFund                       1                       0
      InflationProtectedBonds                     0                       1
    Message
      
      -- Weights --
      
    Output
                              human_capital liabilities
      GlobalStocksIndexFund             0.5         0.5
      InflationProtectedBonds           0.5         0.5
      
    Message
      -- Accounts --
      
    Output
      # A tibble: 3 x 4
        name                    taxable taxadvantaged total
        <chr>                     <dbl>         <dbl> <dbl>
      1 GlobalStocksIndexFund         0             0     0
      2 InflationProtectedBonds       0             0     0
      3 total                         0             0     0
      
    Message
      -- Pre-tax --
      
    Output
      # A tibble: 2 x 7
        name                    turnover income_qualified capital_gains_long_term
        <chr>                      <dbl>            <dbl>                   <dbl>
      1 GlobalStocksIndexFund       0.04                0                       1
      2 InflationProtectedBonds     0.04                0                       1
        income capital_gains cost_basis
         <dbl>         <dbl>      <dbl>
      1      0        0.0461      0.637
      2      0        0.02        0.820
      
    Message
      -- After-tax --
      
      ! After-tax information is not available yet.
      i Use `calc_effective_tax_rate()` to calculate it.
    Output
      
    Message
      
      -- Allocation --
      
      ! Allocation information is not available yet.
      i Use `calc_optimal_asset_allocation()` to calculate it.

