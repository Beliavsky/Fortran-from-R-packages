program saccr_demo
  use saccr
  implicit none

  type(portfolio_result_t) :: result

  call example_ird(result)
  write(*, '(a,f12.2)') "IRD example EAD: ", result%total_ead

  call example_credit(result)
  write(*, '(a,f12.2)') "Credit example EAD: ", result%total_ead

  call example_commodity(result)
  write(*, '(a,f12.2)') "Commodity example EAD: ", result%total_ead

  call example_ird_commodity_margined(result)
  write(*, '(a,f12.2)') "Margined IRD/commodity EAD: ", result%total_ead

  call example_basis_volatility(result)
  write(*, '(a,f12.2)') "Basis/volatility example EAD: ", result%total_ead
  call write_exposure_csv("exposure_per_counterparty.csv", result)
end program saccr_demo
