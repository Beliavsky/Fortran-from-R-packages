! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_tax_aware_portfolio
  use r4good_personal_finances
  implicit none
  type(tax_assumptions) :: assumptions
  type(tax_result) :: taxes
  type(portfolio_result) :: result
  real(dp) :: mu(9), sd(9), corr(9,9), effective(9), hc(9), liab_w(9)
  integer :: i

  mu = [0.0468_dp,0.0501_dp,0.0505_dp,0.0540_dp,0.0269_dp,0.0288_dp,0.0190_dp,0.0329_dp,0.0250_dp]
  sd = [0.1542_dp,0.1795_dp,0.1671_dp,0.2142_dp,0.0379_dp,0.0581_dp,0.03138274_dp,0.0833_dp,0.0055_dp]
  corr = 0.0_dp
  do i = 1, 9
    corr(i,i) = 1.0_dp
  end do
  allocate(assumptions%turnover(9), assumptions%income_qualified(9), assumptions%capital_gains_long_term(9), &
    assumptions%income(9), assumptions%capital_gains(9), assumptions%cost_basis(9))
  assumptions%capital_gains = [0.0349_dp,0.0387_dp,0.0336_dp,0.0388_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
  assumptions%income = mu - assumptions%capital_gains
  assumptions%turnover = [0.3300_dp,0.3652_dp,0.1800_dp,0.3300_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp]
  assumptions%cost_basis = [0.9364_dp,0.9393_dp,0.8750_dp,0.9301_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp]
  assumptions%income_qualified = [0.9762_dp,0.9032_dp,0.7998_dp,0.7387_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
  assumptions%capital_gains_long_term = [0.9502_dp,0.9032_dp,0.8951_dp,0.9023_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
  call effective_tax_rates(mu, sd, assumptions, 0.20_dp, 0.40_dp, taxes)
  call assert_vector_close(taxes%effective_tax_rate, [0.2571344_dp,0.2652787_dp,0.2682702_dp,0.2705776_dp, &
    0.4_dp,0.4_dp,0.4_dp,0.4_dp,0.4_dp], 2.0e-6_dp, "effective tax rates")
  effective = taxes%effective_tax_rate
  effective(7) = 0.0_dp
  call optimize_portfolio(0.35_dp, mu, sd, corr, result, effective_tax_rates_input=effective, &
    fraction_taxable=250000.0_dp/270500.0_dp)
  call assert_vector_close(result%taxable, &
       [0.2712545_dp, 0.2482917_dp, 0.2924056_dp, 0.1122625_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
    5.0e-5_dp, "taxable allocation")
  call assert_vector_close(result%taxadvantaged, [0.0_dp,0.0_dp,0.0_dp,0.0757856_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp], &
    5.0e-5_dp, "tax advantaged allocation")

  hc = [0.0936_dp,0.0468_dp,0.0468_dp,0.0_dp,0.3746_dp,0.3746_dp,0.0_dp,0.0_dp,0.0635_dp]
  liab_w = [0.1263_dp,0.0_dp,0.0_dp,0.0_dp,0.3368_dp,0.3789_dp,0.0_dp,0.0_dp,0.1581_dp]
  mu(7) = -0.10_dp
  call optimize_portfolio(0.35_dp, mu, sd, corr, result, effective_tax_rates_input=effective, &
    fraction_taxable=250000.0_dp/270500.0_dp, financial_wealth=270500.0_dp, human_capital=2767689.0_dp, &
    liabilities=1392064.0_dp, nondiscretionary_consumption=40000.0_dp, discretionary_consumption=46000.0_dp, &
    income=75000.0_dp, life_insurance_premium=(1.0_dp-0.99999_dp)/(1.0_dp+0.025_dp)*1200000.0_dp, &
    human_capital_weights=hc, liabilities_weights=liab_w)
  call assert_vector_close(result%taxable, &
       [0.1322951_dp, 0.0210542_dp, 0.1512539_dp, 0.6196112_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
    8.0e-5_dp, "net worth taxable allocation")
  print *, "test_tax_aware_portfolio: PASS"
contains
  subroutine assert_vector_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    if (size(actual) /= size(expected) .or. maxval(abs(actual - expected)) > tolerance) then
      print *, "FAIL: ", trim(message)
      print *, actual
      print *, expected
      error stop 1
    end if
  end subroutine assert_vector_close
end program test_tax_aware_portfolio
