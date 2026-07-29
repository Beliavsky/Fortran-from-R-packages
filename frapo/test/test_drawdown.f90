! SPDX-License-Identifier: GPL-3.0-or-later
program test_drawdown
  use frapo
  implicit none

  real(dp) :: prices(8, 3)
  type(portfolio_result) :: maxdd, avedd, cdar, mincdar, soft

  prices(:, 1) = [100.0_dp, 102.0_dp, 101.0_dp, 104.0_dp, &
                  106.0_dp, 105.0_dp, 108.0_dp, 110.0_dp]
  prices(:, 2) = [100.0_dp, 100.5_dp, 101.0_dp, 101.5_dp, &
                  102.0_dp, 102.5_dp, 103.0_dp, 103.5_dp]
  prices(:, 3) = [100.0_dp, 99.0_dp, 98.0_dp, 100.0_dp, &
                  101.0_dp, 103.0_dp, 102.0_dp, 104.0_dp]

  maxdd = pmaxdd(prices, max_drawdown=0.08_dp)
  call validate_solution(maxdd, 'PMaxDD')
  call assert_close(maxdd%weights(1), 1.0_dp, 2.0e-6_dp, 'PMaxDD reference weight')
  call assert_true(maxdd%risk_value <= 0.08_dp + 1.0e-8_dp, 'PMaxDD constraint')
  call assert_close(maxdd%terminal_return, 0.1_dp, 2.0e-7_dp, 'PMaxDD terminal return')

  avedd = pavedd(prices, average_drawdown=0.05_dp)
  call validate_solution(avedd, 'PAveDD')
  call assert_true(sum(avedd%drawdowns) / real(size(avedd%drawdowns), dp) <= &
                   0.05_dp + 1.0e-8_dp, 'PAveDD constraint')
  call assert_close(avedd%terminal_return, 0.1_dp, 5.0e-7_dp, 'PAveDD terminal return')

  cdar = pcdar(prices, alpha=0.9_dp, bound=0.08_dp)
  call validate_solution(cdar, 'PCDaR')
  call assert_true(cdar%risk_value <= 0.08_dp + 1.0e-8_dp, 'PCDaR constraint')
  call assert_close(cdar%terminal_return, 0.1_dp, 1.0e-6_dp, 'PCDaR terminal return')

  mincdar = pmincdar(prices, alpha=0.9_dp)
  call validate_solution(mincdar, 'PMinCDaR')
  call assert_true(mincdar%risk_value <= 1.0e-7_dp, 'minimum CDaR reference')
  call assert_close(mincdar%objective, mincdar%risk_value, 1.0e-14_dp, &
                    'minimum CDaR objective')

  soft = pmaxdd(prices, max_drawdown=0.08_dp, soft_budget=.true.)
  call assert_equal_int(soft%status, frapo_ok, 'soft-budget status')
  call assert_true(sum(soft%weights) <= 1.0_dp + 2.0e-7_dp, 'soft-budget inequality')

  print '(a)', 'test_drawdown: PASS'

contains
  subroutine validate_solution(result, label)
    type(portfolio_result), intent(in) :: result
    character(len=*), intent(in) :: label
    call assert_equal_int(result%status, frapo_ok, trim(label)//' status')
    call assert_close(sum(result%weights), 1.0_dp, 2.0e-7_dp, trim(label)//' budget')
    call assert_true(minval(result%weights) >= -2.0e-7_dp, trim(label)//' nonnegative')
    call assert_true(minval(result%drawdowns) >= -2.0e-8_dp, trim(label)//' drawdowns nonnegative')
  end subroutine validate_solution

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2es24.15)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a)') trim(label)//': failed'
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_equal_int(actual, expected, label)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,2i8)') trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_equal_int
end program test_drawdown
