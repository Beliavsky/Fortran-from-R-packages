! SPDX-License-Identifier: GPL-2.0-or-later
module test_grouped_risk_support
  use actuar, only : dp
  implicit none
contains
  pure function exponential_mgf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y = 1.0_dp / (1.0_dp - x)
  end function exponential_mgf
end module test_grouped_risk_support

program test_grouped_risk
  use actuar, only : dp, grouped_mean, grouped_variance, grouped_quantile, &
    ogive_cdf, apply_coverage, ruin_exponential, adjustment_coefficient_poisson
  use test_grouped_risk_support, only : exponential_mgf
  implicit none
  real(dp) :: breaks(4), freq(3), r
  real(dp), allocatable :: paid(:)

  breaks = [0.0_dp, 10.0_dp, 20.0_dp, 30.0_dp]
  freq = [2.0_dp, 3.0_dp, 5.0_dp]
  call assert_close(grouped_mean(breaks, freq), 18.0_dp, 1.0e-14_dp)
  call assert_close(grouped_variance(breaks, freq), 61.0_dp, 1.0e-13_dp)
  call assert_close(ogive_cdf(15.0_dp, breaks, freq), 0.35_dp, 1.0e-14_dp)
  call assert_close(grouped_quantile(0.35_dp, breaks, freq), 15.0_dp, 1.0e-13_dp)
  paid = apply_coverage([5.0_dp, 15.0_dp, 40.0_dp], deductible=10.0_dp, &
    limit=30.0_dp, coinsurance=0.8_dp)
  call assert_close(sum(paid), 20.0_dp, 1.0e-14_dp)
  call assert_close(ruin_exponential(10.0_dp, 1.0_dp, 2.0_dp, 1.0_dp), &
    0.5_dp * exp(-10.0_dp), 1.0e-13_dp)
  r = adjustment_coefficient_poisson(exponential_mgf, 1.0_dp, 1.2_dp, 1.9_dp)
  call assert_close(r, 1.0_dp / 6.0_dp, 1.0e-10_dp)

  print '(a)', 'test_grouped_risk: PASS'
contains
  subroutine assert_close(actual, expected, tol)
    real(dp), intent(in) :: actual, expected, tol
    if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
      print '(a,3es24.15)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_grouped_risk
