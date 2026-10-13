! SPDX-License-Identifier: GPL-2.0-only
program test_properties
  use propcis
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite, ieee_is_nan
  implicit none
  type(ci_result) :: ci, swapped, tail, shortest
  real(dp) :: quiet_nan, positive_inf, v, w, near_one
  integer :: status, k, x
  quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)
  positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
  near_one = nearest(1.0_dp, -1.0_dp)

  ! Invalid domains are explicit status values rather than process termination.
  ci = exactci(-1, 10, 0.95_dp)
  call invalid(ci)
  ci = add4ci(2, 0, 0.95_dp)
  call invalid(ci)
  ci = addz2ci(1, 10, quiet_nan)
  call invalid(ci)
  ci = scoreci(1, 10, positive_inf)
  call invalid(ci)
  ci = scoreci_mp(6, 5, 10, 0.95_dp)
  call invalid(ci)
  ci = oddsratioci_mp(0, 0, 0.95_dp)
  call invalid(ci)
  ci = wald2ci(1, 10, 1, 10, 0.95_dp, 'bad')
  call invalid(ci)
  ci = blakerci(1, 10, 0.95_dp, 0.0_dp)
  call invalid(ci)
  ci = diffci_bayes(1, 10, 1, 10, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call invalid(ci)
  ci = rrci_bayes(1, 10, 1, 10, 1.0_dp, quiet_nan, 1.0_dp, 1.0_dp)
  call invalid(ci)
  ci = orci_bayes(11, 10, 1, 10, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call invalid(ci)
  v = bayes_cdf(9, 0.0_dp, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], status)
  if (status /= ci_invalid .or. .not. ieee_is_nan(v)) error stop 'invalid posterior kind'
  v = acceptbin(1, 10, quiet_nan)
  if (.not. ieee_is_nan(v)) error stop 'invalid binomial probability'
  v = z2stat(quiet_nan, 10, 0.5_dp, 10, 0.0_dp)
  if (.not. ieee_is_nan(v)) error stop 'invalid score probability'
  v = limit(1, 10, 2, 10, 0.95_dp, 0.5_dp, 2)
  if (.not. ieee_is_nan(v)) error stop 'invalid grid direction'

  ! Regression tests for verified upstream defects and boundary roundoff.
  do k = 1, 20
    ci = oddsratioci_mp(0, k, 0.95_dp)
    if (ci%status /= ci_ok .or. ci%upper /= positive_inf) error stop 'paired OR infinite upper bound'
    ci = oddsratioci_mp(k, 0, 0.95_dp)
    if (ci%status /= ci_ok .or. ci%lower /= 0.0_dp) error stop 'paired OR zero lower bound'
  end do
  ci = wald2ci(5, 20, 10, 30, 0.95_dp)
  swapped = wald2ci(5, 20, 10, 30, 0.95_dp, 'AC')
  call close(ci%lower, swapped%lower, 1.0e-14_dp)
  ci = riskscoreci(10, 100, 90, 100, 0.95_dp)
  if (ci%lower < 0.05_dp .or. ci%upper > 0.2_dp) error stop 'RR lower-root bracket regression'
  ci = orci_bayes(5, 20, 30, 30, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  if (ci%status /= ci_ok) error stop 'all-success Bayesian OR'
  call close(ci%lower, 0.000256412816066956_dp, 2.0e-9_dp)
  call close(ci%upper, 0.0588509662053621_dp, 2.0e-9_dp)
  v = limit(5, 10, 3, 10, 0.95_dp, 1.0_dp, 0)
  if (.not. ieee_is_finite(v)) error stop 'OR grid starting at one'

  ! Confidence tails below supported numerical resolution must not report success.
  ci = exactci(1, 100, near_one)
  if (ci%status /= ci_numerical) error stop 'extreme beta quantile status'
  ci = blakerci(1, 100, near_one)
  if (ci%status /= ci_numerical) error stop 'extreme Blaker status'
  ci = midpci(1, 100, near_one)
  if (ci%status /= ci_numerical) error stop 'extreme mid-P status'
  ci = diffci_bayes(0, 0, 0, 0, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, near_one)
  if (ci%status /= ci_numerical) error stop 'extreme Bayesian tail status'
  ci = diffci_bayes(0, 0, 0, 0, 1.0e-300_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  if (ci%status /= ci_numerical) error stop 'extreme shape status'
  ci = addz2ci(1, 100, near_one)
  if (ci%status /= ci_ok .or. .not. ieee_is_finite(ci%upper)) error stop 'normal critical tail stability'

  ! Complement and sample-exchange identities, including zero and all successes.
  do x = 0, 20
    ci = exactci(x, 20, 0.95_dp)
    swapped = exactci(20 - x, 20, 0.95_dp)
    call close(ci%lower, 1.0_dp - swapped%upper, 2.0e-12_dp)
    ci = diffscoreci(x, 20, 7, 20, 0.95_dp)
    swapped = diffscoreci(7, 20, x, 20, 0.95_dp)
    call close(ci%lower, -swapped%upper, 3.0e-7_dp)
    ci = riskscoreci(x, 20, 7, 20, 0.95_dp)
    swapped = riskscoreci(7, 20, x, 20, 0.95_dp)
    if (x == 0) then
      call close(ci%lower, 0.0_dp, 1.0e-14_dp)
    else
      call close(ci%lower, 1.0_dp / swapped%upper, 1.0e-8_dp)
    end if
    call close(ci%upper, 1.0_dp / swapped%lower, 1.0e-8_dp)
  end do

  ! Uniform independent priors give a triangular difference and elementary RR law.
  ci = diffci_bayes(0, 0, 0, 0, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call close(ci%lower, -1.0_dp + sqrt(0.05_dp), 2.0e-9_dp)
  call close(ci%upper, 1.0_dp - sqrt(0.05_dp), 2.0e-9_dp)
  ci = rrci_bayes(0, 0, 0, 0, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call close(ci%lower, 0.05_dp, 2.0e-9_dp)
  call close(ci%upper, 20.0_dp, 2.0e-7_dp)
  shortest = diffci_bayes_hpd(0, 0, 0, 0, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call close(shortest%lower, -shortest%upper, 2.0e-9_dp)
  tail = diffci_bayes(0, 0, 0, 0, 2.5_dp, 11.5_dp, 7.2_dp, 4.1_dp)
  shortest = diffci_bayes_hpd(0, 0, 0, 0, 2.5_dp, 11.5_dp, 7.2_dp, 4.1_dp)
  if (shortest%status /= ci_ok) error stop 'asymmetric connected HPD status'
  if (shortest%upper - shortest%lower > tail%upper - tail%lower + 2.0e-8_dp) error stop 'shortest interval width'
  v = bayes_cdf(bayes_difference, shortest%lower, [2.5_dp, 11.5_dp, 7.2_dp, 4.1_dp], status)
  w = bayes_cdf(bayes_difference, shortest%upper, [2.5_dp, 11.5_dp, 7.2_dp, 4.1_dp], status)
  call close(w - v, 0.95_dp, 2.0e-7_dp)
  swapped = diffci_bayes_hpd(0, 0, 0, 0, 7.2_dp, 4.1_dp, 2.5_dp, 11.5_dp)
  call close(shortest%lower, -swapped%upper, 5.0e-5_dp)
  print *, 'Boundary, invalid-input, known-defect, and distribution-identity checks passed.'
contains
  subroutine invalid(value)
    type(ci_result), intent(in) :: value !! Interval from deliberately invalid input; requires status and NaN bounds.
    if (value%status /= ci_invalid) error stop 'missing invalid-input status'
    if (.not. ieee_is_nan(value%lower) .or. .not. ieee_is_nan(value%upper)) error stop 'invalid bounds must be NaN'
  end subroutine invalid
  subroutine close(actual, expected, tolerance)
    real(dp), intent(in) :: actual !! Computed finite statistic or bound.
    real(dp), intent(in) :: expected !! Analytic or independently computed reference value.
    real(dp), intent(in) :: tolerance !! Allowed absolute error.
    if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
      print *, actual, expected, tolerance
      error stop 'property mismatch'
    end if
  end subroutine close
end program test_properties
