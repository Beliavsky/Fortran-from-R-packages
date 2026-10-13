! SPDX-License-Identifier: GPL-2.0-only
! Generated from independently inverted R beta-product distributions.
program test_bayesian_reference
  use propcis
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(ci_result) :: ci
  real(dp) :: shapes(4), lower, upper, cdf1, cdf2
  integer :: status, count
  count = 0
  shapes = [ &
    1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp ]
  ci = diffci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = -7.76393202248979009e-01_dp
  upper = 7.76393202252456005e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 1')
  cdf1 = bayes_cdf(1, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(1, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp ]
  ci = rrci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 5.00000000000000028e-02_dp
  upper = 2.00000000004095000e+01_dp
  call check(ci, lower, upper, 'Independent posterior case 2')
  cdf1 = bayes_cdf(2, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(2, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 1.00000000000000000e+00_dp ]
  ci = orci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 5.99022342473039000e-03_dp
  upper = 1.66938681430736011e+02_dp
  call check(ci, lower, upper, 'Independent posterior case 3')
  cdf1 = bayes_cdf(3, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(3, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp ]
  ci = diffci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = -9.22992316542723956e-01_dp
  upper = 9.22992316542722957e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 4')
  cdf1 = bayes_cdf(1, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(1, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp ]
  ci = rrci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 3.80182556080956979e-03_dp
  upper = 2.63031531543968981e+02_dp
  call check(ci, lower, upper, 'Independent posterior case 5')
  cdf1 = bayes_cdf(2, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(2, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp, 5.00000000000000000e-01_dp ]
  ci = orci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 1.26310147775409993e-04_dp
  upper = 7.91702026806047979e+03_dp
  call check(ci, lower, upper, 'Independent posterior case 6')
  cdf1 = bayes_cdf(3, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(3, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.50000000000000000e+00_dp, 1.15000000000000000e+01_dp, 7.20000000000000018e+00_dp, 4.09999999999999964e+00_dp ]
  ci = diffci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = -7.58171252870026957e-01_dp
  upper = -1.04094985690640005e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 7')
  cdf1 = bayes_cdf(1, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(1, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.50000000000000000e+00_dp, 1.15000000000000000e+01_dp, 7.20000000000000018e+00_dp, 4.09999999999999964e+00_dp ]
  ci = rrci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 5.11055816093518023e-02_dp
  upper = 7.60382535720829034e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 8')
  cdf1 = bayes_cdf(2, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(2, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.50000000000000000e+00_dp, 1.15000000000000000e+01_dp, 7.20000000000000018e+00_dp, 4.09999999999999964e+00_dp ]
  ci = orci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 1.21367699360339001e-02_dp
  upper = 6.26444147930944051e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 9')
  cdf1 = bayes_cdf(3, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(3, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    6.00000000000000000e+00_dp, 1.60000000000000000e+01_dp, 3.10000000000000000e+01_dp, 1.00000000000000000e+00_dp ]
  ci = diffci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = -8.65125374547799031e-01_dp
  upper = -4.87020612063609026e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 10')
  cdf1 = bayes_cdf(1, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(1, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    6.00000000000000000e+00_dp, 1.60000000000000000e+01_dp, 3.10000000000000000e+01_dp, 1.00000000000000000e+00_dp ]
  ci = rrci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 1.16277832059997005e-01_dp
  upper = 4.89136277440124001e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 11')
  cdf1 = bayes_cdf(2, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(2, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    6.00000000000000000e+00_dp, 1.60000000000000000e+01_dp, 3.10000000000000000e+01_dp, 1.00000000000000000e+00_dp ]
  ci = orci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 2.56412816066956006e-04_dp
  upper = 5.88509662053620999e-02_dp
  call check(ci, lower, upper, 'Independent posterior case 12')
  cdf1 = bayes_cdf(3, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(3, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.00000000000000000e+01_dp, 2.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 9.00000000000000000e+00_dp ]
  ci = diffci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 5.45455036897044998e-01_dp
  upper = 9.60366678925729023e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 13')
  cdf1 = bayes_cdf(1, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(1, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.00000000000000000e+01_dp, 2.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 9.00000000000000000e+00_dp ]
  ci = rrci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 2.67716228418291990e+00_dp
  upper = 3.23603451375579994e+02_dp
  call check(ci, lower, upper, 'Independent posterior case 14')
  cdf1 = bayes_cdf(2, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(2, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.00000000000000000e+01_dp, 2.00000000000000000e+00_dp, 1.00000000000000000e+00_dp, 9.00000000000000000e+00_dp ]
  ci = orci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 1.45105308873310008e+01_dp
  upper = 6.77533127532362960e+03_dp
  call check(ci, lower, upper, 'Independent posterior case 15')
  cdf1 = bayes_cdf(3, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(3, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.00000000000000011e-01_dp, 4.00000000000000000e+00_dp, 3.00000000000000000e+00_dp, 2.99999999999999989e-01_dp ]
  ci = diffci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = -9.99836803801021023e-01_dp
  upper = -3.96689082399425019e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 16')
  cdf1 = bayes_cdf(1, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(1, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.00000000000000011e-01_dp, 4.00000000000000000e+00_dp, 3.00000000000000000e+00_dp, 2.99999999999999989e-01_dp ]
  ci = rrci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 1.96620280287568984e-09_dp
  upper = 3.90819817685369997e-01_dp
  call check(ci, lower, upper, 'Independent posterior case 17')
  cdf1 = bayes_cdf(2, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(2, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  shapes = [ &
    2.00000000000000011e-01_dp, 4.00000000000000000e+00_dp, 3.00000000000000000e+00_dp, 2.99999999999999989e-01_dp ]
  ci = orci_bayes(0, 0, 0, 0, shapes(1), shapes(2), shapes(3), shapes(4))
  lower = 2.92678890277155999e-12_dp
  upper = 8.00464192468516994e-02_dp
  call check(ci, lower, upper, 'Independent posterior case 18')
  cdf1 = bayes_cdf(3, ci%lower, shapes, status)
  if (status /= ci_ok) error stop 'lower CDF status'
  cdf2 = bayes_cdf(3, ci%upper, shapes, status, upper_tail=.true.)
  if (status /= ci_ok) error stop 'upper CDF status'
  if (abs(cdf1 - 0.025_dp) > 2.0e-7_dp .or. abs(cdf2 - 0.025_dp) > 2.0e-7_dp) error stop 'posterior tail mass'
  print *, "Independent posterior assertions passed:", count
contains
  subroutine check(actual, lo, hi, label)
    type(ci_result), intent(in) :: actual !! Computed equal-tail posterior interval and numerical status.
    real(dp), intent(in) :: lo !! Independent lower quantile from R integrate plus uniroot.
    real(dp), intent(in) :: hi !! Independent upper quantile from R integrate plus uniroot.
    character(len=*), intent(in) :: label !! Posterior reference case identifier.
    count = count + 1
    if (actual%status /= ci_ok .or. .not. ieee_is_finite(actual%lower) .or. .not. ieee_is_finite(actual%upper)) then
      print *, label, actual
      error stop 'Independent posterior status'
    end if
    if (abs(actual%lower - lo) > 2.0e-6_dp * max(abs(lo), 1.0e-14_dp) .or. &
        abs(actual%upper - hi) > 2.0e-6_dp * max(abs(hi), 1.0e-14_dp)) then
      print *, label, actual%lower, lo, actual%upper, hi
      error stop 'Independent posterior quantiles'
    end if
  end subroutine check
end program test_bayesian_reference
