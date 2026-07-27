! SPDX-License-Identifier: MIT
program test_conversions
  use yieldcurves
  implicit none
  real(dp), parameter :: m(5) = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 10.0_dp]
  real(dp), parameter :: par(5) = [0.040_dp, 0.042_dp, 0.043_dp, 0.044_dp, 0.045_dp]
  real(dp), parameter :: zero_ref(5) = [0.04_dp, 0.0420420833185258_dp, 0.0430724866011889_dp, &
    0.0441362559737479_dp, 0.0452818493700884_dp]
  type(series_t) :: zeros, back, discount, forward, semi_zero, semi_back
  type(duration_result_t) :: zero_duration
  type(curve_t) :: curve
  type(bond_duration_result_t) :: bond

  zeros = yc_par_to_zero(m, par)
  call require(zeros%ok, 'par-to-zero status')
  call close_vector(zeros%y, zero_ref, 2.0e-13_dp, 'par-to-zero')
  back = yc_zero_to_par(m, zeros%y)
  call close_vector(back%y, par, 1.0e-4_dp, 'zero-to-par round trip')

  semi_zero = yc_par_to_zero([0.5_dp,1.0_dp,1.5_dp,2.0_dp,3.0_dp,5.0_dp], &
    [0.040_dp,0.042_dp,0.043_dp,0.044_dp,0.046_dp,0.048_dp], 2)
  semi_back = yc_zero_to_par(semi_zero%x, semi_zero%y, 2)
  call close_vector(semi_back%y, [0.040_dp,0.042_dp,0.043_dp,0.044_dp,0.046_dp,0.048_dp], &
    1.0e-3_dp, 'semi-annual round trip')

  curve = yc_curve([1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp], &
    [0.045_dp, 0.043_dp, 0.042_dp, 0.040_dp])
  discount = yc_discount(curve, [1.0_dp, 5.0_dp], 'continuous')
  call close_vector(discount%y, [exp(-0.045_dp), exp(-0.042_dp*5.0_dp)], 1.0e-14_dp, 'discount factors')

  curve = yc_nelson_siegel([0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp], &
    [0.052_dp, 0.050_dp, 0.048_dp, 0.045_dp, 0.042_dp, 0.040_dp])
  zero_duration = yc_duration(curve, [1.0_dp,5.0_dp], 'annual')
  call close_vector(zero_duration%macaulay, [1.0_dp,5.0_dp], 1.0e-15_dp, 'zero Macaulay duration')

  forward = yc_forward(curve, [1.0_dp, 5.0_dp])
  call close_scalar(forward%y(1), ns_forward_scalar(1.0_dp, curve%beta0, curve%beta1, &
    curve%beta2, curve%tau), 1.0e-13_dp, 'instantaneous forward')

  bond = yc_bond_duration(100.0_dp, 0.05_dp, 2.0_dp, 0.04_dp, 2, 'semi_annual')
  call require(bond%ok, 'bond duration status')
  call close_scalar(bond%price, 101.90386434933714_dp, 1.0e-12_dp, 'bond price')
  call close_scalar(bond%macaulay_duration, 1.928782921964222_dp, 1.0e-12_dp, 'Macaulay duration')
  call close_scalar(bond%modified_duration, 1.890963648984531_dp, 1.0e-12_dp, 'modified duration')
  call close_scalar(bond%convexity, 4.578046580359792_dp, 1.0e-12_dp, 'bond convexity')
  print '(a)', 'test_conversions: PASS'
contains
  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'failed: '//trim(label)
      error stop 1
    end if
  end subroutine require
  subroutine close_scalar(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual-expected) > tolerance) then
      write(*,'(a,3es24.15)') trim(label)//': ', actual, expected, abs(actual-expected)
      error stop 1
    end if
  end subroutine close_scalar
  subroutine close_vector(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: label
    if (size(actual) /= size(expected) .or. maxval(abs(actual-expected)) > tolerance) then
      write(*,'(a,es24.15)') trim(label)//' maximum error: ', maxval(abs(actual-expected))
      error stop 1
    end if
  end subroutine close_vector
end program test_conversions
