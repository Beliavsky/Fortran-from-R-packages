program test_robust_prewhiten
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpeif, only : dp, robust_clean, robust_location_scale, ar_prewhiten
  implicit none
  real(dp) :: x(12), location, scale, intercept
  real(dp), allocatable :: cleaned(:), residuals(:), coefficients(:)
  integer :: status, i

  x = [(-0.02_dp + 0.004_dp * real(i - 1, dp), i = 1, 11), 1.5_dp]
  call robust_location_scale(x, location, scale, 'mopt', 0.95_dp)
  call assert_true(scale > 0.0_dp, 'robust scale')
  call robust_clean(x, cleaned, 0.95_dp, 'mopt', status)
  call assert_true(all(ieee_is_finite(cleaned)), 'cleaned finite')
  call assert_true(cleaned(12) < 0.5_dp, 'outlier clipped')

  do i = 2, size(x)
    x(i) = 0.7_dp * x(i - 1) + 0.01_dp * sin(real(i, dp))
  end do
  call ar_prewhiten(x, 1, residuals, coefficients, intercept, status)
  call assert_true(status == 0, 'AR prewhiten status')
  call assert_true(size(coefficients) == 1, 'AR coefficient size')
  call assert_true(all(ieee_is_finite(residuals)), 'AR residual finite')
  call assert_true(abs(sum(residuals(2:))) < 1.0e-10_dp, 'AR residual mean')

  print '(a)', 'test_robust_prewhiten: PASS'
contains
  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label)//' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_robust_prewhiten
