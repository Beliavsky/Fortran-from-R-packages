! SPDX-License-Identifier: MIT
program test_risk_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use yieldcurves
  implicit none
  real(dp), parameter :: m(5) = [1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 30.0_dp]
  real(dp), parameter :: r(5) = [0.040_dp, 0.042_dp, 0.044_dp, 0.045_dp, 0.046_dp]
  real(dp), parameter :: krd_ref(5) = [0.0439513644309006_dp, 0.216559482794440_dp, &
    0.797625910195198_dp, 6.62179002412338_dp, 0.0_dp]
  type(curve_t) :: curve, short_curve
  type(zspread_result_t) :: zs
  type(series_t) :: krd
  type(carry_result_t) :: carry
  type(slope_result_t) :: slopes
  type(factor_result_t) :: factors
  real(dp) :: price
  real(dp), allocatable :: times(:), cfs(:), z(:)
  integer :: i

  curve = yc_curve(m, r)
  allocate(times(20), cfs(20), z(20))
  do i = 1, 20
    times(i) = real(i,dp)/2.0_dp
    cfs(i) = 2.0_dp
    z(i) = interp(m,r,times(i))
  end do
  cfs(20) = cfs(20) + 100.0_dp
  price = sum(cfs*(1.0_dp+z)**(-times))
  zs = yc_zspread(price, 0.04_dp, 10.0_dp, curve)
  call require(zs%ok, 'z-spread status')
  call close_scalar(zs%zspread, 0.0_dp, 1.0e-10_dp, 'zero z-spread')
  call close_scalar(zs%model_price, price, 2.0e-10_dp, 'z-spread model price')

  krd = yc_key_rate_duration(0.05_dp, 10.0_dp, curve)
  call require(krd%ok, 'key-rate-duration status')
  call close_vector(krd%y, krd_ref, 2.0e-10_dp, 'key-rate durations')

  curve = yc_nelson_siegel(m, [0.04_dp,0.04_dp,0.04_dp,0.04_dp,0.04_dp])
  carry = yc_carry(curve, [2.0_dp,5.0_dp,10.0_dp])
  call close_vector(carry%rolldown, [0.0_dp,0.0_dp,0.0_dp], 1.0e-8_dp, 'flat-curve roll-down')
  call close_vector(carry%total, carry%carry+carry%rolldown, 1.0e-15_dp, 'carry identity')

  short_curve = yc_curve([1.0_dp,2.0_dp,5.0_dp,10.0_dp], [0.04_dp,0.042_dp,0.043_dp,0.045_dp])
  slopes = yc_slope(short_curve)
  call require(ieee_is_nan(slopes%spread_2s30s), 'out-of-range slope should be NaN')
  factors = yc_level_slope_curvature(short_curve)
  call close_scalar(factors%level, 0.0425_dp, 1.0e-15_dp, 'empirical level')
  call close_scalar(factors%slope, -0.005_dp, 1.0e-15_dp, 'empirical slope')
  print '(a)', 'test_risk_analysis: PASS'
contains
  pure real(dp) function interp(x,y,q) result(v)
    real(dp),intent(in)::x(:),y(:),q
    integer::j
    if(q<=x(1))then;v=y(1);return;end if
    if(q>=x(size(x)))then;v=y(size(y));return;end if
    do j=1,size(x)-1
      if(q<=x(j+1))then
        v=y(j)+(y(j+1)-y(j))*(q-x(j))/(x(j+1)-x(j));return
      end if
    end do
    v=y(size(y))
  end function interp
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
end program test_risk_analysis
