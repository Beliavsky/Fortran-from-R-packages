! SPDX-License-Identifier: MIT
program test_models
  use yieldcurves
  implicit none
  real(dp), parameter :: m(10) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
    5.0_dp, 7.0_dp, 10.0_dp, 20.0_dp, 30.0_dp]
  real(dp) :: rates(10)
  real(dp), parameter :: knots(5) = [1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 30.0_dp]
  real(dp), parameter :: krates(5) = [0.045_dp, 0.043_dp, 0.042_dp, 0.040_dp, 0.043_dp]
  real(dp), parameter :: query(4) = [3.0_dp, 7.0_dp, 15.0_dp, 35.0_dp]
  real(dp), parameter :: natural_ref(4) = [0.0420458067632850_dp, 0.0415626121739130_dp, &
    0.0386245108695652_dp, 0.0452682065217391_dp]
  real(dp), parameter :: fmm_ref(4) = [0.0421317408906883_dp, 0.0416547141700405_dp, &
    0.0359113360323886_dp, 0.0624433198380567_dp]
  type(curve_t) :: fit
  type(series_t) :: pred
  integer :: i

  do i = 1, size(m)
    rates(i) = ns_rate_scalar(m(i), 0.04_dp, -0.02_dp, 0.03_dp, 2.0_dp)
  end do
  fit = yc_nelson_siegel(m, rates)
  call require(fit%ok, 'Nelson-Siegel fit status')
  call close_scalar(fit%tau, 2.0_dp, 1.0e-6_dp, 'Nelson-Siegel tau')
  call close_scalar(maxval(abs(fit%residuals)), 0.0_dp, 1.0e-10_dp, 'Nelson-Siegel residuals')

  do i = 1, size(m)
    rates(i) = sv_rate_scalar(m(i), 0.04_dp, -0.02_dp, 0.03_dp, -0.01_dp, 2.0_dp, 7.0_dp)
  end do
  fit = yc_svensson(m, rates)
  call require(fit%ok, 'Svensson fit status')
  call close_scalar(fit%tau1, 2.0_dp, 1.0e-5_dp, 'Svensson tau1')
  call close_scalar(fit%tau2, 7.0_dp, 1.0e-5_dp, 'Svensson tau2')
  call close_scalar(maxval(abs(fit%residuals)), 0.0_dp, 1.0e-10_dp, 'Svensson residuals')

  fit = yc_cubic_spline(knots, krates, 'natural')
  pred = yc_predict(fit, query)
  call close_vector(pred%y, natural_ref, 5.0e-14_dp, 'natural spline')

  fit = yc_cubic_spline(knots, krates, 'fmm')
  pred = yc_predict(fit, query)
  call close_vector(pred%y, fmm_ref, 5.0e-14_dp, 'FMM spline')

  fit = yc_curve(knots, krates)
  pred = yc_interpolate(fit, [3.0_dp, 7.0_dp], 'linear')
  call close_vector(pred%y, [0.0426666666666667_dp, 0.0412_dp], 1.0e-14_dp, 'linear interpolation')
  pred = yc_interpolate(fit, [3.0_dp], 'log_linear')
  call close_scalar(pred%y(1), 0.0424444444444444_dp, 1.0e-14_dp, 'log-linear interpolation')
  pred = yc_interpolate(fit, query, 'cubic')
  call close_vector(pred%y, natural_ref, 5.0e-14_dp, 'observed cubic interpolation')

  fit = yc_fit(knots, krates, 'cubic_spline')
  call require(fit%ok .and. trim(fit%method) == 'cubic_spline', 'fit dispatch')
  print '(a)', 'test_models: PASS'
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
end program test_models
