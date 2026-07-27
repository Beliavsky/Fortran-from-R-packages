! SPDX-License-Identifier: MIT
program curve_fitting
  use yieldcurves
  implicit none
  real(dp), parameter :: m(10) = [0.25_dp,0.5_dp,1.0_dp,2.0_dp,3.0_dp,5.0_dp,7.0_dp,10.0_dp,20.0_dp,30.0_dp]
  real(dp), parameter :: r(10) = [0.052_dp,0.050_dp,0.048_dp,0.045_dp,0.043_dp,0.042_dp,0.041_dp,0.040_dp,0.042_dp,0.043_dp]
  type(curve_t) :: ns, sv, spline
  type(series_t) :: p

  ns = yc_nelson_siegel(m,r)
  sv = yc_svensson(m,r)
  spline = yc_cubic_spline(m,r,'natural')
  if (.not. ns%ok) error stop trim(ns%message)
  if (.not. sv%ok) error stop trim(sv%message)
  if (.not. spline%ok) error stop trim(spline%message)

  p = yc_predict(ns,[3.0_dp,8.0_dp,15.0_dp])
  print '(a,3f12.7)', 'NS rates:       ',p%y
  p = yc_predict(sv,[3.0_dp,8.0_dp,15.0_dp])
  print '(a,3f12.7)', 'Svensson rates: ',p%y
  p = yc_predict(spline,[3.0_dp,8.0_dp,15.0_dp])
  print '(a,3f12.7)', 'Spline rates:   ',p%y
end program curve_fitting
