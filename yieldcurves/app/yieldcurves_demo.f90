! SPDX-License-Identifier: MIT
program yieldcurves_demo
  use yieldcurves
  implicit none
  real(dp), parameter :: maturities(10) = [0.25_dp,0.5_dp,1.0_dp,2.0_dp,3.0_dp, &
    5.0_dp,7.0_dp,10.0_dp,20.0_dp,30.0_dp]
  real(dp), parameter :: rates(10) = [0.052_dp,0.050_dp,0.048_dp,0.045_dp,0.043_dp, &
    0.042_dp,0.041_dp,0.040_dp,0.042_dp,0.043_dp]
  type(curve_t) :: fit
  type(series_t) :: forecast
  type(carry_result_t) :: carry
  integer :: i

  fit = yc_nelson_siegel(maturities, rates)
  if (.not. fit%ok) error stop trim(fit%message)
  forecast = yc_predict(fit, [1.0_dp,2.0_dp,5.0_dp,10.0_dp,30.0_dp])
  carry = yc_carry(fit, [2.0_dp,5.0_dp,10.0_dp])

  print '(a)', 'Nelson-Siegel parameters'
  print '(a,f10.6)', 'beta0 = ', fit%beta0
  print '(a,f10.6)', 'beta1 = ', fit%beta1
  print '(a,f10.6)', 'beta2 = ', fit%beta2
  print '(a,f10.6)', 'tau   = ', fit%tau
  print '(a)', 'Predicted zero rates'
  do i=1,size(forecast%x)
    print '(f7.2,2x,f10.6)', forecast%x(i), forecast%y(i)
  end do
  print '(a)', 'Carry and roll-down'
  do i=1,size(carry%maturity)
    print '(f7.2,3(2x,f10.6))', carry%maturity(i), carry%carry(i), carry%rolldown(i), carry%total(i)
  end do
end program yieldcurves_demo
