! SPDX-License-Identifier: GPL-3.0-only
program gaussian_ar1_example
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use imputefin
  implicit none
  real(dp) :: y(12)
  type(ar1_fit_result) :: fit
  type(imputation_result) :: imp
  type(imputation_options) :: io
  y=[1.00_dp,1.12_dp,1.18_dp,ieee_value(0.0_dp,ieee_quiet_nan), &
     ieee_value(0.0_dp,ieee_quiet_nan),1.45_dp,1.50_dp,1.62_dp, &
     1.70_dp,1.76_dp,1.82_dp,1.91_dp]
  call fit_ar1_gaussian(y,fit,return_conditional=.true.)
  io%seed=1234_8
  call impute_ar1_gaussian(y,imp,impute_options_in=io)
  write(*,'(a,3f12.6)')'phi0, phi1, sigma2: ',fit%phi0,fit%phi1,fit%sigma2
  write(*,'(a,*(f9.4,1x))')'imputed: ',imp%values(:,1,1)
end program gaussian_ar1_example
