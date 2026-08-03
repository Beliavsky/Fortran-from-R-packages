! SPDX-License-Identifier: GPL-3.0-only
program student_t_ar1_example
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use imputefin
  implicit none
  real(dp) :: y(15)
  type(ar1_fit_result) :: fit
  type(imputation_result) :: imp
  type(ar1_options) :: opt
  type(imputation_options) :: io
  y=[0.0_dp,0.3_dp,0.1_dp,0.4_dp,8.0_dp,0.2_dp, &
     ieee_value(0.0_dp,ieee_quiet_nan),0.1_dp,-0.1_dp,0.0_dp, &
     0.2_dp,0.1_dp,-0.2_dp,0.0_dp,0.1_dp]
  opt%remove_outliers=.true.;opt%outlier_prob_th=1.0e-3_dp
  call fit_ar1_t(y,fit,opt)
  io%n_burn=50;io%seed=4321_8
  call impute_ar1_t(y,imp,opt,io)
  write(*,'(a,4f12.6)')'phi0, phi1, sigma2, nu: ',fit%phi0,fit%phi1,fit%sigma2,fit%nu
  write(*,'(a,*(f9.4,1x))')'cleaned/imputed: ',imp%values(:,1,1)
end program student_t_ar1_example
