program test_bootstrap
  use mbbefd, only : dp, dr_fit_result, dr_boot_result, fit_dr, boot_dr, roiunif
  implicit none
  type(dr_fit_result)::fit
  type(dr_boot_result)::boot
  real(dp)::x(100)
  call roiunif(x,0.2_dp);call fit_dr(x,'oiunif',fit)
  call boot_dr(fit,x,20,boot,parametric=.true.)
  if(boot%status/=0.or.boot%successful<15) error stop 'bootstrap'
  if(any(boot%confidence_interval(:,2)>boot%confidence_interval(:,3))) error stop 'CI order'
  print '(a)', 'test_bootstrap: PASS'
end program test_bootstrap
