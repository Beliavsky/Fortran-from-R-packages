program test_bootstrap_custom
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_positive_inf
  use fitdistrplus
  use test_support
  implicit none
  type(distribution_model)::dist,custom
  type(fit_result)::fit
  type(bootstrap_result)::boot
  type(censored_sample)::sample
  type(npmle_result)::npmle
  type(fit_control)::ctl
  real(dp)::x(40),neg_inf,pos_inf
  real(dp),allocatable::band(:,:)
  integer::i,status

  call make_exponential(dist)
  do i=1,size(x);x(i)=dist%quantile((real(i,dp)-0.5_dp)/real(size(x),dp),[2.0_dp]);end do
  ctl%max_iterations=2500;ctl%calculate_vcov=.false.
  call mledist(x,dist,[1.0_dp],fit,ctl)
  call bootdist(x,dist,fit,12,boot,.true.,12345,ctl)
  call assert_true(boot%status==fit_success .and. boot%successful>=10,"parametric bootstrap")
  call assert_close(boot%mean(1),2.0_dp,0.25_dp,"bootstrap mean")
  call cdf_bootstrap_band([0.25_dp,0.50_dp],dist,boot,[0.025_dp,0.975_dp],band,status)
  call assert_true(status==fit_success .and. all(band>=0.0_dp) .and. all(band<=1.0_dp),"CDF bootstrap band")

  custom%name="laplace";custom%npar=2;custom%discrete=.false.
  custom%logpdf=>laplace_logpdf_test;custom%cdf=>laplace_cdf_test
  allocate(custom%default_lower(2),custom%default_upper(2),custom%parameter_names(2))
  custom%default_lower=[-huge(1.0_dp),tiny(1.0_dp)]
  custom%default_upper=[huge(1.0_dp),huge(1.0_dp)]
  custom%parameter_names=[character(len=24)::"location","scale"]
  do i=1,size(x);x(i)=0.3_dp+0.8_dp*sign(1.0_dp,real(2*mod(i,2)-1,dp))*log(1.0_dp+real(i,dp)/40.0_dp);end do
  call mledist(x,custom,[0.0_dp,1.0_dp],fit,ctl)
  call assert_true(fit%convergence==fit_success,"custom callback MLE")
  call assert_close(fit%estimate(1),0.3_dp,0.15_dp,"custom Laplace location")

  neg_inf=ieee_value(0.0_dp,ieee_negative_inf);pos_inf=ieee_value(0.0_dp,ieee_positive_inf)
  sample%left=[neg_inf,1.0_dp,2.0_dp,2.0_dp,3.0_dp]
  sample%right=[1.0_dp,1.0_dp,2.0_dp,3.0_dp,pos_inf]
  call turnbull_npmle(sample,npmle)
  call assert_true(npmle%status==fit_success,"Turnbull NPMLE")
  call assert_close(sum(npmle%probability),1.0_dp,1.0e-12_dp,"NPMLE probability sum")
  write(*,'(a)')"test_bootstrap_custom: PASS"
end program test_bootstrap_custom
