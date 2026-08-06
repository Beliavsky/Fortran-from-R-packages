program test_stats
  use fitdistrplus
  use test_support
  implicit none
  type(distribution_model)::dist
  type(fit_result)::fit
  type(descriptive_result)::description
  type(gof_result)::gof
  real(dp)::x(10),breaks(3)
  integer::i,status
  real(dp),allocatable::lower(:),upper(:)

  x=[(real(i,dp),i=1,10)]
  call descdist(x,description)
  call assert_true(description%status==fit_success,"descdist status")
  call assert_close(description%mean,5.5_dp,1.0e-14_dp,"descdist mean")
  call assert_close(description%median,5.5_dp,1.0e-14_dp,"descdist median")

  call make_uniform(dist)
  fit%estimate=[0.0_dp,11.0_dp];fit%aic=10.0_dp;fit%bic=11.0_dp
  breaks=[3.0_dp,6.0_dp,8.0_dp]
  call gofstat(x,dist,fit,gof,breaks)
  call assert_true(gof%status==fit_success,"gof status")
  call assert_true(gof%ks<0.12_dp .and. gof%cvm<0.1_dp,"uniform gof")
  call assert_true(gof%chi_square_df==1,"chi-square degrees of freedom")

  call make_gamma(dist);call detectbound(dist,lower,upper,status)
  call assert_true(status==fit_success .and. lower(1)>0.0_dp,"detectbound")
  write(*,'(a)')"test_stats: PASS"
end program test_stats
