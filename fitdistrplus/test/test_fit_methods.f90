program test_fit_methods
  use fitdistrplus
  use test_support
  implicit none
  type(distribution_model)::dist
  type(fit_result)::fit
  type(fit_control)::ctl
  real(dp),allocatable::x(:),start(:)
  real(dp)::probs(2)
  integer::i,status,orders(2)

  ctl%max_iterations=4000;ctl%tolerance=1.0e-9_dp
  allocate(x(80))
  call make_normal(dist)
  do i=1,size(x)
    x(i)=dist%quantile((real(i,dp)-0.5_dp)/real(size(x),dp),[0.5_dp,1.7_dp])
  end do
  start=[0.0_dp,1.0_dp]
  call mledist(x,dist,start,fit,ctl)
  call assert_true(fit%convergence==fit_success,"normal MLE convergence")
  call assert_close(fit%estimate(1),0.5_dp,2.0e-3_dp,"normal MLE mean")
  call assert_close(fit%estimate(2),1.67_dp,2.0e-2_dp,"normal MLE sd")
  call assert_true(size(fit%covariance,1)==2,"normal MLE covariance")

  probs=[0.25_dp,0.75_dp]
  call qmedist(x,dist,start,probs,fit,ctl)
  call assert_true(fit%convergence==fit_success,"normal QME convergence")
  call assert_close(fit%estimate(1),0.5_dp,3.0e-3_dp,"normal QME mean")
  call assert_close(fit%estimate(2),1.7_dp,2.0e-2_dp,"normal QME sd")

  orders=[1,2]
  call mmedist(x,dist,start,orders,fit,ctl)
  call assert_true(fit%convergence==fit_success,"normal MME convergence")
  call assert_close(fit%estimate(1),sum(x)/real(size(x),dp),2.0e-3_dp,"normal MME mean")

  call mgedist(x,dist,start,fit,gof_cvm,ctl)
  call assert_true(fit%convergence==fit_success,"normal MGE convergence")
  call assert_close(fit%estimate(1),0.5_dp,1.0e-2_dp,"normal MGE mean")

  call msedist(x,dist,start,fit,phi_kl,control=ctl)
  call assert_true(fit%convergence==fit_success,"normal MSE convergence")
  call assert_close(fit%estimate(1),0.5_dp,3.0e-2_dp,"normal MSE mean")

  call default_start(x,dist,start,status)
  call assert_true(status==fit_success .and. size(start)==2,"default start")
  write(*,'(a)')"test_fit_methods: PASS"
end program test_fit_methods
