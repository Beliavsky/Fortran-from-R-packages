program test_censored
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_positive_inf
  use fitdistrplus
  use test_support
  implicit none
  type(distribution_model)::dist
  type(censored_sample)::sample
  type(fit_result)::fit
  type(fit_control)::ctl
  real(dp)::neg_inf,pos_inf,value
  integer::i,n

  neg_inf=ieee_value(0.0_dp,ieee_negative_inf)
  pos_inf=ieee_value(0.0_dp,ieee_positive_inf)
  n=120;allocate(sample%left(n),sample%right(n));call make_exponential(dist)
  do i=1,n
    value=dist%quantile((real(i,dp)-0.5_dp)/real(n,dp),[1.5_dp])
    if(value<0.20_dp)then
      sample%left(i)=neg_inf;sample%right(i)=0.20_dp
    else if(value>1.20_dp)then
      sample%left(i)=1.20_dp;sample%right(i)=pos_inf
    else if(mod(i,7)==0)then
      sample%left(i)=max(0.0_dp,value-0.025_dp);sample%right(i)=value+0.025_dp
    else
      sample%left(i)=value;sample%right(i)=value
    end if
  end do
  ctl%max_iterations=3000;ctl%tolerance=1.0e-9_dp
  call mledist_censored(sample,dist,[1.0_dp],fit,ctl)
  call assert_true(fit%convergence==fit_success,"censored MLE convergence")
  call assert_close(fit%estimate(1),1.5_dp,6.0e-2_dp,"censored exponential rate")
  write(*,'(a)')"test_censored: PASS"
end program test_censored
