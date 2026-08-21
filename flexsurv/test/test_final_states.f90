program test_final_states
  use flexsurv_kinds, only : dp
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_fit, only : initialize_spec, initial_theta
  use flexsurv_multistate, only : flexsurv_transition, msm_parametric
  use flexsurv_final_states, only : final_state_summary, simfinal_fmsm, &
    simfinal_fmsm_ci, pfinal_fmsm
  implicit none
  type(flexsurv_transition)::tr(2)
  type(final_state_summary)::fs
  real(dp)::pf(3)
  real(dp),allocatable::plo(:),phi(:),mlo(:),mhi(:),qlo(:,:),qhi(:,:)
  integer::i2,i3
  call make_transition(tr(1),1,2,0.4_dp)
  call make_transition(tr(2),1,3,0.1_dp)
  call pfinal_fmsm(tr,3,1,100.0_dp,pf)
  if(abs(pf(2)-0.8_dp)>2.0e-4_dp.or.abs(pf(3)-0.2_dp)>2.0e-4_dp) &
    error stop 'pfinal probabilities'
  call simfinal_fmsm(tr,3,1,50.0_dp,8000,[0.5_dp],fs,123)
  i2=findidx(fs%state,2);i3=findidx(fs%state,3)
  if(i2==0.or.i3==0)error stop 'absorbing states'
  if(abs(fs%probability(i2)-0.8_dp)>0.035_dp)error stop 'simfinal p2'
  if(abs(fs%probability(i3)-0.2_dp)>0.035_dp)error stop 'simfinal p3'
  if(abs(fs%mean_time(i2)-2.0_dp)>0.15_dp)error stop 'simfinal mean'
  if(abs(fs%quantile(i2,1)-log(2.0_dp)/0.5_dp)>0.15_dp)error stop 'simfinal median'
  call simfinal_fmsm_ci(tr,3,1,50.0_dp,800,[0.5_dp],12,0.95_dp,fs, &
    plo,phi,mlo,mhi,qlo,qhi,456)
  if(plo(i2)>fs%probability(i2).or.phi(i2)<fs%probability(i2)) &
    error stop 'simfinal ci'
  print *, 'test_final_states: PASS'
contains
  subroutine make_transition(t,from,to,rate)
    type(flexsurv_transition),intent(out)::t
    integer,intent(in)::from,to
    real(dp),intent(in)::rate
    t%from=from;t%to=to;t%model_kind=msm_parametric;t%row=1
    call initialize_spec(t%spec,dist_exponential,1,[rate])
    t%theta=initial_theta(t%spec)
    allocate(t%covariance(1,1));t%covariance(1,1)=0.0025_dp
  end subroutine make_transition
  integer function findidx(x,v) result(k)
    integer,intent(in)::x(:),v
    integer::j
    k=0;do j=1,size(x);if(x(j)==v)then;k=j;return;end if;end do
  end function findidx
end program test_final_states
