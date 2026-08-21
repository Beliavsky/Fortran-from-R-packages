program test_multistate
  use flexsurv
  implicit none
  type(flexsurv_transition)::tr(1)
  real(dp)::times(3),elos(2)
  real(dp),allocatable::p(:,:,:)
  integer::fails
  fails=0;times=[0.0_dp,1.0_dp,2.0_dp]
  tr(1)%from=1;tr(1)%to=2;tr(1)%model_kind=msm_parametric
  call initialize_spec(tr(1)%spec,dist_exponential,1,[0.5_dp]);allocate(tr(1)%theta(1));tr(1)%theta=log(0.5_dp)
  call pmatrix_flexsurv(tr,2,times,p)
  if(abs(p(1,1,2)-exp(-0.5_dp))>3e-5_dp)then;print *,'p11 ',p(1,1,2);fails=fails+1;end if
  if(abs(p(1,2,2)-(1.0_dp-exp(-0.5_dp)))>3e-5_dp)fails=fails+1
  if(abs(p(2,2,3)-1.0_dp)>1e-10_dp)fails=fails+1
  call expected_length_of_stay(tr,2,1,0.0_dp,2.0_dp,elos,401)
  if(abs(elos(1)-2.0_dp*(1.0_dp-exp(-1.0_dp)))>2e-5_dp)then;print *,'elos ',elos;fails=fails+1;end if
  if(fails>0)error stop 1
  print *,'test_multistate: PASS'
end program test_multistate
