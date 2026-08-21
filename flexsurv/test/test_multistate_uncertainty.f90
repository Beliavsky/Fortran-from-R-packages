program test_multistate_uncertainty
  use flexsurv_kinds, only : dp
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_fit, only : initialize_spec
  use flexsurv_multistate, only : flexsurv_transition, msm_parametric
  use flexsurv_multistate_uncertainty, only : pmatrix_simfs, totlos_simfs, &
    pmatrix_simfs_ci, pmatrix_flexsurv_ci
  implicit none
  type(flexsurv_transition)::tr(1)
  real(dp),allocatable::p(:,:,:),est(:,:,:),lo(:,:,:),hi(:,:,:)
  real(dp)::times(2),elos(2),truth
  integer::fails
  fails=0
  tr(1)%from=1;tr(1)%to=2;tr(1)%model_kind=msm_parametric;tr(1)%row=1
  call initialize_spec(tr(1)%spec,dist_exponential,1)
  allocate(tr(1)%theta(1),tr(1)%covariance(1,1));tr(1)%theta=log(0.5_dp);tr(1)%covariance=0.02_dp
  times=[0.0_dp,2.0_dp]
  call pmatrix_simfs(tr,2,times,18000,p,seed=1234)
  truth=exp(-1.0_dp)
  if(abs(p(1,1,2)-truth)>0.025_dp)then;print *,'FAIL pmatrix sim',p(1,1,2),truth;fails=fails+1;end if
  if(abs(sum(p(1,:,2))-1.0_dp)>1.0e-12_dp)then;print *,'FAIL rowsum';fails=fails+1;end if
  call totlos_simfs(tr,2,1,2.0_dp,18000,elos,seed=4321)
  truth=(1.0_dp-exp(-1.0_dp))/0.5_dp
  if(abs(elos(1)-truth)>0.035_dp)then;print *,'FAIL los1',elos(1),truth;fails=fails+1;end if
  if(abs(sum(elos)-2.0_dp)>0.02_dp)then;print *,'FAIL los sum',elos;fails=fails+1;end if
  call pmatrix_flexsurv_ci(tr,2,times,35,0.90_dp,est,lo,hi,seed=77)
  if(any(lo>hi))then;print *,'FAIL ODE ci ordering';fails=fails+1;end if
  if(.not.(lo(1,1,2)<est(1,1,2).and.hi(1,1,2)>est(1,1,2)))then
    print *,'FAIL ODE ci coverage',lo(1,1,2),est(1,1,2),hi(1,1,2);fails=fails+1
  end if
  call pmatrix_simfs_ci(tr,2,[2.0_dp],2000,25,0.90_dp,est,lo,hi,seed=99)
  if(any(lo>hi))then;print *,'FAIL sim ci ordering';fails=fails+1;end if
  if(fails==0)then
    print *,'test_multistate_uncertainty: PASS'
  else
    print *,'test_multistate_uncertainty: FAIL',fails;error stop 1
  end if
end program test_multistate_uncertainty
