program test_shared_multistate
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_spec,flexsurv_result,initialize_spec
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_shared_multistate
  implicit none
  type(flexsurv_spec)::sp
  type(flexsurv_result)::fit
  type(flexsurv_shared_msm)::sm
  real(dp),allocatable::cv(:,:,:)
  real(dp)::h1,h2,want,t(1)
  integer::st
  call initialize_spec(sp,dist_exponential,2,[0.1_dp])
  deallocate(sp%reg(1)%x);allocate(sp%reg(1)%x(2,1));sp%reg(1)%x(:,1)=[0.0_dp,1.0_dp]
  fit%theta=[log(0.1_dp),0.3_dp]
  allocate(fit%covariance(2,2));fit%covariance=reshape([0.01_dp,0.005_dp,0.005_dp,0.04_dp],[2,2])
  call make_shared_msm(sm,3,sp,fit,[1,1],[2,3],[1,2],status=st)
  if(st/=0)error stop 'make shared'
  t=[2.0_dp];call shared_cumhaz_covariance(sm,t,cv,status=st)
  if(st/=0)error stop 'shared covariance status'
  h1=2.0_dp*exp(fit%theta(1));h2=2.0_dp*exp(sum(fit%theta))
  want=h1*h2*(fit%covariance(1,1)+fit%covariance(1,2))
  if(abs(cv(1,2,1)-want)>2.0e-7_dp)error stop 'shared cross covariance'
  if(cv(1,2,1)<=0.0_dp)error stop 'shared covariance nonzero'
  print *,'test_shared_multistate: PASS'
end program
