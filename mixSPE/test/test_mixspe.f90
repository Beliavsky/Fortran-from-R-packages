program test_mixspe
  use mixspe, only : dp, dpe, log_dpe, rpe, dspe, log_dspe, rspe, cov_pe, spe_model, em_fit
  implicit none
  real(dp), allocatable :: x(:,:), x1(:,:), x2(:,:), cov(:,:)
  real(dp) :: mu(2), s(2,2), psi(2), m1(2)
  type(spe_model) :: fit
  integer :: correct, map1, map2
  mu=0.0_dp; s=0.0_dp; s(1,1)=1.0_dp; s(2,2)=1.0_dp; psi=0.0_dp
  call check(abs(dpe(mu,mu,s,1.0_dp)-1.0_dp/(2.0_dp*acos(-1.0_dp)))<1.0e-12_dp,'normal density')
  call check(abs(dspe(mu,mu,s,psi,0.8_dp)-dpe(mu,mu,s,0.8_dp))<1.0e-12_dp,'zero-skew identity')
  x=rpe(30000,1.0_dp,mu,s)
  m1=sum(x,dim=1)/real(size(x,1),dp)
  call check(maxval(abs(m1))<0.04_dp,'rpe mean')
  cov=cov_pe(s,1.0_dp)
  call check(maxval(abs(cov-s))<1.0e-12_dp,'normal covariance identity')
  x1=rpe(250,2.0_dp,[0.0_dp,0.0_dp],s)
  x2=rpe(250,5.0_dp,[3.0_dp,0.0_dp],s)
  deallocate(x)
  allocate(x(500,2)); x(1:250,:)=x1; x(251:500,:)=x2
  call em_fit(x,2,'EIIV',fit,max_iter=150,tol=1.0e-5_dp)
  call check(fit%loglik > -huge(1.0_dp)/10.0_dp,'finite mixture loglik')
  call check(all(fit%pi>0.05_dp) .and. all(fit%pi<0.95_dp),'mixture proportions')
  map1=nint(sum(real(fit%map(1:250),dp))/250.0_dp)
  map2=nint(sum(real(fit%map(251:500),dp))/250.0_dp)
  correct=count(fit%map(1:250)==map1)+count(fit%map(251:500)==map2)
  call check(map1/=map2 .and. correct>400,'two-group clustering')
  x=rspe(250,[0.0_dp,0.0_dp],s,1.0_dp,[2.0_dp,0.0_dp],burnin=200)
  call check(all(abs(sum(x*x,dim=2))<huge(1.0_dp)),'rspe finite')
  print '(a)', 'test_mixspe: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print '(a)', 'FAIL: '//trim(msg)
      error stop 1
    end if
  end subroutine
end program
