! SPDX-License-Identifier: GPL-2.0-or-later
program test_regression
  use mixtools
  implicit none
  type(em_control) :: ctl
  type(regression_mixture_result) :: fit, lfit, pfit, hfit, mfit
  type(rng_state) :: rng
  real(dp), allocatable :: x(:,:), y(:), yb(:), yp(:), re(:,:)
  integer, allocatable :: groups(:)
  integer :: i, n

  ctl%max_iterations=500;ctl%tolerance=1.0e-7_dp;call rng_seed(rng,24680);n=240
  allocate(x(n,1),y(n),yb(n),yp(n),groups(n));do i=1,n;x(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
    if(i<=n/2)then;y(i)=-1.0_dp+0.5_dp*x(i,1)+0.25_dp*random_normal(rng)
    else;y(i)=2.0_dp-0.7_dp*x(i,1)+0.3_dp*random_normal(rng);end if
    yb(i)=real(random_binomial(rng,1,logistic(-0.3_dp+1.2_dp*x(i,1))),dp)
    yp(i)=real(random_poisson(rng,exp(0.2_dp+0.4_dp*x(i,1))),dp);groups(i)=1+mod(i-1,12)
  end do
  call regmixEM(y,x,2,fit,ctl,.true.)
  call check(fit%status==0.and.size(fit%beta,2)==2,"linear regression mixture")
  call logisregmixEM(yb,x,1,lfit,control=ctl,addintercept=.true.)
  call check(lfit%status==0.and.lfit%beta(2,1)>0.4_dp,"logistic regression mixture")
  call poisregmixEM(yp,x,1,pfit,ctl,.true.)
  call check(pfit%status==0.and.pfit%beta(2,1)>0.1_dp,"Poisson regression mixture")
  call hmeEM(y,x,hfit,ctl,.true.)
  call check(size(hfit%auxiliary,1)==2,"mixture of experts")
  call regmixEM_mixed(y,x,groups,2,mfit,re,ctl,.true.)
  call check((mfit%status==0.or.mfit%status==MIXTOOLS_NOT_CONVERGED).and.size(re,1)==12,"mixed regression mixture")
  print '(a)', 'test_regression: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine check
end program test_regression
