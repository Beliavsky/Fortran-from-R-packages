! SPDX-License-Identifier: GPL-2.0-or-later
program regression_mixture_example
  use mixtools
  implicit none
  type(rng_state) :: rng
  type(regression_mixture_result) :: fit
  type(em_control) :: control
  real(dp), allocatable :: x(:,:),y(:)
  integer :: i,n

  n=300;allocate(x(n,1),y(n));call rng_seed(rng,91)
  do i=1,n
    x(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
    if(i<=n/2)then
      y(i)=-1.0_dp+0.5_dp*x(i,1)+0.25_dp*random_normal(rng)
    else
      y(i)=2.0_dp-0.8_dp*x(i,1)+0.30_dp*random_normal(rng)
    end if
  end do
  control%max_iterations=1000
  call regmixEM(y,x,2,fit,control,.true.)
  print '(a,*(f10.5,1x))','weights: ',fit%lambda
  print '(a)','component coefficients (intercept, slope):'
  do i=1,2
    print '(i3,2x,*(f10.5,1x))',i,fit%beta(:,i)
  end do
  print '(a,*(f10.5,1x))','scales: ',fit%sigma
end program regression_mixture_example
