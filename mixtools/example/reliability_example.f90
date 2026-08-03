! SPDX-License-Identifier: GPL-2.0-or-later
program reliability_example
  use mixtools
  implicit none
  type(rng_state) :: rng
  type(reliability_mixture_result) :: fit
  type(em_control) :: control
  real(dp), allocatable :: time(:)
  integer, allocatable :: event(:)
  integer :: i,n

  n=300;allocate(time(n),event(n));call rng_seed(rng,731)
  do i=1,n/2
    time(i)=random_exponential(rng,1.5_dp);event(i)=1
  end do
  do i=n/2+1,n
    time(i)=random_weibull(rng,2.2_dp,2.0_dp)
    event(i)=merge(1,0,time(i)<2.5_dp);time(i)=min(time(i),2.5_dp)
  end do
  control%max_iterations=500
  call weibullRMM_SEM(time,event,2,fit,control)
  print '(a,*(f10.5,1x))','weights: ',fit%lambda
  print '(a,*(f10.5,1x))','shapes:  ',fit%shape
  print '(a,*(f10.5,1x))','scales:  ',fit%scale
end program reliability_example
