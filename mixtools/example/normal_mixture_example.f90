! SPDX-License-Identifier: GPL-2.0-or-later
program normal_mixture_example
  use mixtools
  implicit none
  type(rng_state) :: rng
  type(mixture_result) :: fit
  type(em_control) :: control
  real(dp), allocatable :: x(:)

  call rng_seed(rng,2026)
  allocate(x(500))
  call rnormmix(rng,500,[0.35_dp,0.65_dp],[-1.5_dp,2.0_dp],[0.6_dp,0.9_dp],x)
  control%max_iterations=1000
  call normalmixEM(x,2,fit,control)
  print '(a,*(f10.5,1x))','weights: ',fit%lambda
  print '(a,*(f10.5,1x))','means:   ',fit%mu
  print '(a,*(f10.5,1x))','scales:  ',fit%sigma
  print '(a,f14.6)','loglik:  ',fit%loglik
end program normal_mixture_example
