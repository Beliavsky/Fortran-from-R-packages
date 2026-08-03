! SPDX-License-Identifier: GPL-2.0-or-later
program demo_mixtools
  use mixtools
  implicit none
  type(rng_state) :: rng
  type(mixture_result) :: normal_fit
  type(model_selection_result) :: selection
  type(em_control) :: control
  real(dp), allocatable :: x(:)

  call rng_seed(rng,1234);allocate(x(400))
  call rnormmix(rng,400,[0.4_dp,0.6_dp],[-2.0_dp,2.0_dp],[0.5_dp,0.8_dp],x)
  control%max_iterations=800
  call normalmixEM(x,2,normal_fit,control)
  call normalmix_model_selection(x,[1,2,3],selection,control)
  print '(a,*(f9.4,1x))','estimated means: ',normal_fit%mu
  print '(a,i0)','BIC-selected components: ',selection%best_bic
end program demo_mixtools
