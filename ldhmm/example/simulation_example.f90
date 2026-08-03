! SPDX-License-Identifier: Artistic-2.0
program simulation_example
   use ldhmm
   implicit none
   type(ldhmm_model) :: model, simulated
   real(dp) :: param(2,3), gamma_matrix(2,2), delta(2)
   real(dp), allocatable :: acf(:)
   integer :: status

   param(1,:) = [0.002_dp,0.015_dp,1.0_dp]
   param(2,:) = [-0.003_dp,0.040_dp,1.4_dp]
   gamma_matrix(1,:) = [0.97_dp,0.03_dp]
   gamma_matrix(2,:) = [0.08_dp,0.92_dp]
   delta = [0.70_dp,0.30_dp]
   model = ldhmm_create(2,param,gamma_matrix,delta)

   call seed_random(20260801)
   simulated = ldhmm_simulate_state_transition(model,init=20,status=status)
   write(*,'(a)') 'states:'
   write(*,'(*(i0,1x))') simulated%states_local
   write(*,'(a)') 'observations:'
   write(*,'(*(es11.3,1x))') simulated%observations

   acf = ldhmm_simulate_abs_acf(model,n=5000,lag_max=5,status=status)
   write(*,'(a,*(f9.5,1x))') 'simulated absolute-return ACF: ', acf
end program simulation_example
