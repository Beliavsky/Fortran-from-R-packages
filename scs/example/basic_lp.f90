! SPDX-License-Identifier: GPL-3.0-only
program basic_lp
   use scs_kinds, only : dp
   use scs_types
   use scs_sparse, only : dense_to_csc
   use scs_solver, only : scs
   implicit none

   type(scs_data) :: data
   type(scs_cone) :: cone
   type(scs_settings) :: settings
   type(scs_solution) :: solution
   type(scs_info) :: info
   real(dp) :: a(2,1)

   ! Minimize x subject to x = 1 (represented by two equality rows).
   a(:,1) = [1.0_dp, 1.0_dp]
   data%m = 2
   data%n = 1
   call dense_to_csc(a, data%A)
   data%b = [1.0_dp, 1.0_dp]
   data%c = [1.0_dp]
   cone%z = 2

   settings%max_iters = 1000
   call scs(data, cone, settings, solution, info)

   write(*,'(a,a)') 'status: ', trim(info%status)
   write(*,'(a,es16.8)') 'x: ', solution%x(1)
   write(*,'(a,es16.8)') 'primal objective: ', info%pobj
end program basic_lp
