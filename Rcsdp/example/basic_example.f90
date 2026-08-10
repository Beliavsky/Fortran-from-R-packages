program basic_example
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_problem, csdp_solution, csdp_control, csdp_matrix, csdp_diag
   use rcsdp_problem_mod, only : init_problem, set_c_matrix_block, set_c_diag_block, &
      set_a_matrix_block, set_a_diag_block
   use rcsdp_solver, only : csdp
   implicit none
   type(csdp_problem) :: p
   type(csdp_solution) :: s
   type(csdp_control) :: ctrl
   real(dp) :: c1(2,2), c2(3,3), a11(2,2), a22(3,3)

   call init_problem(p,[csdp_matrix,csdp_matrix,csdp_diag],[2,3,2],2)
   c1=reshape([2.0_dp,1.0_dp,1.0_dp,2.0_dp],[2,2])
   c2=reshape([3.0_dp,0.0_dp,1.0_dp,0.0_dp,2.0_dp,0.0_dp,1.0_dp,0.0_dp,3.0_dp],[3,3])
   call set_c_matrix_block(p,1,c1)
   call set_c_matrix_block(p,2,c2)
   call set_c_diag_block(p,3,[0.0_dp,0.0_dp])

   a11=reshape([3.0_dp,1.0_dp,1.0_dp,3.0_dp],[2,2])
   call set_a_matrix_block(p,1,1,a11)
   call set_a_matrix_block(p,1,2,0.0_dp*c2)
   call set_a_diag_block(p,1,3,[1.0_dp,0.0_dp])

   a22=reshape([3.0_dp,0.0_dp,1.0_dp,0.0_dp,4.0_dp,0.0_dp,1.0_dp,0.0_dp,5.0_dp],[3,3])
   call set_a_matrix_block(p,2,1,0.0_dp*c1)
   call set_a_matrix_block(p,2,2,a22)
   call set_a_diag_block(p,2,3,[0.0_dp,1.0_dp])
   p%b=[1.0_dp,2.0_dp]

   ctrl%printlevel=1
   call csdp(p,s,ctrl)
   write(*,'("status = ",i0)') s%status
   write(*,'("pobj   = ",es16.8)') s%pobj
   write(*,'("dobj   = ",es16.8)') s%dobj
   write(*,'("gap    = ",es16.8)') s%relgap
   write(*,'("pinf   = ",es16.8)') s%pinfeas
   write(*,'("dinf   = ",es16.8)') s%dinfeas
end program basic_example
