program test_theta
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_problem, csdp_solution, csdp_control
   use rcsdp_io, only : read_sdpa_sparse
   use rcsdp_solver, only : csdp
   implicit none
   type(csdp_problem) :: p
   type(csdp_solution) :: s
   type(csdp_control) :: ctrl
   integer :: info
   call read_sdpa_sparse('data/theta1.dat-s',p,info)
   if (info/=0) error stop 'read_sdpa_sparse failed'
   ctrl%printlevel=0
   ctrl%maxiter=50
   call csdp(p,s,ctrl)
   write(*,'("status=",i0," pobj=",es16.8," dobj=",es16.8," pinf=",es12.4," dinf=",es12.4," gap=",es12.4)') &
      s%status,s%pobj,s%dobj,s%pinfeas,s%dinfeas,s%relgap
   if (s%status/=0) error stop 'theta did not converge'
   if (abs(s%pobj-23.0_dp)>2.0e-5_dp) error stop 'theta objective mismatch'
end program test_theta
