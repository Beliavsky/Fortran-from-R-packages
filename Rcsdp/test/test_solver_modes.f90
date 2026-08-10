program test_solver_modes
   use rcsdp
   implicit none
   type(csdp_problem) :: p
   type(csdp_solution) :: snew, snofill, snoscale
   type(csdp_control) :: ctrl
   integer :: info

   call read_sdpa_sparse('data/theta1.dat-s',p,info)
   if (info/=0) error stop 'test_solver_modes: read failed'
   ctrl%printlevel=0
   ctrl%maxiter=50

   call csdp(p,snew,ctrl)
   if (snew%status/=csdp_success) error stop 'default v0.3 mode failed'

   ctrl%use_fill_products=.false.
   call csdp(p,snofill,ctrl)
   if (snofill%status/=csdp_success) error stop 'no-fill mode failed'

   ctrl%use_fill_products=.true.
   ctrl%use_schur_scaling=.false.
   call csdp(p,snoscale,ctrl)
   if (snoscale%status/=csdp_success) error stop 'unscaled Schur mode failed'

   if (abs(snew%pobj-snofill%pobj)>2.0e-7_dp) error stop 'fill mode objective mismatch'
   if (abs(snew%pobj-snoscale%pobj)>2.0e-7_dp) error stop 'Schur scaling objective mismatch'
   if (snew%fill_nnz<=0) error stop 'fill statistics missing'
   print *, 'test_solver_modes: PASS'
end program test_solver_modes
