program test_reference_data
   use rdsdp
   implicit none
   call check('data/control1.dat-s',-17.7846267761_dp,5.0e-5_dp)
   call check('data/truss1.dat-s',   8.99999631296_dp,5.0e-5_dp)
   call check('data/vibra1.dat-s', -40.8190124025_dp,8.0e-5_dp)
   call check_mcp100
   print *, 'test_reference_data: PASS'
contains
   subroutine check(file,reference,tol)
      character(len=*), intent(in) :: file
      real(dp), intent(in) :: reference,tol
      type(dsdp_solution) :: sol
      type(dsdp_control) :: ctrl
      ctrl=dsdp_control(); ctrl%gaptol=2.0e-7_dp; ctrl%pinfeastol=1.0e-7_dp
      call dsdp_readsdpa(file,sol,control=ctrl)
      write(*,'(a,1x,es16.8,1x,a,es10.2)') trim(file),sol%dobj,'error=',abs(sol%dobj-reference)
      if (abs(sol%dobj-reference)>tol) error stop 'reference objective mismatch'
   end subroutine check
   subroutine check_mcp100
      type(dsdp_problem) :: p
      type(dsdp_solution) :: sol
      type(dsdp_control) :: ctrl
      call read_sdpa('data/mcp100.dat-s',p)
      if (p%m/=100 .or. size(p%block)/=1 .or. p%block(1)%n/=100) error stop 'mcp100 parse mismatch'
      ctrl=dsdp_control(); ctrl%gaptol=1.0e-7_dp; ctrl%pinfeastol=1.0e-7_dp
      call dsdp_solve(p,sol,ctrl)
      write(*,'(a,1x,es16.8,1x,a,i0)') 'data/mcp100.dat-s',sol%dobj,'sparse-pairs=',sol%sparse_pair_evals
      if (abs(sol%dobj-(-226.15735223_dp))>8.0e-5_dp) error stop 'mcp100 objective mismatch'
      if (sol%sparse_pair_evals<=0) error stop 'mcp100 did not use sparse Schur path'
   end subroutine check_mcp100
end program test_reference_data
