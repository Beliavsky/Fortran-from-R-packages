! Public facade for the modern Fortran Rdsdp/DSDP translation.
! DSDP copyright/license notice retained in LICENSE and licenses/DSDP-LICENSE.
module rdsdp
   use rdsdp_kinds, only : dp
   use rdsdp_types
   use rdsdp_problem_mod, only : dsdp_from_sedumi, validate_problem, symmetrize_problem
   use rdsdp_solver, only : dsdp_solve
   use rdsdp_io, only : read_sdpa, write_sdpa, flatten_primal
   use rdsdp_options, only : read_dsdp_options
   use rdsdp_data, only : get_data_dense, set_constraint_lowrank, set_objective_lowrank, data_nnz
   implicit none
   public
contains

   subroutine dsdp(amat,b,c,l,s,solution,control)
      ! Fortran analogue of Rdsdp::dsdp(A,b,C,K,OPTIONS).
      real(dp), intent(in) :: amat(:,:), b(:), c(:)
      integer, intent(in) :: l, s(:)
      type(dsdp_solution), intent(out) :: solution
      type(dsdp_control), intent(in), optional :: control
      type(dsdp_problem) :: problem
      call dsdp_from_sedumi(amat,b,c,l,s,problem)
      if (present(control)) then
         call dsdp_solve(problem,solution,control)
      else
         call dsdp_solve(problem,solution)
      end if
   end subroutine dsdp

   subroutine dsdp_readsdpa(filename,solution,options_filename,control)
      ! Fortran analogue of Rdsdp::dsdp.readsdpa().
      character(len=*), intent(in) :: filename
      type(dsdp_solution), intent(out) :: solution
      character(len=*), intent(in), optional :: options_filename
      type(dsdp_control), intent(in), optional :: control
      type(dsdp_problem) :: problem
      type(dsdp_control) :: ctrl
      ctrl=dsdp_control()
      if (present(control)) ctrl=control
      if (present(options_filename)) then
         if (len_trim(options_filename)>0) call read_dsdp_options(trim(options_filename),ctrl)
      end if
      call read_sdpa(filename,problem)
      call dsdp_solve(problem,solution,ctrl)
   end subroutine dsdp_readsdpa

end module rdsdp
