module dslnex_problem
   use nleqslv_fortran, only : dp
   implicit none
contains
   subroutine equations(x, f)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
      f(1) = x(1)**2 + x(2)**2 - 2.0_dp
      f(2) = exp(x(1)-1.0_dp) + x(2)**3 - 2.0_dp
   end subroutine equations
end module dslnex_problem

program dslnex_example
   use nleqslv_fortran
   use dslnex_problem
   implicit none
   type(nleq_options) :: opt
   type(nleq_result) :: sol
   real(dp) :: x0(2)

   x0 = [2.0_dp, 0.5_dp]
   opt = nleq_options()
   opt%method = NLEQ_BROYDEN
   opt%global = NLEQ_DBLDOG
   call solve_nleqslv(x0, equations, sol, opt)

   print '(a,2f16.10)', 'solution: ', sol%x
   print '(a,2es16.6)', 'f(x):     ', sol%fvec
   print '(a,i0,2a)', 'termcd: ', sol%termcd, '  ', trim(sol%message)
end program dslnex_example
