module example_poo_problem
   use oor, only : dp
   implicit none
contains
   function objective(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = -(x - 0.35_dp)**2
   end function objective
end module example_poo_problem

program poo_quadratic
   use oor, only : dp, poo, poo_result, set_random_seed
   use example_poo_problem, only : objective
   implicit none
   type(poo_result) :: result

   call set_random_seed(12)
   call poo(objective, 300, 0.0_dp, result, rhomax=10, nu=1.0_dp)
   print '(a,f12.8)', "x = ", result%par
   print '(a,es16.8)', "f(x) = ", result%value
   print '(a,f12.8)', "best rho = ", result%best_rho
end program poo_quadratic
