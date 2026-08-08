module example_stosoo_problem
   use oor, only : dp, guirland
   implicit none
contains
   function objective(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = -guirland(x(1))
   end function objective
end module example_stosoo_problem

program stosoo_guirland
   use oor, only : dp, stosoo, stosoo_options, stosoo_result
   use example_stosoo_problem, only : objective
   implicit none
   type(stosoo_options) :: options
   type(stosoo_result) :: result
   real(dp) :: lower(1), upper(1)

   lower = 0.0_dp
   upper = 1.0_dp
   options%stochastic = .false.
   call stosoo(objective, lower, upper, 500, result, options)
   print '(a,f12.8)', "x = ", result%par(1)
   print '(a,es16.8)', "f(x) = ", result%value
end program stosoo_guirland
