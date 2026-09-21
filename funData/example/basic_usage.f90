program basic_usage
   use funData_api
   implicit none

   type(fun_data) :: object
   real(dp), allocatable :: integrals(:)
   real(dp) :: x(2, 5)
   logical :: ok

   x(1, :) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   x(2, :) = [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
   call create_fun_data_1d([0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp], x, object, ok)
   if (.not. ok) error stop "could not construct functional data"
   call integrate_fun_data(object, integrals, ok = ok)
   if (.not. ok) error stop "could not integrate functional data"
   print '(a,2f10.4)', "integrals:", integrals
end program basic_usage
