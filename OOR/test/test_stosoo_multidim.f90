module test_stosoo_multidim_callbacks
   use oor, only : dp
   implicit none
contains
   function quadratic(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = (x(1) - 0.2_dp)**2 + (x(2) - 0.7_dp)**2
   end function quadratic
end module test_stosoo_multidim_callbacks

program test_stosoo_multidim
   use oor, only : dp, stosoo, stosoo_options, stosoo_result
   use test_stosoo_multidim_callbacks, only : quadratic
   implicit none
   type(stosoo_options) :: options
   type(stosoo_result) :: result
   real(dp) :: lower(2), upper(2)

   lower = 0.0_dp
   upper = 1.0_dp
   options%stochastic = .false.
   options%keep_tree = .true.
   call stosoo(quadratic, lower, upper, 500, result, options)
   if (maxval(abs(result%par - [0.2_dp, 0.7_dp])) > 1.0e-3_dp) then
      error stop "multidimensional StoSOO parameter failure"
   end if
   if (result%value > 1.0e-6_dp) error stop "multidimensional StoSOO value failure"
   if (.not. allocated(result%tree)) error stop "keep_tree failure"
   if (size(result%xs, 1) /= 2) error stop "history dimension failure"
   print *, "test_stosoo_multidim: PASS"
end program test_stosoo_multidim
