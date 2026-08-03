program test_lp
   use highs
   implicit none
   real(dp) :: a(3,2), c(2), lo(2), up(2), lhs(3), rhs(3)
   type(highs_solution) :: sol
   integer(highs_int) :: status
   if (.not. highs_backend_available()) then
      print *, "test_lp: SKIP (run scripts/build_backend first)"
      stop
   end if
   c = [1.0_dp, 1.0_dp]
   lo = [0.0_dp, 1.0_dp]
   up = [4.0_dp, highs_default_infinity]
   a = reshape([0.0_dp,1.0_dp,3.0_dp, 1.0_dp,2.0_dp,2.0_dp], [3,2])
   lhs = [-highs_default_infinity, 5.0_dp, 6.0_dp]
   rhs = [7.0_dp, 15.0_dp, highs_default_infinity]
   call highs_solve(c, lo, up, sol, status, a=a, lhs=lhs, rhs=rhs, offset=3.0_dp)
   if (status /= highs_status_ok) error stop "LP solve call failed"
   if (sol%model_status /= highs_model_optimal) error stop "LP not optimal"
   if (maxval(abs(sol%col_value - [0.5_dp,2.25_dp])) > 1.0e-7_dp) error stop "wrong LP solution"
   if (abs(sol%objective_value - 5.75_dp) > 1.0e-7_dp) error stop "wrong LP objective"
   print *, "test_lp: PASS"
end program test_lp
