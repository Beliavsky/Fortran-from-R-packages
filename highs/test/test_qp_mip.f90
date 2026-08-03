program test_qp_mip
   use highs
   implicit none
   type(highs_solution) :: sol
   integer(highs_int) :: status
   real(dp) :: q(2,2), a(1,2)

   if (.not. highs_backend_available()) then
      print *, "test_qp_mip: SKIP (run scripts/build_backend first)"
      stop
   end if

   q = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
   call highs_solve([-1.0_dp,-2.0_dp], [0.0_dp,0.0_dp], &
      [highs_default_infinity,highs_default_infinity], sol, status, q=q)
   if (status /= highs_status_ok .or. sol%model_status /= highs_model_optimal) &
      error stop "QP failed"
   if (maxval(abs(sol%col_value - [1.0_dp,2.0_dp])) > 2.0e-6_dp) error stop "wrong QP solution"
   if (abs(sol%objective_value + 2.5_dp) > 2.0e-6_dp) error stop "wrong QP objective"

   a = reshape([2.0_dp,1.0_dp], [1,2])
   call highs_solve([3.0_dp,2.0_dp], [0.0_dp,0.0_dp], [1.0_dp,1.0_dp], &
      sol, status, a=a, lhs=[-highs_default_infinity], rhs=[2.0_dp], &
      vartype=[highs_var_integer,highs_var_integer], maximum=.true.)
   if (status /= highs_status_ok .or. sol%model_status /= highs_model_optimal) &
      error stop "MIP failed"
   if (abs(sol%objective_value - 3.0_dp) > 1.0e-7_dp) error stop "wrong MIP objective"
   if (maxval(abs(sol%col_value - [1.0_dp,0.0_dp])) > 1.0e-7_dp) error stop "wrong MIP solution"
   print *, "test_qp_mip: PASS"
end program test_qp_mip
