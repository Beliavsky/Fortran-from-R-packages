program lp_example
   use highs
   implicit none
   real(dp) :: a(3,2)
   type(highs_solution) :: solution
   integer(highs_int) :: status

   if (.not. highs_backend_available()) then
      print '(a)', "HiGHS backend not found. Run scripts/build_backend.bat (Windows) or scripts/build_backend.sh."
      stop
   end if

   a = reshape([0.0_dp,1.0_dp,3.0_dp, 1.0_dp,2.0_dp,2.0_dp], [3,2])
   call highs_solve([1.0_dp,1.0_dp], [0.0_dp,1.0_dp], &
      [4.0_dp,highs_default_infinity], solution, status, a=a, &
      lhs=[-highs_default_infinity,5.0_dp,6.0_dp], &
      rhs=[7.0_dp,15.0_dp,highs_default_infinity], offset=3.0_dp)
   print '(a,i0)', "model status: ", solution%model_status
   print '(a,a)', "message: ", trim(solution%status_message)
   print '(a,*(f10.5,1x))', "x: ", solution%col_value
   print '(a,f12.6)', "objective: ", solution%objective_value
end program lp_example
