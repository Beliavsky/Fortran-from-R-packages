program test_edge_cases
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use quarks
   implicit none
   real(dp) :: x(5), loss(5), var(5)
   type(risk_result) :: tail
   type(coverage_result) :: coverage
   type(rollcast_result) :: final_forecast, invalid

   x = 0.01_dp
   tail = hs(x, 0.95_dp, method_plain)
   if (tail%status /= quarks_empty_tail) error stop 'empty tail not reported'
   if (.not. ieee_is_nan(tail%es)) error stop 'empty-tail ES should be NaN'

   loss = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
   var = 2.0_dp
   coverage = cvgtest(loss, var, 0.95_dp)
   if (coverage%status /= quarks_no_violations) error stop 'no-violation status'

   x = [0.01_dp, -0.01_dp, 0.02_dp, -0.02_dp, 0.005_dp]
   final_forecast = rollcast(x, p=0.80_dp, method=method_plain, nout=0, nwin=5)
   if (size(final_forecast%var) /= 1 .or. size(final_forecast%xout) /= 0) then
      error stop 'nout=0 behavior'
   end if

   invalid = rollcast(x, p=0.80_dp, method=method_plain, nout=3, nwin=4)
   if (invalid%status /= quarks_invalid_input) error stop 'invalid size not caught'
   print *, 'test_edge_cases: PASS'
end program test_edge_cases
