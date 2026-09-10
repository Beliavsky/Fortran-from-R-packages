! SPDX-License-Identifier: GPL-2.0-only
program dfa_high_level
   use marss_api
   implicit none

   real(dp) :: y(2, 8)
   type(marss_dfa_spec) :: spec
   type(marss_constraints) :: constraints
   type(marss_fit_result) :: fit
   integer :: info

   y(1, :) = [1.0_dp, 1.2_dp, 1.4_dp, 1.5_dp, 1.7_dp, 1.8_dp, 2.0_dp, 2.1_dp]
   y(2, :) = [2.0_dp, 1.9_dp, 2.1_dp, 2.3_dp, 2.2_dp, 2.4_dp, 2.5_dp, 2.7_dp]

   spec%ntrends = 1
   spec%demean = .true.
   spec%z_score = .true.
   spec%r = "diagonal and equal"

   call marss_from_data(y, spec, 3, 1.0e-4_dp, fit, constraints, "kem", info)
   if (info /= 0) error stop "high-level DFA fit failed"

   write (*, '(a,f12.6)') "DFA logLik = ", fit%loglik
   write (*, '(a,i0)') "DFA free parameters = ", size(fit%free_parameters)
end program dfa_high_level
