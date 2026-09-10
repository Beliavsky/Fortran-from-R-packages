! SPDX-License-Identifier: GPL-2.0-only
program high_level_model
   use marss_api
   implicit none

   real(dp) :: y(2, 8)
   type(marss_model_spec) :: spec
   type(marss_constraints) :: constraints
   type(marss_fit_result) :: fit
   integer :: info

   y(1, :) = [0.2_dp, 0.4_dp, 0.5_dp, 0.7_dp, 0.8_dp, 1.0_dp, 1.1_dp, 1.3_dp]
   y(2, :) = [-0.1_dp, 0.0_dp, 0.2_dp, 0.1_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.7_dp]

   spec%b = "diagonal and equal"
   spec%a = "zero"
   spec%v0 = "zero"

   call marss_from_data(y, spec, 8, 1.0e-5_dp, fit, constraints, "kem", info)
   if (info /= 0) error stop "high-level MARSS fit failed"

   write (*, '(a,f12.6)') "logLik = ", fit%loglik
   write (*, '(a,i0)') "free parameters = ", size(fit%free_parameters)
end program high_level_model
