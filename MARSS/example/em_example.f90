! SPDX-License-Identifier: GPL-2.0-only
program em_example
   use marss_api
   implicit none

   type(marss_model) :: model
   type(marss_fit_result) :: fit

   allocate(model%y(1, 12), model%b(1, 1), model%u(1), model%q(1, 1))
   allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
   model%y(1, :) = [0.5_dp, 0.7_dp, 0.2_dp, 1.0_dp, 1.1_dp, 0.9_dp, &
      1.3_dp, 1.0_dp, 1.4_dp, 1.6_dp, 1.5_dp, 1.8_dp]
   model%b = 0.5_dp
   model%u = 0.0_dp
   model%q = 0.5_dp
   model%z = 1.0_dp
   model%a = 0.0_dp
   model%r = 0.5_dp
   model%x0 = 0.0_dp
   model%v0 = 1.0_dp
   model%tinitx = 0

   call marss_kem(model, 40, 1.0e-7_dp, fit, estimate_z=.false., estimate_a=.false., estimate_v0=.false.)
   if (fit%info /= 0) error stop 'EM fit failed'
   write (*, '(a,f12.6)') 'fit log likelihood: ', fit%loglik
   write (*, '(a,f10.5)') 'B: ', fit%model%b(1, 1)
   write (*, '(a,f10.5)') 'U: ', fit%model%u(1)
   write (*, '(a,f10.5)') 'Q: ', fit%model%q(1, 1)
   write (*, '(a,f10.5)') 'R: ', fit%model%r(1, 1)
end program em_example
