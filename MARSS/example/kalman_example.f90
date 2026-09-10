! SPDX-License-Identifier: GPL-2.0-only
program kalman_example
   use marss_api
   implicit none

   type(marss_model) :: model
   type(marss_kf_result) :: kf

   allocate(model%y(1, 8), model%b(1, 1), model%u(1), model%q(1, 1))
   allocate(model%z(1, 1), model%a(1), model%r(1, 1), model%x0(1), model%v0(1, 1))
   model%y(1, :) = [0.2_dp, 0.4_dp, 0.3_dp, 0.7_dp, 0.8_dp, 0.6_dp, 1.0_dp, 1.1_dp]
   model%b = 0.85_dp
   model%u = 0.05_dp
   model%q = 0.08_dp
   model%z = 1.0_dp
   model%a = 0.0_dp
   model%r = 0.12_dp
   model%x0 = 0.0_dp
   model%v0 = 0.5_dp
   model%tinitx = 0

   call marss_kfss(model, kf)
   if (.not. kf%ok) error stop 'Kalman filter failed'
   write (*, '(a,f12.6)') 'log likelihood: ', kf%loglik
   write (*, '(a,*(f9.4,1x))') 'smoothed state: ', kf%x_smooth(1, :)
end program kalman_example
