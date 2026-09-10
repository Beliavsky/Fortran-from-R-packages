! SPDX-License-Identifier: GPL-2.0-only
program constraints_example
   use marss_api
   implicit none

   type(marss_model) :: model
   type(marss_constraints) :: constraints
   type(marss_fit_result) :: fit

   allocate(model%y(2, 8), model%b(2, 2), model%u(2), model%q(2, 2))
   allocate(model%z(2, 2), model%a(2), model%r(2, 2), model%x0(2), model%v0(2, 2))
   model%y(1, :) = [0.2_dp, 0.4_dp, 0.5_dp, 0.7_dp, 0.8_dp, 1.0_dp, 1.1_dp, 1.3_dp]
   model%y(2, :) = [-0.1_dp, 0.0_dp, 0.2_dp, 0.1_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.7_dp]
   model%b = 0.0_dp
   model%b(1, 1) = 0.4_dp
   model%b(2, 2) = 0.4_dp
   model%u = 0.0_dp
   model%q = reshape([0.15_dp, 0.0_dp, 0.0_dp, 0.15_dp], [2, 2])
   model%z = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
   model%a = 0.0_dp
   model%r = reshape([0.10_dp, 0.0_dp, 0.0_dp, 0.10_dp], [2, 2])
   model%x0 = 0.0_dp
   model%v0 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
   model%tinitx = 0

   allocate(constraints%b%fixed(4), constraints%b%design(4, 1), constraints%b%start(1))
   constraints%b%fixed = 0.0_dp
   constraints%b%design(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
   constraints%b%start = 0.4_dp

   call marss_optim_linear(model, constraints, 20, 1.0e-5_dp, fit)
   if (fit%info /= 0) error stop "constrained fit failed"
   write (*, '(a,f10.5)') 'shared B diagonal: ', fit%model%b(1, 1)
   write (*, '(a,f12.6)') 'log likelihood: ', fit%loglik
end program constraints_example
