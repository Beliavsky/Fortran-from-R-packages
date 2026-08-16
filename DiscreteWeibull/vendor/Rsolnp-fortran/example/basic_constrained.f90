! SPDX-License-Identifier: GPL-2.0-only
module basic_constrained_callbacks
   use rsolnp, only : dp
   implicit none
contains
   subroutine objective(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      if (present(data)) continue
      value = (x(1) - 2.0_dp) ** 2 + (x(2) - 1.0_dp) ** 2
   end subroutine objective

   subroutine gradient(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      class(*), intent(in), optional :: data
      if (present(data)) continue
      value = [2.0_dp * (x(1) - 2.0_dp), 2.0_dp * (x(2) - 1.0_dp)]
   end subroutine gradient

   subroutine equality(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      class(*), intent(in), optional :: data
      if (present(data)) continue
      value(1) = x(1) + x(2)
   end subroutine equality

   subroutine equality_jacobian(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:, :)
      class(*), intent(in), optional :: data
      if (present(data)) continue
      if (size(x) < 0) continue
      value(1, :) = [1.0_dp, 1.0_dp]
   end subroutine equality_jacobian
end module basic_constrained_callbacks

program basic_constrained
   use rsolnp, only : dp, solnp_problem, solnp_result, solnp_control, csolnp
   use basic_constrained_callbacks, only : objective, gradient, equality, equality_jacobian
   implicit none

   type(solnp_problem) :: problem
   type(solnp_result) :: result
   type(solnp_control) :: control

   problem%name = 'quadratic with equality'
   problem%n = 2
   problem%n_eq = 1
   problem%fn => objective
   problem%gr => gradient
   problem%eq_fn => equality
   problem%eq_jac => equality_jacobian
   allocate(problem%start(2), problem%lower(2), problem%upper(2), problem%eq_b(1))
   problem%start = [0.0_dp, 0.0_dp]
   problem%lower = -10.0_dp
   problem%upper = 10.0_dp
   problem%eq_b = 1.0_dp
   control%trace = 1

   call csolnp(problem, result, control)
   write(*, '(a,i0)') 'convergence: ', result%convergence
   write(*, '(a,es16.8)') 'objective:   ', result%objective
   write(*, '(a,2f14.8)') 'parameters:  ', result%pars
   write(*, '(a,es12.4)') 'eq violation:', result%kkt%eq_violation
end program basic_constrained
