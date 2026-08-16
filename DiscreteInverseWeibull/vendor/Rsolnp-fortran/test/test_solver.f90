! SPDX-License-Identifier: GPL-2.0-only
module test_solver_callbacks
   use rsolnp_kinds, only : dp
   implicit none
contains
   subroutine rosenbrock(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      if (present(data)) continue
      value = 100.0_dp * (x(2) - x(1) ** 2) ** 2 + (1.0_dp - x(1)) ** 2
   end subroutine rosenbrock

   subroutine rosenbrock_gradient(x, gradient, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      class(*), intent(in), optional :: data
      if (present(data)) continue
      gradient(1) = -400.0_dp * x(1) * (x(2) - x(1) ** 2) - 2.0_dp * (1.0_dp - x(1))
      gradient(2) = 200.0_dp * (x(2) - x(1) ** 2)
   end subroutine rosenbrock_gradient
end module test_solver_callbacks

program test_solver
   use rsolnp, only : dp, solnp_problem, solnp_result, solnp_control, csolnp, solnp_success
   use test_solver_callbacks, only : rosenbrock, rosenbrock_gradient
   implicit none

   type(solnp_problem) :: problem
   type(solnp_result) :: result
   type(solnp_control) :: control

   problem%name = 'Rosenbrock'
   problem%n = 2
   problem%fn => rosenbrock
   problem%gr => rosenbrock_gradient
   allocate(problem%start(2), problem%lower(2), problem%upper(2))
   problem%start = [-1.2_dp, 1.0_dp]
   problem%lower = -5.0_dp
   problem%upper = 5.0_dp
   control%max_iter = 20
   control%min_iter = 500
   control%tol = 1.0e-10_dp

   call csolnp(problem, result, control)
   call check(result%convergence == solnp_success, 'analytic Rosenbrock convergence')
   call check(result%objective < 1.0e-10_dp, 'analytic Rosenbrock objective')
   call check(maxval(abs(result%pars - 1.0_dp)) < 1.0e-4_dp, 'analytic Rosenbrock parameters')
   call check(result%kkt%stationarity < 2.0e-5_dp, 'analytic Rosenbrock KKT')

   nullify(problem%gr)
   problem%start = [-1.2_dp, 1.0_dp]
   control%delta = 1.0e-6_dp
   call csolnp(problem, result, control)
   call check(result%objective < 1.0e-7_dp, 'finite-difference Rosenbrock objective')
   call check(maxval(abs(result%pars - 1.0_dp)) < 1.0e-3_dp, 'finite-difference parameters')

   print '(a)', 'test_solver: PASS'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check
end program test_solver
