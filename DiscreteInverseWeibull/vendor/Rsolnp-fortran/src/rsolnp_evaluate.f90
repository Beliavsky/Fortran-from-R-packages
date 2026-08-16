! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_evaluate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rsolnp_kinds, only : dp
   use rsolnp_types, only : solnp_problem
   implicit none
   private

   public :: evaluate_objective, evaluate_gradient
   public :: evaluate_equalities, evaluate_eq_jacobian
   public :: evaluate_inequalities, evaluate_ineq_jacobian
   public :: finite_difference_gradient, finite_difference_jacobian
   public :: constraint_values, constraint_jacobian

contains

   subroutine evaluate_objective(problem, x, value, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      logical, intent(out) :: ok

      if (.not. associated(problem%fn)) then
         value = huge(1.0_dp)
         ok = .false.
         return
      end if
      if (allocated(problem%data)) then
         call problem%fn(x, value, problem%data)
      else
         call problem%fn(x, value)
      end if
      ok = ieee_is_finite(value)
   end subroutine evaluate_objective

   subroutine evaluate_gradient(problem, x, gradient, delta, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      real(dp), intent(in) :: delta
      logical, intent(out) :: ok

      if (associated(problem%gr)) then
         if (allocated(problem%data)) then
            call problem%gr(x, gradient, problem%data)
         else
            call problem%gr(x, gradient)
         end if
         ok = all(ieee_is_finite(gradient))
      else
         call finite_difference_gradient(problem, x, gradient, delta, ok)
      end if
   end subroutine evaluate_gradient

   subroutine finite_difference_gradient(problem, x, gradient, delta, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      real(dp), intent(in) :: delta
      logical, intent(out) :: ok

      real(dp), allocatable :: xp(:), xm(:)
      real(dp) :: fp, fm, h
      logical :: okp, okm
      integer :: i

      allocate(xp(size(x)), xm(size(x)))
      ok = .true.
      do i = 1, size(x)
         h = delta * max(1.0_dp, abs(x(i)))
         if (h <= 0.0_dp) h = sqrt(epsilon(1.0_dp))
         xp = x
         xm = x
         xp(i) = x(i) + h
         xm(i) = x(i) - h
         call evaluate_objective(problem, xp, fp, okp)
         call evaluate_objective(problem, xm, fm, okm)
         if (.not. okp .or. .not. okm) then
            gradient = 0.0_dp
            ok = .false.
            return
         end if
         gradient(i) = (fp - fm) / (2.0_dp * h)
      end do
   end subroutine finite_difference_gradient

   subroutine raw_equalities(problem, x, value, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      logical, intent(out) :: ok

      if (size(value) == 0) then
         ok = .true.
         return
      end if
      if (.not. associated(problem%eq_fn)) then
         value = 0.0_dp
         ok = .false.
         return
      end if
      if (allocated(problem%data)) then
         call problem%eq_fn(x, value, problem%data)
      else
         call problem%eq_fn(x, value)
      end if
      ok = all(ieee_is_finite(value))
   end subroutine raw_equalities

   subroutine raw_inequalities(problem, x, value, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      logical, intent(out) :: ok

      if (size(value) == 0) then
         ok = .true.
         return
      end if
      if (.not. associated(problem%ineq_fn)) then
         value = 0.0_dp
         ok = .false.
         return
      end if
      if (allocated(problem%data)) then
         call problem%ineq_fn(x, value, problem%data)
      else
         call problem%ineq_fn(x, value)
      end if
      ok = all(ieee_is_finite(value))
   end subroutine raw_inequalities

   subroutine evaluate_equalities(problem, x, value, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: raw(:)

      if (problem%n_eq == 0) then
         ok = .true.
         return
      end if
      allocate(raw(problem%raw_n_eq))
      call raw_equalities(problem, x, raw, ok)
      if (.not. ok) return
      if (problem%standard_form) then
         value = raw - problem%standard_eq_shift
      else
         value = raw - problem%eq_b
      end if
   end subroutine evaluate_equalities

   subroutine evaluate_inequalities(problem, x, value, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: raw(:)
      integer :: i, k

      if (problem%n_ineq == 0) then
         ok = .true.
         return
      end if
      allocate(raw(problem%raw_n_ineq))
      call raw_inequalities(problem, x, raw, ok)
      if (.not. ok) return
      if (.not. problem%standard_form) then
         value = raw
         return
      end if
      k = 0
      do i = 1, problem%raw_n_ineq
         if (ieee_is_finite(problem%standard_ineq_lower(i))) then
            k = k + 1
            value(k) = problem%standard_ineq_lower(i) - raw(i)
         end if
         if (ieee_is_finite(problem%standard_ineq_upper(i))) then
            k = k + 1
            value(k) = raw(i) - problem%standard_ineq_upper(i)
         end if
      end do
   end subroutine evaluate_inequalities

   subroutine evaluate_eq_jacobian(problem, x, jacobian, delta, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: jacobian(:, :)
      real(dp), intent(in) :: delta
      logical, intent(out) :: ok
      real(dp), allocatable :: raw(:, :)

      if (problem%n_eq == 0) then
         ok = .true.
         return
      end if
      allocate(raw(problem%raw_n_eq, problem%n))
      if (associated(problem%eq_jac)) then
         if (allocated(problem%data)) then
            call problem%eq_jac(x, raw, problem%data)
         else
            call problem%eq_jac(x, raw)
         end if
         ok = all(ieee_is_finite(raw))
      else
         call finite_difference_jacobian(problem, x, raw, delta, .true., ok)
      end if
      if (ok) jacobian = raw
   end subroutine evaluate_eq_jacobian

   subroutine evaluate_ineq_jacobian(problem, x, jacobian, delta, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: jacobian(:, :)
      real(dp), intent(in) :: delta
      logical, intent(out) :: ok
      real(dp), allocatable :: raw(:, :)
      integer :: i, k

      if (problem%n_ineq == 0) then
         ok = .true.
         return
      end if
      allocate(raw(problem%raw_n_ineq, problem%n))
      if (associated(problem%ineq_jac)) then
         if (allocated(problem%data)) then
            call problem%ineq_jac(x, raw, problem%data)
         else
            call problem%ineq_jac(x, raw)
         end if
         ok = all(ieee_is_finite(raw))
      else
         call finite_difference_jacobian(problem, x, raw, delta, .false., ok)
      end if
      if (.not. ok) return
      if (.not. problem%standard_form) then
         jacobian = raw
         return
      end if
      k = 0
      do i = 1, problem%raw_n_ineq
         if (ieee_is_finite(problem%standard_ineq_lower(i))) then
            k = k + 1
            jacobian(k, :) = -raw(i, :)
         end if
         if (ieee_is_finite(problem%standard_ineq_upper(i))) then
            k = k + 1
            jacobian(k, :) = raw(i, :)
         end if
      end do
   end subroutine evaluate_ineq_jacobian

   subroutine finite_difference_jacobian(problem, x, jacobian, delta, equality, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: jacobian(:, :)
      real(dp), intent(in) :: delta
      logical, intent(in) :: equality
      logical, intent(out) :: ok

      real(dp), allocatable :: xp(:), xm(:), fp(:), fm(:)
      real(dp) :: h
      logical :: okp, okm
      integer :: i, m

      m = size(jacobian, 1)
      allocate(xp(size(x)), xm(size(x)), fp(m), fm(m))
      ok = .true.
      do i = 1, size(x)
         h = delta * max(1.0_dp, abs(x(i)))
         if (h <= 0.0_dp) h = sqrt(epsilon(1.0_dp))
         xp = x
         xm = x
         xp(i) = x(i) + h
         xm(i) = x(i) - h
         if (equality) then
            call raw_equalities(problem, xp, fp, okp)
            call raw_equalities(problem, xm, fm, okm)
         else
            call raw_inequalities(problem, xp, fp, okp)
            call raw_inequalities(problem, xm, fm, okm)
         end if
         if (.not. okp .or. .not. okm) then
            jacobian = 0.0_dp
            ok = .false.
            return
         end if
         jacobian(:, i) = (fp - fm) / (2.0_dp * h)
      end do
   end subroutine finite_difference_jacobian

   subroutine constraint_values(problem, x, slack, constraints, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: slack(:)
      real(dp), intent(out) :: constraints(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: eq(:), ineq(:)
      logical :: ok1, ok2

      allocate(eq(problem%n_eq), ineq(problem%n_ineq))
      call evaluate_equalities(problem, x, eq, ok1)
      call evaluate_inequalities(problem, x, ineq, ok2)
      ok = ok1 .and. ok2
      if (.not. ok) then
         constraints = 0.0_dp
         return
      end if
      if (problem%n_eq > 0) constraints(1:problem%n_eq) = eq
      if (problem%n_ineq > 0) then
         constraints(problem%n_eq + 1:) = ineq - slack
      end if
   end subroutine constraint_values

   subroutine constraint_jacobian(problem, x, jacobian, delta, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: jacobian(:, :)
      real(dp), intent(in) :: delta
      logical, intent(out) :: ok
      real(dp), allocatable :: eq_j(:, :), ineq_j(:, :)
      logical :: ok1, ok2

      allocate(eq_j(problem%n_eq, problem%n), ineq_j(problem%n_ineq, problem%n))
      call evaluate_eq_jacobian(problem, x, eq_j, delta, ok1)
      call evaluate_ineq_jacobian(problem, x, ineq_j, delta, ok2)
      ok = ok1 .and. ok2
      if (.not. ok) then
         jacobian = 0.0_dp
         return
      end if
      if (problem%n_eq > 0) jacobian(1:problem%n_eq, :) = eq_j
      if (problem%n_ineq > 0) jacobian(problem%n_eq + 1:, :) = ineq_j
   end subroutine constraint_jacobian

end module rsolnp_evaluate
