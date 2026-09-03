! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Deterministic bounded optimization helpers used by translated ape workflows.
! This module is translation infrastructure and is not copied from upstream ape.
module ape_optimize
   use r_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   type, abstract, public :: bounded_problem
   contains
      procedure(problem_value_interface), deferred :: value
      procedure(problem_gradient_interface), deferred :: gradient
   end type bounded_problem

   abstract interface
      function problem_value_interface(self, x) result(value)
         import :: bounded_problem, dp
         class(bounded_problem), intent(inout) :: self !! Objective object holding any data needed to evaluate the problem.
         real(dp), intent(in) :: x(:) !! Parameter vector at which the minimization objective is evaluated.
         real(dp) :: value
      end function problem_value_interface

      subroutine problem_gradient_interface(self, x, gradient)
         import :: bounded_problem, dp
         class(bounded_problem), intent(inout) :: self !! Objective object holding any data needed for derivatives.
         real(dp), intent(in) :: x(:) !! Parameter vector at which the objective gradient is evaluated.
         real(dp), intent(out) :: gradient(:) !! Gradient of the minimization objective, with the same size as `x`.
      end subroutine problem_gradient_interface
   end interface

   public :: bounded_bfgs
   public :: finite_difference_hessian

contains

   subroutine bounded_bfgs(problem, x, lower, upper, objective, info, iterations, max_iter, tolerance)
      !! Minimizes a differentiable objective subject to independent box constraints.
      class(bounded_problem), intent(inout) :: problem !! Objective and gradient implementation for the fitted model.
      real(dp), intent(inout) :: x(:) !! Initial parameters on input and final bounded estimate on output.
      real(dp), intent(in) :: lower(:) !! Inclusive lower bound for every parameter.
      real(dp), intent(in) :: upper(:) !! Inclusive upper bound for every parameter.
      real(dp), intent(out) :: objective !! Objective value at the returned parameter vector.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid input, nonfinite values, or line-search failure.
      integer, intent(out), optional :: iterations !! Number of accepted quasi-Newton iterations.
      integer, intent(in), optional :: max_iter !! Maximum accepted iterations; default 500.
      real(dp), intent(in), optional :: tolerance !! Projected-gradient and relative-objective tolerance; default `1e-8`.
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: gradient_new(:)
      real(dp), allocatable :: projected(:)
      real(dp), allocatable :: direction(:)
      real(dp), allocatable :: x_new(:)
      real(dp), allocatable :: x_trial(:)
      real(dp), allocatable :: step(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: hessian_inverse(:, :)
      real(dp), allocatable :: identity(:, :)
      real(dp), allocatable :: transform(:, :)
      real(dp), allocatable :: work(:, :)
      real(dp) :: alpha
      real(dp) :: armijo
      real(dp) :: curvature
      real(dp) :: f_new
      real(dp) :: f_old
      real(dp) :: f_trial
      real(dp) :: rho
      real(dp) :: tol
      integer :: i
      integer :: iteration
      integer :: limit
      integer :: ls
      integer :: n
      logical :: accepted

      info = 0
      if (present(iterations)) iterations = 0
      n = size(x)
      if (size(lower) /= n .or. size(upper) /= n) then
         objective = huge(1.0_dp)
         info = 1
         return
      end if
      if (any(lower > upper)) then
         objective = huge(1.0_dp)
         info = 2
         return
      end if
      if (n == 0) then
         objective = problem%value(x)
         if (.not. ieee_is_finite(objective)) info = 3
         return
      end if

      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      limit = 500
      if (present(max_iter)) limit = max(1, max_iter)
      allocate(gradient(n), gradient_new(n), projected(n), direction(n), x_new(n), x_trial(n), step(n), y(n))
      allocate(hessian_inverse(n, n), identity(n, n), transform(n, n), work(n, n))
      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      hessian_inverse = identity
      x = max(lower, min(upper, x))
      objective = problem%value(x)
      if (.not. ieee_is_finite(objective)) then
         info = 3
         return
      end if
      call problem%gradient(x, gradient)
      if (.not. all(ieee_is_finite(gradient))) then
         info = 4
         return
      end if

      do iteration = 1, limit
         call projected_gradient(x, lower, upper, gradient, tol, projected)
         if (maxval(abs(projected)) <= tol * (1.0_dp + abs(objective))) then
            if (present(iterations)) iterations = iteration - 1
            return
         end if

         direction = -matmul(hessian_inverse, projected)
         call restrict_direction(x, lower, upper, direction, tol)
         if (dot_product(gradient, direction) >= -sqrt(epsilon(1.0_dp)) * max(1.0_dp, norm2(direction))) then
            direction = -projected
            call restrict_direction(x, lower, upper, direction, tol)
            hessian_inverse = identity
         end if
         if (maxval(abs(direction)) <= tol) then
            if (present(iterations)) iterations = iteration - 1
            return
         end if

         f_old = objective
         alpha = 1.0_dp
         armijo = 1.0e-4_dp
         accepted = .false.
         f_new = huge(1.0_dp)
         do ls = 1, 50
            x_trial = max(lower, min(upper, x + alpha * direction))
            step = x_trial - x
            if (maxval(abs(step)) <= epsilon(1.0_dp) * (1.0_dp + maxval(abs(x)))) then
               alpha = 0.5_dp * alpha
               cycle
            end if
            f_trial = problem%value(x_trial)
            if (ieee_is_finite(f_trial)) then
               if (f_trial <= objective + armijo * dot_product(gradient, step)) then
                  if (.not. accepted .or. f_trial < f_new) then
                     accepted = .true.
                     f_new = f_trial
                     x_new = x_trial
                  end if
               end if
            end if
            alpha = 0.5_dp * alpha
         end do
         if (accepted) step = x_new - x
         if (.not. accepted) then
            info = 5
            if (present(iterations)) iterations = iteration - 1
            return
         end if

         call problem%gradient(x_new, gradient_new)
         if (.not. all(ieee_is_finite(gradient_new))) then
            info = 6
            if (present(iterations)) iterations = iteration - 1
            return
         end if
         y = gradient_new - gradient
         curvature = dot_product(y, step)
         if (curvature > sqrt(epsilon(1.0_dp)) * max(1.0_dp, norm2(y) * norm2(step))) then
            rho = 1.0_dp / curvature
            transform = identity - rho * outer_product(step, y)
            work = matmul(transform, hessian_inverse)
            hessian_inverse = matmul(work, transpose(transform)) + rho * outer_product(step, step)
            hessian_inverse = 0.5_dp * (hessian_inverse + transpose(hessian_inverse))
         else
            hessian_inverse = identity
         end if

         x = x_new
         objective = f_new
         gradient = gradient_new
         if (present(iterations)) iterations = iteration
         if (abs(f_old - objective) <= tol * (1.0_dp + abs(f_old))) then
            call projected_gradient(x, lower, upper, gradient, tol, projected)
            if (maxval(abs(projected)) <= sqrt(tol) * (1.0_dp + abs(objective))) return
         end if
      end do
      info = 7
   end subroutine bounded_bfgs

   subroutine finite_difference_hessian(problem, x, hessian, info, relative_step, lower, upper)
      !! Computes a symmetric finite-difference Hessian of an objective at `x`.
      class(bounded_problem), intent(inout) :: problem !! Objective whose local Hessian is required.
      real(dp), intent(in) :: x(:) !! Parameter vector around which second derivatives are evaluated.
      real(dp), allocatable, intent(out) :: hessian(:, :) !! Symmetric Hessian matrix with shape `(size(x), size(x))`.
      integer, intent(out) :: info !! Zero on success or nonzero if a finite-difference evaluation is nonfinite.
      real(dp), intent(in), optional :: relative_step !! Relative finite-difference scale; default `epsilon**0.25`.
      real(dp), intent(in), optional :: lower(:) !! Optional lower parameter bounds used to keep perturbations feasible.
      real(dp), intent(in), optional :: upper(:) !! Optional upper parameter bounds used to keep perturbations feasible.
      real(dp), allocatable :: xp(:)
      real(dp), allocatable :: xm(:)
      real(dp), allocatable :: xpp(:)
      real(dp), allocatable :: xpm(:)
      real(dp), allocatable :: xmp(:)
      real(dp), allocatable :: xmm(:)
      real(dp) :: f0
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: fmm
      real(dp) :: fmp
      real(dp) :: fpm
      real(dp) :: fpp
      real(dp) :: hi
      real(dp) :: hj
      real(dp) :: scale
      integer :: i
      integer :: j
      integer :: n

      info = 0
      n = size(x)
      allocate(hessian(n, n))
      hessian = 0.0_dp
      if (n == 0) return
      if (present(lower)) then
         if (size(lower) /= n) then
            info = 1
            return
         end if
      end if
      if (present(upper)) then
         if (size(upper) /= n) then
            info = 1
            return
         end if
      end if
      scale = epsilon(1.0_dp)**0.25_dp
      if (present(relative_step)) scale = max(relative_step, 100.0_dp * epsilon(1.0_dp))
      allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n))
      f0 = problem%value(x)
      if (.not. ieee_is_finite(f0)) then
         info = 2
         return
      end if

      do i = 1, n
         hi = scale * max(1.0_dp, abs(x(i)))
         if (present(lower)) hi = min(hi, 0.45_dp * max(0.0_dp, x(i) - lower(i)))
         if (present(upper)) hi = min(hi, 0.45_dp * max(0.0_dp, upper(i) - x(i)))
         if (hi <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)))) then
            hessian(i, i) = huge(1.0_dp)
            cycle
         end if
         xp = x
         xm = x
         xp(i) = x(i) + hi
         xm(i) = x(i) - hi
         fp = problem%value(xp)
         fm = problem%value(xm)
         if (.not. ieee_is_finite(fp) .or. .not. ieee_is_finite(fm)) then
            info = 3
            return
         end if
         hessian(i, i) = (fp - 2.0_dp * f0 + fm) / hi**2

         do j = 1, i - 1
            hj = scale * max(1.0_dp, abs(x(j)))
            if (present(lower)) hj = min(hj, 0.45_dp * max(0.0_dp, x(j) - lower(j)))
            if (present(upper)) hj = min(hj, 0.45_dp * max(0.0_dp, upper(j) - x(j)))
            if (hj <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(j)))) cycle
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = x(i) + hi
            xpp(j) = x(j) + hj
            xpm(i) = x(i) + hi
            xpm(j) = x(j) - hj
            xmp(i) = x(i) - hi
            xmp(j) = x(j) + hj
            xmm(i) = x(i) - hi
            xmm(j) = x(j) - hj
            fpp = problem%value(xpp)
            fpm = problem%value(xpm)
            fmp = problem%value(xmp)
            fmm = problem%value(xmm)
            if (.not. ieee_is_finite(fpp) .or. .not. ieee_is_finite(fpm) .or. &
               .not. ieee_is_finite(fmp) .or. .not. ieee_is_finite(fmm)) then
               info = 4
               return
            end if
            hessian(i, j) = (fpp - fpm - fmp + fmm) / (4.0_dp * hi * hj)
            hessian(j, i) = hessian(i, j)
         end do
      end do
   end subroutine finite_difference_hessian

   pure subroutine projected_gradient(x, lower, upper, gradient, tolerance, projected)
      !! Applies first-order optimality conditions for a box-constrained minimization problem.
      real(dp), intent(in) :: x(:) !! Current bounded parameter vector.
      real(dp), intent(in) :: lower(:) !! Lower bounds corresponding to `x`.
      real(dp), intent(in) :: upper(:) !! Upper bounds corresponding to `x`.
      real(dp), intent(in) :: gradient(:) !! Unconstrained objective gradient.
      real(dp), intent(in) :: tolerance !! Boundary tolerance used to identify active constraints.
      real(dp), intent(out) :: projected(:) !! Gradient with components blocked by active constraints set to zero.
      integer :: i

      projected = gradient
      do i = 1, size(x)
         if (x(i) <= lower(i) + tolerance * (1.0_dp + abs(lower(i))) .and. gradient(i) > 0.0_dp) then
            projected(i) = 0.0_dp
         else if (x(i) >= upper(i) - tolerance * (1.0_dp + abs(upper(i))) .and. gradient(i) < 0.0_dp) then
            projected(i) = 0.0_dp
         end if
      end do
   end subroutine projected_gradient

   pure subroutine restrict_direction(x, lower, upper, direction, tolerance)
      !! Removes direction components that point outward through active box constraints.
      real(dp), intent(in) :: x(:) !! Current bounded parameter vector.
      real(dp), intent(in) :: lower(:) !! Lower bounds corresponding to `x`.
      real(dp), intent(in) :: upper(:) !! Upper bounds corresponding to `x`.
      real(dp), intent(inout) :: direction(:) !! Search direction modified in place to respect active constraints.
      real(dp), intent(in) :: tolerance !! Boundary tolerance used to identify active constraints.
      integer :: i

      do i = 1, size(x)
         if (x(i) <= lower(i) + tolerance * (1.0_dp + abs(lower(i))) .and. direction(i) < 0.0_dp) direction(i) = 0.0_dp
         if (x(i) >= upper(i) - tolerance * (1.0_dp + abs(upper(i))) .and. direction(i) > 0.0_dp) direction(i) = 0.0_dp
      end do
   end subroutine restrict_direction

   pure function outer_product(a, b) result(matrix)
      !! Forms the rank-one outer product `a*b^T`.
      real(dp), intent(in) :: a(:) !! Left vector of length `m`.
      real(dp), intent(in) :: b(:) !! Right vector of length `n`.
      real(dp) :: matrix(size(a), size(b))
      integer :: j

      do j = 1, size(b)
         matrix(:, j) = a * b(j)
      end do
   end function outer_product

end module ape_optimize
