! SPDX-License-Identifier: GPL-2.0-only
module gensa_local
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gensa_kinds, only : dp
   implicit none
   private

   abstract interface
      subroutine local_evaluate(x, value, halt)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: value
         logical, intent(out) :: halt
      end subroutine local_evaluate
   end interface

   public :: projected_bfgs, bounded_pattern_search

contains

   subroutine projected_bfgs(evaluate, x, lower, upper, maxit, tolerance, value, halted)
      procedure(local_evaluate) :: evaluate
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: lower(:), upper(:)
      integer, intent(in) :: maxit
      real(dp), intent(in) :: tolerance
      real(dp), intent(out) :: value
      logical, intent(out) :: halted
      integer :: n, iter, ls
      real(dp) :: alpha, directional, ys, rho, ftrial
      real(dp), allocatable :: gradient(:), gradient_new(:), projected(:)
      real(dp), allocatable :: direction(:), trial(:), s(:), y(:)
      real(dp), allocatable :: hessian_inv(:, :), identity(:, :), left(:, :), right(:, :)
      logical :: stop_now, accepted

      n = size(x)
      allocate(gradient(n), gradient_new(n), projected(n), direction(n), trial(n), s(n), y(n))
      allocate(hessian_inv(n, n), identity(n, n), left(n, n), right(n, n))
      identity = 0.0_dp
      do iter = 1, n
         identity(iter, iter) = 1.0_dp
      end do
      hessian_inv = identity
      halted = .false.

      call evaluate(x, value, stop_now)
      if (stop_now) then
         halted = .true.
         return
      end if
      call numerical_gradient(evaluate, x, lower, upper, gradient, stop_now)
      if (stop_now) then
         halted = .true.
         return
      end if

      do iter = 1, maxit
         projected = gradient
         where (x <= lower + 10.0_dp * epsilon(1.0_dp) .and. gradient > 0.0_dp) projected = 0.0_dp
         where (x >= upper - 10.0_dp * epsilon(1.0_dp) .and. gradient < 0.0_dp) projected = 0.0_dp
         if (maxval(abs(projected)) <= tolerance) exit

         direction = -matmul(hessian_inv, gradient)
         if (dot_product(direction, gradient) >= -epsilon(1.0_dp)) direction = -gradient

         alpha = 1.0_dp
         accepted = .false.
         do ls = 1, 30
            trial = min(upper, max(lower, x + alpha * direction))
            s = trial - x
            if (maxval(abs(s)) <= epsilon(1.0_dp)) then
               alpha = 0.5_dp * alpha
               cycle
            end if
            directional = dot_product(gradient, s)
            call evaluate(trial, ftrial, stop_now)
            if (stop_now) then
               halted = .true.
               return
            end if
            if (ftrial <= value + 1.0e-4_dp * directional) then
               accepted = .true.
               exit
            end if
            alpha = 0.5_dp * alpha
         end do
         if (.not. accepted) exit

         call numerical_gradient(evaluate, trial, lower, upper, gradient_new, stop_now)
         if (stop_now) then
            halted = .true.
            return
         end if
         s = trial - x
         y = gradient_new - gradient
         ys = dot_product(y, s)
         if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, sqrt(dot_product(s, s) * dot_product(y, y)))) then
            rho = 1.0_dp / ys
            left = identity - rho * outer_product(s, y)
            right = identity - rho * outer_product(y, s)
            hessian_inv = matmul(left, matmul(hessian_inv, right)) + rho * outer_product(s, s)
         else
            hessian_inv = identity
         end if

         x = trial
         value = ftrial
         gradient = gradient_new
      end do
   end subroutine projected_bfgs

   subroutine bounded_pattern_search(evaluate, x, lower, upper, maxit, tolerance, value, halted)
      procedure(local_evaluate) :: evaluate
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: lower(:), upper(:)
      integer, intent(in) :: maxit
      real(dp), intent(in) :: tolerance
      real(dp), intent(out) :: value
      logical, intent(out) :: halted
      real(dp), allocatable :: step(:), trial(:), best_x(:)
      real(dp) :: trial_value, best_value
      integer :: iter, i, sign_index
      logical :: improved, stop_now

      allocate(step(size(x)), trial(size(x)), best_x(size(x)))
      step = max(0.2_dp * (upper - lower), 10.0_dp * tolerance)
      call evaluate(x, value, stop_now)
      if (stop_now) then
         halted = .true.
         return
      end if
      halted = .false.

      do iter = 1, maxit
         improved = .false.
         best_value = value
         best_x = x
         do i = 1, size(x)
            do sign_index = -1, 1, 2
               trial = x
               trial(i) = min(upper(i), max(lower(i), x(i) + real(sign_index, dp) * step(i)))
               if (abs(trial(i) - x(i)) <= epsilon(1.0_dp)) cycle
               call evaluate(trial, trial_value, stop_now)
               if (stop_now) then
                  halted = .true.
                  return
               end if
               if (trial_value < best_value) then
                  best_value = trial_value
                  best_x = trial
                  improved = .true.
               end if
            end do
         end do
         if (improved) then
            x = best_x
            value = best_value
         else
            step = 0.5_dp * step
         end if
         if (maxval(step) <= tolerance * max(1.0_dp, maxval(upper - lower))) exit
      end do
   end subroutine bounded_pattern_search

   subroutine numerical_gradient(evaluate, x, lower, upper, gradient, halted)
      procedure(local_evaluate) :: evaluate
      real(dp), intent(in) :: x(:), lower(:), upper(:)
      real(dp), intent(out) :: gradient(:)
      logical, intent(out) :: halted
      real(dp), allocatable :: left(:), right(:)
      real(dp) :: fleft, fright, h, denominator
      integer :: i
      logical :: stop_now

      allocate(left(size(x)), right(size(x)))
      halted = .false.
      do i = 1, size(x)
         h = max(1.0e-6_dp * max(1.0_dp, abs(x(i))), sqrt(epsilon(1.0_dp)) * (upper(i) - lower(i)))
         left = x
         right = x
         left(i) = max(lower(i), x(i) - h)
         right(i) = min(upper(i), x(i) + h)
         denominator = right(i) - left(i)
         if (denominator <= epsilon(1.0_dp)) then
            gradient(i) = 0.0_dp
            cycle
         end if
         call evaluate(left, fleft, stop_now)
         if (stop_now) then
            halted = .true.
            return
         end if
         call evaluate(right, fright, stop_now)
         if (stop_now) then
            halted = .true.
            return
         end if
         gradient(i) = (fright - fleft) / denominator
         if (.not. ieee_is_finite(gradient(i))) gradient(i) = 0.0_dp
      end do
   end subroutine numerical_gradient

   pure function outer_product(a, b) result(matrix)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: matrix(size(a), size(b))
      matrix = spread(a, 2, size(b)) * spread(b, 1, size(a))
   end function outer_product

end module gensa_local
