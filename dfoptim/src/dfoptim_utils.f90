! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use dfoptim_kinds, only : dp
   use dfoptim_interfaces, only : dfoptim_objective, dfoptim_monitor
   implicit none
   private
   public :: evaluate_objective, call_monitor, sort_simplex, max_abs, status_message

contains

   function evaluate_objective(fn, x, maximize, user_data, finite_ok) result(value)
      procedure(dfoptim_objective) :: fn
      real(dp), intent(in) :: x(:)
      logical, intent(in) :: maximize
      class(*), intent(inout), optional :: user_data
      logical, intent(out) :: finite_ok
      real(dp) :: value

      if (present(user_data)) then
         value = fn(x, user_data)
      else
         value = fn(x)
      end if
      finite_ok = ieee_is_finite(value)
      if (.not. finite_ok) then
         value = huge(1.0_dp)
      else if (maximize) then
         value = -value
      end if
   end function evaluate_objective

   subroutine call_monitor(monitor, x, value, iteration, evaluations, stop, user_data)
      procedure(dfoptim_monitor), optional :: monitor
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: value
      integer, intent(in) :: iteration, evaluations
      logical, intent(out) :: stop
      class(*), intent(inout), optional :: user_data

      stop = .false.
      if (.not. present(monitor)) return
      if (present(user_data)) then
         call monitor(x, value, iteration, evaluations, stop, user_data)
      else
         call monitor(x, value, iteration, evaluations, stop)
      end if
   end subroutine call_monitor

   subroutine sort_simplex(v, f)
      real(dp), intent(inout) :: v(:, :)
      real(dp), intent(inout) :: f(:)
      integer :: i, j, best
      real(dp) :: tf
      real(dp), allocatable :: tx(:)

      allocate(tx(size(v, 1)))
      do i = 1, size(f) - 1
         best = i
         do j = i + 1, size(f)
            if (f(j) < f(best)) best = j
         end do
         if (best /= i) then
            tf = f(i); f(i) = f(best); f(best) = tf
            tx = v(:, i); v(:, i) = v(:, best); v(:, best) = tx
         end if
      end do
   end subroutine sort_simplex

   pure function max_abs(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = maxval(abs(x))
      end if
   end function max_abs

   pure function status_message(code) result(message)
      integer, intent(in) :: code
      character(len=:), allocatable :: message
      select case (code)
      case (0)
         message = 'Successful convergence'
      case (1)
         message = 'Maximum number of function evaluations exceeded'
      case (2)
         message = 'Stagnation or restart limit reached'
      case (-1)
         message = 'Invalid input'
      case (-2)
         message = 'Objective returned a non-finite value at the starting point'
      case (-3)
         message = 'Cancelled by monitor callback'
      case default
         message = 'Unknown status'
      end select
   end function status_message

end module dfoptim_utils
