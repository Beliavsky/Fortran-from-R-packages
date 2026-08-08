! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim_hooke_jeeves
   use dfoptim_kinds, only : dp
   use dfoptim_interfaces, only : hj_control_t, dfoptim_result_t, dfoptim_objective, &
      dfoptim_monitor, dfoptim_success, dfoptim_max_evaluations, dfoptim_invalid_input, &
      dfoptim_nonfinite_objective, dfoptim_cancelled
   use dfoptim_rng, only : rng_t
   use dfoptim_utils, only : evaluate_objective, call_monitor, status_message
   implicit none
   private
   public :: hjk, hjkb

contains

   function hjk(par, fn, control, user_data, monitor) result(result)
      real(dp), intent(in) :: par(:)
      procedure(dfoptim_objective) :: fn
      type(hj_control_t), intent(in), optional :: control
      class(*), intent(inout), optional :: user_data
      procedure(dfoptim_monitor), optional :: monitor
      type(dfoptim_result_t) :: result
      real(dp), allocatable :: lower(:), upper(:)

      allocate(lower(size(par)), upper(size(par)))
      lower = -huge(1.0_dp)
      upper = huge(1.0_dp)
      result = hooke_jeeves_impl(par, fn, lower, upper, .false., control, user_data, monitor)
   end function hjk

   function hjkb(par, fn, lower, upper, control, user_data, monitor) result(result)
      real(dp), intent(in) :: par(:)
      procedure(dfoptim_objective) :: fn
      real(dp), intent(in) :: lower(:), upper(:)
      type(hj_control_t), intent(in), optional :: control
      class(*), intent(inout), optional :: user_data
      procedure(dfoptim_monitor), optional :: monitor
      type(dfoptim_result_t) :: result

      result = hooke_jeeves_impl(par, fn, lower, upper, .true., control, user_data, monitor)
   end function hjkb

   function hooke_jeeves_impl(par, fn, lower, upper, bounded, control_in, user_data, monitor) result(result)
      real(dp), intent(in) :: par(:), lower(:), upper(:)
      procedure(dfoptim_objective) :: fn
      logical, intent(in) :: bounded
      type(hj_control_t), intent(in), optional :: control_in
      class(*), intent(inout), optional :: user_data
      procedure(dfoptim_monitor), optional :: monitor
      type(dfoptim_result_t) :: result
      type(hj_control_t) :: control
      type(rng_t) :: rng
      real(dp), allocatable :: x(:), direction(:, :)
      real(dp) :: fx, h
      integer :: n, nsteps, ns, calls
      logical :: finite_ok, cancelled, improved

      control = hj_control_t()
      if (present(control_in)) control = control_in
      n = size(par)
      allocate(result%x(n))
      result%x = par

      if (n < 2 .or. size(lower) /= n .or. size(upper) /= n .or. &
          control%tol <= 0.0_dp .or. control%tol >= 1.0_dp .or. &
          control%maxfeval < 1 .or. any(lower > upper) .or. &
          (bounded .and. (any(par < lower) .or. any(par > upper)))) then
         result%convergence = dfoptim_invalid_input
         result%message = status_message(result%convergence)
         return
      end if

      allocate(x(n), direction(n, n))
      x = par
      direction = 0.0_dp
      do ns = 1, n
         direction(ns, ns) = 1.0_dp
      end do
      call rng%seed(control%seed)

      fx = evaluate_objective(fn, x, control%maximize, user_data, finite_ok)
      result%feval = 1
      if (.not. finite_ok) then
         result%convergence = dfoptim_nonfinite_objective
         result%message = status_message(result%convergence)
         return
      end if

      nsteps = floor(log(1.0_dp / control%tol) / log(2.0_dp))
      nsteps = max(1, nsteps)
      ns = 0
      cancelled = .false.

      if (control%trace) write(*, '(a)') 'step  feval                 value  first_parameter'
      do while (ns < nsteps .and. result%feval < control%maxfeval .and. abs(fx) < control%target)
         ns = ns + 1
         h = 2.0_dp ** (-(ns - 1))
         call hj_search(x, fx, h, direction, lower, upper, bounded, fn, control, rng, &
            result%feval, calls, improved, user_data)
         result%feval = result%feval + calls

         if (control%trace) then
            write(*, '(i5,1x,i7,1x,es22.13,1x,es18.9)') ns, result%feval, &
               merge(-fx, fx, control%maximize), x(1)
         end if
         call call_monitor(monitor, x, merge(-fx, fx, control%maximize), ns, &
            result%feval, cancelled, user_data)
         if (cancelled) exit
      end do

      result%x = x
      result%value = merge(-fx, fx, control%maximize)
      result%niter = ns
      if (cancelled) then
         result%convergence = dfoptim_cancelled
      else if (result%feval >= control%maxfeval .and. ns < nsteps) then
         result%convergence = dfoptim_max_evaluations
      else if (abs(fx) >= control%target) then
         result%convergence = dfoptim_max_evaluations
      else
         result%convergence = dfoptim_success
      end if
      result%message = status_message(result%convergence)
   end function hooke_jeeves_impl

   subroutine hj_search(x, fx, h, direction, lower, upper, bounded, fn, control, rng, &
      feval_before, calls, improved, user_data)
      real(dp), intent(inout) :: x(:)
      real(dp), intent(inout) :: fx
      real(dp), intent(in) :: h, direction(:, :), lower(:), upper(:)
      logical, intent(in) :: bounded
      procedure(dfoptim_objective) :: fn
      type(hj_control_t), intent(in) :: control
      type(rng_t), intent(inout) :: rng
      integer, intent(in) :: feval_before
      integer, intent(out) :: calls
      logical, intent(out) :: improved
      class(*), intent(inout), optional :: user_data
      real(dp), allocatable :: xb(:), xc(:), trial_x(:)
      real(dp) :: trial_f, fb
      integer :: local_calls
      logical :: local_improved

      allocate(xb(size(x)), xc(size(x)), trial_x(size(x)))
      xb = x
      xc = x
      calls = 0
      call hj_explore(xb, xc, fx, h, direction, lower, upper, bounded, fn, control, rng, &
         feval_before + calls, local_calls, trial_x, trial_f, local_improved, user_data)
      calls = calls + local_calls
      x = trial_x
      fx = trial_f
      improved = local_improved

      do while (improved)
         xc = x + (x - xb)
         if (bounded) xc = max(lower, min(upper, xc))
         xb = x
         fb = fx
         call hj_explore(xb, xc, fb, h, direction, lower, upper, bounded, fn, control, rng, &
            feval_before + calls, local_calls, trial_x, trial_f, local_improved, user_data)
         calls = calls + local_calls
         x = trial_x
         fx = trial_f
         improved = local_improved

         if (.not. improved .and. feval_before + calls < control%maxfeval) then
            call hj_explore(xb, xb, fb, h, direction, lower, upper, bounded, fn, control, rng, &
               feval_before + calls, local_calls, trial_x, trial_f, local_improved, user_data)
            calls = calls + local_calls
            x = trial_x
            fx = trial_f
            improved = local_improved
         end if
         if (feval_before + calls >= control%maxfeval .or. abs(fx) >= control%target) exit
      end do
   end subroutine hj_search

   subroutine hj_explore(xb, xc, fbold, h, direction, lower, upper, bounded, fn, control, rng, &
      feval_before, calls, x, fx, improved, user_data)
      real(dp), intent(in) :: xb(:), xc(:), fbold, h, direction(:, :), lower(:), upper(:)
      logical, intent(in) :: bounded
      procedure(dfoptim_objective) :: fn
      type(hj_control_t), intent(in) :: control
      type(rng_t), intent(inout) :: rng
      integer, intent(in) :: feval_before
      integer, intent(out) :: calls
      real(dp), intent(out) :: x(:), fx
      logical, intent(out) :: improved
      class(*), intent(inout), optional :: user_data
      integer, allocatable :: permutation(:)
      real(dp), allocatable :: xt(:), p1(:), p2(:)
      real(dp) :: best, f1, f2
      integer :: j, k
      logical :: finite_ok

      allocate(permutation(size(xb)), xt(size(xb)), p1(size(xb)), p2(size(xb)))
      call rng%permutation(permutation)
      xt = xc
      best = fbold
      calls = 0
      improved = .false.

      do j = 1, size(permutation)
         k = permutation(j)
         p1 = xt + h * direction(:, k)
         p2 = xt - h * direction(:, k)
         f1 = best
         f2 = best

         if ((.not. bounded .or. all(p1 >= lower .and. p1 <= upper)) .and. &
             feval_before + calls < control%maxfeval) then
            f1 = evaluate_objective(fn, p1, control%maximize, user_data, finite_ok)
            calls = calls + 1
         end if
         if ((.not. bounded .or. all(p2 >= lower .and. p2 <= upper)) .and. &
             feval_before + calls < control%maxfeval) then
            f2 = evaluate_objective(fn, p2, control%maximize, user_data, finite_ok)
            calls = calls + 1
         end if

         if (min(f1, f2) < best) then
            improved = .true.
            if (f1 < f2) then
               xt = p1
               best = f1
            else
               xt = p2
               best = f2
            end if
         end if
         if (feval_before + calls >= control%maxfeval) exit
      end do

      if (improved) then
         x = xt
         fx = best
      else
         x = xb
         fx = fbold
      end if
   end subroutine hj_explore

end module dfoptim_hooke_jeeves
