! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim_mads
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf, ieee_negative_inf
   use dfoptim_kinds, only : dp
   use dfoptim_interfaces, only : mads_control_t, dfoptim_result_t, dfoptim_objective, &
      dfoptim_monitor, mads_poll_lite, mads_poll_full, dfoptim_success, &
      dfoptim_max_evaluations, dfoptim_invalid_input, dfoptim_nonfinite_objective, &
      dfoptim_cancelled
   use dfoptim_rng, only : rng_t
   use dfoptim_utils, only : evaluate_objective, call_monitor, status_message
   implicit none
   private
   public :: mads

contains

   function mads(par, fn, lower, upper, scale, control_in, user_data, monitor) result(result)
      real(dp), intent(in) :: par(:)
      procedure(dfoptim_objective) :: fn
      real(dp), intent(in), optional :: lower(:), upper(:), scale(:)
      type(mads_control_t), intent(in), optional :: control_in
      class(*), intent(inout), optional :: user_data
      procedure(dfoptim_monitor), optional :: monitor
      type(dfoptim_result_t) :: result
      type(mads_control_t) :: control
      type(rng_t) :: rng
      real(dp), allocatable :: lo(:), up(:), user_scale(:), solver_scale(:), offset(:), span(:)
      real(dp), allocatable :: best_x(:), current_x(:), poll_x(:, :), poll_y(:), direction(:)
      real(dp), allocatable :: line_x(:), temp_x(:)
      integer, allocatable :: fibonacci(:), log_eval(:), log_search(:)
      real(dp), allocatable :: log_delta(:), log_value(:)
      real(dp) :: best_f, candidate_f, delta, zoom
      integer :: nvar, poll_size, evaluations, iteration, i, k, best_index, search_success
      integer :: log_count, valid_poll
      logical :: bounded, found_better, finite_ok, cancelled, all_lo_finite, all_up_finite

      control = mads_control_t()
      if (present(control_in)) control = control_in
      nvar = size(par)
      allocate(result%x(nvar))
      result%x = par
      allocate(lo(nvar), up(nvar), user_scale(nvar), solver_scale(nvar), offset(nvar), span(nvar))

      lo = ieee_value(1.0_dp, ieee_negative_inf)
      up = ieee_value(1.0_dp, ieee_positive_inf)
      user_scale = 1.0_dp
      if (present(lower)) then
         if (size(lower) /= nvar) then
            call invalid_result(); return
         end if
         lo = lower
      end if
      if (present(upper)) then
         if (size(upper) /= nvar) then
            call invalid_result(); return
         end if
         up = upper
      end if
      if (present(scale)) then
         if (size(scale) /= nvar) then
            call invalid_result(); return
         end if
         user_scale = scale
      end if

      all_lo_finite = all(ieee_is_finite(lo))
      all_up_finite = all(ieee_is_finite(up))
      bounded = all_lo_finite .and. all_up_finite
      if ((all_lo_finite .neqv. all_up_finite) .or. &
          ((any(ieee_is_finite(lo)) .or. any(ieee_is_finite(up))) .and. .not. bounded)) then
         call invalid_result(); return
      end if

      if (nvar < 2 .or. any(lo > up) .or. any(par < lo) .or. any(par > up) .or. &
          any(user_scale <= control%tol) .or. control%tol <= 0.0_dp .or. &
          control%maxfeval < 1 .or. control%delta_init <= control%tol .or. &
          control%delta_init > 1.0_dp .or. control%expand <= 1.0_dp .or. &
          (control%poll_style /= mads_poll_lite .and. control%poll_style /= mads_poll_full)) then
         call invalid_result(); return
      end if
      if (bounded) then
         if (any(up <= lo) .or. any(user_scale >= up - lo)) then
            call invalid_result(); return
         end if
         offset = lo
         span = up - lo
         solver_scale = 2.0_dp * user_scale / span
      else
         offset = 0.0_dp
         span = 1.0_dp
         solver_scale = user_scale
      end if

      allocate(best_x(nvar), current_x(nvar), direction(nvar), line_x(nvar), temp_x(nvar))
      if (bounded) then
         best_x = 2.0_dp * (par - offset) / span - 1.0_dp
      else
         best_x = par
      end if

      call rng%seed(control%seed)
      delta = control%delta_init
      zoom = control%expand
      poll_size = merge(nvar + 1, 2 * nvar, control%poll_style == mads_poll_lite)
      allocate(poll_x(poll_size, nvar), poll_y(poll_size))

      if (control%line_search > 2) then
         allocate(fibonacci(control%line_search))
         fibonacci(1:2) = [1, 2]
         do i = 3, size(fibonacci)
            if (fibonacci(i - 1) > huge(0) - fibonacci(i - 2)) then
               fibonacci(i) = huge(0)
            else
               fibonacci(i) = fibonacci(i - 1) + fibonacci(i - 2)
            end if
         end do
      else
         allocate(fibonacci(0))
      end if

      allocate(log_eval(control%maxfeval + 2), log_search(control%maxfeval + 2))
      allocate(log_delta(control%maxfeval + 2), log_value(control%maxfeval + 2))

      call evaluate_point(best_x, best_f, finite_ok)
      evaluations = 1
      if (.not. finite_ok) then
         result%convergence = dfoptim_nonfinite_objective
         result%feval = 1
         result%message = status_message(result%convergence)
         return
      end if
      iteration = 0
      log_count = 1
      log_eval(1) = 1
      log_delta(1) = delta
      log_search(1) = 0
      log_value(1) = merge(-best_f, best_f, control%maximize)
      cancelled = .false.

      do while (delta > control%tol .and. evaluations < control%maxfeval)
         iteration = iteration + 1
         found_better = .false.
         search_success = 0
         current_x = best_x

         call build_poll_set(best_x, poll_size, delta, solver_scale, rng, poll_x)
         poll_y = huge(1.0_dp)
         valid_poll = 0
         do i = 1, poll_size
            if (evaluations >= control%maxfeval) exit
            if (bounded .and. any(poll_x(i, :) < -1.0_dp .or. poll_x(i, :) > 1.0_dp)) cycle
            call evaluate_point(poll_x(i, :), poll_y(i), finite_ok)
            evaluations = evaluations + 1
            valid_poll = valid_poll + 1
         end do

         if (valid_poll > 0) then
            best_index = minloc(poll_y, dim=1)
            if (poll_y(best_index) < best_f) then
               found_better = .true.
               best_f = poll_y(best_index)
               best_x = poll_x(best_index, :)
            end if
         end if

         if (found_better .and. size(fibonacci) > 0) then
            direction = best_x - current_x
            do k = 1, size(fibonacci)
               if (evaluations >= control%maxfeval) exit
               line_x = best_x + delta * solver_scale * real(fibonacci(k), dp) * direction
               if (bounded .and. any(line_x < -1.0_dp .or. line_x > 1.0_dp)) exit
               call evaluate_point(line_x, candidate_f, finite_ok)
               evaluations = evaluations + 1
               if (.not. finite_ok) exit
               if (candidate_f < best_f) then
                  best_f = candidate_f
                  best_x = line_x
                  search_success = fibonacci(k)
               else
                  exit
               end if
            end do
         end if

         if (found_better) then
            delta = min(1.0_dp, delta * zoom)
         else
            delta = delta / sqrt(zoom)
         end if

         log_count = log_count + 1
         log_eval(log_count) = evaluations
         log_delta(log_count) = delta
         log_search(log_count) = search_success
         log_value(log_count) = merge(-best_f, best_f, control%maximize)

         if (bounded) then
            temp_x = offset + 0.5_dp * (best_x + 1.0_dp) * span
         else
            temp_x = best_x
         end if
         if (control%trace) then
            write(*, '(a,i0,2x,a,i0,2x,a,es15.7,2x,a,es12.4)') &
               'iter=', iteration, 'feval=', evaluations, 'f=', &
               merge(-best_f, best_f, control%maximize), 'delta=', delta
         end if
         call call_monitor(monitor, temp_x, merge(-best_f, best_f, control%maximize), &
            iteration, evaluations, cancelled, user_data)
         if (cancelled) exit
      end do

      if (bounded) then
         result%x = offset + 0.5_dp * (best_x + 1.0_dp) * span
      else
         result%x = best_x
      end if
      result%value = merge(-best_f, best_f, control%maximize)
      result%feval = evaluations
      result%niter = iteration
      result%final_mesh = delta
      if (cancelled) then
         result%convergence = dfoptim_cancelled
      else if (delta <= control%tol) then
         result%convergence = dfoptim_success
      else
         result%convergence = dfoptim_max_evaluations
      end if
      result%message = status_message(result%convergence)
      allocate(result%log%evaluations(log_count), result%log%delta(log_count))
      allocate(result%log%search_success(log_count), result%log%value(log_count))
      result%log%evaluations = log_eval(1:log_count)
      result%log%delta = log_delta(1:log_count)
      result%log%search_success = log_search(1:log_count)
      result%log%value = log_value(1:log_count)

   contains

      subroutine invalid_result()
         result%convergence = dfoptim_invalid_input
         result%message = status_message(result%convergence)
      end subroutine invalid_result

      subroutine evaluate_point(z, value, ok)
         real(dp), intent(in) :: z(:)
         real(dp), intent(out) :: value
         logical, intent(out) :: ok
         real(dp) :: original(size(z))
         if (bounded) then
            original = offset + 0.5_dp * (z + 1.0_dp) * span
         else
            original = z
         end if
         value = evaluate_objective(fn, original, control%maximize, user_data, ok)
      end subroutine evaluate_point

   end function mads

   subroutine build_poll_set(center, npoints, mesh_size, custom_scale, rng, points)
      real(dp), intent(in) :: center(:), mesh_size, custom_scale(:)
      integer, intent(in) :: npoints
      type(rng_t), intent(inout) :: rng
      real(dp), intent(out) :: points(:, :)
      real(dp), allocatable :: temp(:, :), directions(:, :), average(:)
      integer, allocatable :: p(:)
      integer :: m, i, j, step_size
      real(dp) :: u

      m = size(center)
      allocate(temp(m, m), average(m), p(m))
      temp = 0.0_dp
      if (mesh_size > 1.0_dp) then
         step_size = 1
      else
         step_size = max(1, floor(real(m, dp) / mesh_size))
      end if

      do i = 1, m
         u = rng%uniform()
         temp(i, i) = real(step_size, dp) * merge(1.0_dp, -1.0_dp, u < 0.5_dp)
         do j = 1, i - 1
            temp(i, j) = floor(real(step_size, dp) * (1.0_dp - 2.0_dp * rng%uniform()))
         end do
      end do

      ! The R implementation applies two sample() shuffles.  Row and column
      ! permutations preserve the lower-triangular integer basis distribution.
      call rng%permutation(p)
      temp = temp(p, :)
      call rng%permutation(p)
      temp = temp(:, p)

      allocate(directions(npoints, m))
      if (npoints == m + 1) then
         directions(1:m, :) = mesh_size * temp
         average = sum(directions(1:m, :), dim=1) / real(m, dp)
         do j = 1, m
            if (average(j) > 0.0_dp) then
               directions(m + 1, j) = -1.0_dp
            else if (average(j) < 0.0_dp) then
               directions(m + 1, j) = 1.0_dp
            else
               directions(m + 1, j) = 0.0_dp
            end if
         end do
      else
         directions(1:m, :) = mesh_size * temp
         directions(m + 1:2 * m, :) = -mesh_size * temp
      end if

      do i = 1, npoints
         points(i, :) = center + mesh_size * custom_scale * directions(i, :)
      end do
   end subroutine build_poll_set

end module dfoptim_mads
