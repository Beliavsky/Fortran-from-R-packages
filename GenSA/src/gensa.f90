! SPDX-License-Identifier: GPL-2.0-only
module gensa
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gensa_kinds, only : dp, i8
   use gensa_rng, only : ran2_state, gensa_visit
   use gensa_types
   use gensa_local, only : projected_bfgs, bounded_pattern_search
   implicit none
   private

   public :: dp, i8
   public :: gensa_control, gensa_trace, gensa_result
   public :: gensa_objective, gensa_constraint
   public :: gensa_minimize
   public :: gensa_success, gensa_max_iterations, gensa_max_calls
   public :: gensa_max_time, gensa_threshold_reached, gensa_no_improvement
   public :: gensa_invalid_input, gensa_no_feasible_start

contains

   subroutine gensa_minimize(objective, lower, upper, result, control, initial, constraint)
      procedure(gensa_objective) :: objective
      real(dp), intent(in) :: lower(:), upper(:)
      type(gensa_result), intent(out) :: result
      type(gensa_control), intent(in), optional :: control
      real(dp), intent(in), optional :: initial(:)
      procedure(gensa_constraint), optional :: constraint

      type(gensa_control) :: settings
      type(ran2_state) :: rng
      real(dp), allocatable :: current(:), candidate(:), backup(:), best(:), range(:)
      real(dp), allocatable :: markov_best(:)
      real(dp) :: current_value, candidate_value, best_value, markov_best_value
      real(dp) :: temperature, temp_accept, pbase, probability, u, visit
      real(dp) :: elapsed, clock_start, local_value
      integer :: n, markov_length, local_maxit, iteration, j, k, age
      integer :: attempts, no_improvement, status_code, last_completed
      logical :: feasible, stop_now, improved_outer, halted, accepted
      logical :: initialized, local_improved

      settings = gensa_control()
      if (present(control)) settings = control
      call initialize_empty_result(result)

      n = size(lower)
      if (n <= 0 .or. size(upper) /= n) then
         call fail_result(result, gensa_invalid_input, 'lower and upper must have equal positive size')
         return
      end if
      if (any(.not. ieee_is_finite(lower)) .or. any(.not. ieee_is_finite(upper)) .or. any(lower >= upper)) then
         call fail_result(result, gensa_invalid_input, 'bounds must be finite and satisfy lower < upper')
         return
      end if
      if (present(initial)) then
         if (size(initial) /= n .or. any(.not. ieee_is_finite(initial))) then
            call fail_result(result, gensa_invalid_input, 'initial point has invalid size or values')
            return
         end if
      end if
      if (settings%maxit < 1 .or. settings%max_call < 1 .or. settings%temperature <= 0.0_dp) then
         call fail_result(result, gensa_invalid_input, 'iteration, call, and temperature controls must be positive')
         return
      end if
      if (settings%visiting_param <= 1.0_dp .or. settings%visiting_param >= 3.0_dp) then
         call fail_result(result, gensa_invalid_input, 'visiting_param must lie strictly between 1 and 3')
         return
      end if
      if (settings%acceptance_param >= 1.0_dp) then
         call fail_result(result, gensa_invalid_input, 'acceptance_param must be less than 1')
         return
      end if

      allocate(current(n), candidate(n), backup(n), best(n), range(n), markov_best(n))
      range = upper - lower
      markov_length = settings%markov_length
      if (markov_length <= 0) markov_length = 2 * n
      if (mod(markov_length, n) /= 0) then
         call fail_result(result, gensa_invalid_input, 'markov_length must be a multiple of the dimension')
         return
      end if
      local_maxit = settings%local_maxit
      if (local_maxit <= 0) local_maxit = min(1000, max(100, 6 * n))

      call rng%seed(settings%seed)
      call cpu_time(clock_start)
      status_code = gensa_max_iterations
      last_completed = 0

      initialized = .false.
      do attempts = 1, settings%max_constraint_attempts
         if (present(initial) .and. attempts == 1) then
            current = min(upper, max(lower, initial))
         else
            do k = 1, n
               current(k) = lower(k) + rng%uniform() * range(k)
            end do
         end if
         feasible = is_feasible(current)
         if (.not. feasible) cycle
         call evaluate(current, current_value, stop_now)
         if (stop_now .and. result%counts == 0) cycle
         if (current_value < huge(1.0_dp) / 10.0_dp) then
            initialized = .true.
            exit
         end if
      end do
      if (.not. initialized) then
         call fail_result(result, gensa_no_feasible_start, 'could not find a finite feasible starting point')
         return
      end if

      best = current
      best_value = current_value

      if (settings%local_search) then
         candidate = best
         call run_local(candidate, local_value, halted)
         if (local_value < best_value) then
            best = candidate
            best_value = local_value
            current = candidate
            current_value = local_value
         end if
         if (halted) then
            call determine_stop(status_code)
            call finalize_result(last_completed)
            return
         end if
      end if

      if (settings%trace) then
         allocate(result%trace%step(settings%maxit + 1))
         allocate(result%trace%temperature(settings%maxit + 1))
         allocate(result%trace%current_value(settings%maxit + 1))
         allocate(result%trace%best_value(settings%maxit + 1))
         result%trace%n = 1
         result%trace%step(1) = 0
         result%trace%temperature(1) = settings%temperature
         result%trace%current_value(1) = current_value
         result%trace%best_value(1) = best_value
      end if

      if (threshold_met()) then
         status_code = gensa_threshold_reached
         call finalize_result(0)
         return
      end if

      age = 1
      no_improvement = 0
      do iteration = 1, settings%maxit
         last_completed = iteration
         temperature = annealing_temperature(age)
         if (temperature < settings%temp_restart) then
            age = 1
            temperature = annealing_temperature(age)
         end if
         temp_accept = temperature / real(max(1, age), dp)
         improved_outer = .false.
         markov_best = current
         markov_best_value = current_value

         do j = 1, markov_length
            backup = current
            candidate = current
            accepted = .false.
            do attempts = 1, settings%max_constraint_attempts
               candidate = backup
               if (j <= n) then
                  do k = 1, n
                     visit = clipped_visit(temperature)
                     candidate(k) = wrap_bound(backup(k) + visit, lower(k), range(k))
                  end do
               else
                  k = modulo(j - n - 1, n) + 1
                  visit = clipped_visit(temperature)
                  candidate(k) = wrap_bound(backup(k) + visit, lower(k), range(k))
               end if
               if (is_feasible(candidate)) exit
            end do
            if (attempts > settings%max_constraint_attempts) cycle

            call evaluate(candidate, candidate_value, stop_now)
            if (stop_now) then
               call determine_stop(status_code)
               call finalize_result(iteration - 1)
               return
            end if

            if (candidate_value < current_value) then
               accepted = .true.
            else
               pbase = 1.0_dp + (settings%acceptance_param - 1.0_dp) &
                  * (candidate_value - current_value) / max(temp_accept, tiny(1.0_dp))
               if (pbase > 0.0_dp) then
                  probability = exp(log(pbase) / (1.0_dp - settings%acceptance_param))
               else
                  probability = 0.0_dp
               end if
               u = rng%uniform()
               accepted = u <= min(1.0_dp, probability)
            end if

            if (accepted) then
               current = candidate
               current_value = candidate_value
            end if
            if (current_value < markov_best_value) then
               markov_best = current
               markov_best_value = current_value
            end if
            if (candidate_value < best_value) then
               best = candidate
               best_value = candidate_value
               improved_outer = .true.
            end if

            if (threshold_met()) then
               status_code = gensa_threshold_reached
               call finalize_result(iteration)
               return
            end if
         end do

         local_improved = .false.
         if (settings%local_search .and. improved_outer) then
            candidate = best
            call run_local(candidate, local_value, halted)
            if (local_value < best_value) then
               best = candidate
               best_value = local_value
               current = candidate
               current_value = local_value
               local_improved = .true.
            end if
            if (halted) then
               call determine_stop(status_code)
               call record_trace(iteration, temperature)
               call finalize_result(iteration)
               return
            end if
         else if (settings%local_search .and. no_improvement >= merge(n, 1000, settings%simple_function)) then
            candidate = markov_best
            call run_local(candidate, local_value, halted)
            no_improvement = 0
            if (local_value < best_value) then
               best = candidate
               best_value = local_value
               current = candidate
               current_value = local_value
               local_improved = .true.
            end if
            if (halted) then
               call determine_stop(status_code)
               call record_trace(iteration, temperature)
               call finalize_result(iteration)
               return
            end if
         end if

         if (improved_outer .or. local_improved) then
            no_improvement = 0
         else
            no_improvement = no_improvement + 1
         end if

         call record_trace(iteration, temperature)
         if (settings%verbose .and. mod(iteration, max(1, settings%report)) == 0) then
            write(*, '(a,i0,a,es14.6,a,es14.6,a,i0)') 'iteration ', iteration, &
               ', current=', current_value, ', best=', best_value, ', calls=', result%counts
         end if

         if (settings%no_improvement_stop >= 0 .and. no_improvement >= settings%no_improvement_stop) then
            status_code = gensa_no_improvement
            call finalize_result(iteration)
            return
         end if
         call cpu_time(elapsed)
         if (elapsed - clock_start >= settings%max_time) then
            status_code = gensa_max_time
            call finalize_result(iteration)
            return
         end if
         if (result%counts >= settings%max_call) then
            status_code = gensa_max_calls
            call finalize_result(iteration)
            return
         end if
         age = age + 1
      end do

      status_code = gensa_max_iterations
      call finalize_result(settings%maxit)

   contains

      subroutine evaluate(x, value, halt)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: value
         logical, intent(out) :: halt
         real(dp) :: now

         halt = .false.
         if (result%counts >= settings%max_call) then
            value = huge(1.0_dp)
            halt = .true.
            return
         end if
         call cpu_time(now)
         if (now - clock_start >= settings%max_time) then
            value = huge(1.0_dp)
            halt = .true.
            return
         end if
         if (.not. is_feasible(x)) then
            value = huge(1.0_dp)
            return
         end if
         value = objective(x)
         result%counts = result%counts + 1
         if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
      end subroutine evaluate

      logical function is_feasible(x)
         real(dp), intent(in) :: x(:)
         is_feasible = all(x >= lower) .and. all(x <= upper)
         if (is_feasible .and. present(constraint)) is_feasible = constraint(x)
      end function is_feasible

      subroutine run_local(x, value, halt)
         real(dp), intent(inout) :: x(:)
         real(dp), intent(out) :: value
         logical, intent(out) :: halt
         if (settings%smooth) then
            call projected_bfgs(evaluate, x, lower, upper, local_maxit, &
               settings%local_tolerance, value, halt)
         else
            call bounded_pattern_search(evaluate, x, lower, upper, local_maxit, &
               settings%local_tolerance, value, halt)
         end if
      end subroutine run_local

      real(dp) function annealing_temperature(cycle_age)
         integer, intent(in) :: cycle_age
         real(dp) :: numerator, denominator
         numerator = exp((settings%visiting_param - 1.0_dp) * log(2.0_dp)) - 1.0_dp
         denominator = exp((settings%visiting_param - 1.0_dp) &
            * log(real(cycle_age + 1, dp))) - 1.0_dp
         annealing_temperature = settings%temperature * numerator / denominator
      end function annealing_temperature

      real(dp) function clipped_visit(temp)
         real(dp), intent(in) :: temp
         clipped_visit = gensa_visit(settings%visiting_param, temp, rng)
         if (clipped_visit > 1.0e8_dp) clipped_visit = 1.0e8_dp * rng%uniform()
         if (clipped_visit < -1.0e8_dp) clipped_visit = -1.0e8_dp * rng%uniform()
      end function clipped_visit

      real(dp) function wrap_bound(x, low, width)
         real(dp), intent(in) :: x, low, width
         wrap_bound = low + modulo(x - low, width)
         if (abs(wrap_bound - low) < 1.0e-10_dp) wrap_bound = min(low + 1.0e-10_dp, low + width)
      end function wrap_bound

      logical function threshold_met()
         threshold_met = settings%has_threshold .and. best_value <= settings%threshold_stop
      end function threshold_met

      subroutine determine_stop(code)
         integer, intent(out) :: code
         real(dp) :: now
         call cpu_time(now)
         if (threshold_met()) then
            code = gensa_threshold_reached
         else if (result%counts >= settings%max_call) then
            code = gensa_max_calls
         else if (now - clock_start >= settings%max_time) then
            code = gensa_max_time
         else
            code = gensa_success
         end if
      end subroutine determine_stop

      subroutine record_trace(step, temp)
         integer, intent(in) :: step
         real(dp), intent(in) :: temp
         integer :: index
         if (.not. settings%trace) return
         index = result%trace%n + 1
         if (index > size(result%trace%step)) return
         result%trace%n = index
         result%trace%step(index) = step
         result%trace%temperature(index) = temp
         result%trace%current_value(index) = current_value
         result%trace%best_value(index) = best_value
      end subroutine record_trace

      subroutine finalize_result(completed)
         integer, intent(in) :: completed
         result%par = best
         result%value = best_value
         result%iterations = completed
         result%status = status_code
         result%message = status_message(status_code)
         if (settings%trace .and. result%trace%n > 0) call trim_trace(result%trace)
      end subroutine finalize_result

   end subroutine gensa_minimize

   subroutine initialize_empty_result(result)
      type(gensa_result), intent(out) :: result
      result%value = huge(1.0_dp)
      result%counts = 0
      result%iterations = 0
      result%status = gensa_invalid_input
      result%message = 'not initialized'
      result%trace%n = 0
   end subroutine initialize_empty_result

   subroutine fail_result(result, code, message)
      type(gensa_result), intent(inout) :: result
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      result%status = code
      result%message = message
   end subroutine fail_result

   pure function status_message(code) result(message)
      integer, intent(in) :: code
      character(len=:), allocatable :: message
      select case (code)
      case (gensa_success)
         message = 'search completed'
      case (gensa_max_iterations)
         message = 'maximum iterations reached'
      case (gensa_max_calls)
         message = 'maximum objective-function calls reached'
      case (gensa_max_time)
         message = 'maximum CPU time reached'
      case (gensa_threshold_reached)
         message = 'objective threshold reached'
      case (gensa_no_improvement)
         message = 'no-improvement limit reached'
      case (gensa_no_feasible_start)
         message = 'no feasible finite starting point found'
      case default
         message = 'invalid input'
      end select
   end function status_message

   subroutine trim_trace(trace)
      type(gensa_trace), intent(inout) :: trace
      integer, allocatable :: steps(:)
      real(dp), allocatable :: temperatures(:), currents(:), bests(:)
      integer :: n
      n = trace%n
      allocate(steps(n), temperatures(n), currents(n), bests(n))
      steps = trace%step(:n)
      temperatures = trace%temperature(:n)
      currents = trace%current_value(:n)
      bests = trace%best_value(:n)
      call move_alloc(steps, trace%step)
      call move_alloc(temperatures, trace%temperature)
      call move_alloc(currents, trace%current_value)
      call move_alloc(bests, trace%best_value)
   end subroutine trim_trace

end module gensa
