! SPDX-License-Identifier: GPL-2.0-or-later
module jdmbs_model
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jdmbs_kinds, only : dp, int64
   use jdmbs_status, only : jdmbs_success, jdmbs_invalid_argument, &
      jdmbs_nonfinite_input, jdmbs_numerical_warning
   use jdmbs_rng, only : rng_state, rng_seed, rng_uniform, rng_normal, rng_integer
   implicit none
   private

   type, public :: jdmbs_control
      integer :: day = 180
      integer :: monte_carlo = 1000
      integer(int64) :: seed = 123456789_int64
      real(dp) :: days_per_year = 365.0_dp
      real(dp) :: discount_rate = 0.0_dp
      logical :: legacy_mode = .true.
      logical :: store_paths = .false.
   end type jdmbs_control

   type, public :: jdmbs_result
      real(dp), allocatable :: call_price(:)
      real(dp), allocatable :: put_price(:)
      real(dp), allocatable :: call_se(:)
      real(dp), allocatable :: put_se(:)
      real(dp), allocatable :: terminal_price(:, :)
      real(dp), allocatable :: paths(:, :, :)
      integer :: jump_events = 0
      integer :: clipped_exponentials = 0
      integer :: status = jdmbs_success
      character(len=160) :: message = 'success'
   end type jdmbs_result

   public :: normal_bs, jdm_bs, jdm_new_bs

contains

   subroutine normal_bs(start_price, mu, sigma, strike, result, control)
      real(dp), intent(in) :: start_price(:), mu(:), sigma(:), strike(:)
      type(jdmbs_result), intent(out) :: result
      type(jdmbs_control), intent(in), optional :: control
      type(jdmbs_control) :: ctrl
      type(rng_state) :: rng
      integer :: i, j, k, n
      real(dp) :: dt, root_dt, w, exponent, terminal

      ctrl = jdmbs_control()
      if (present(control)) ctrl = control
      if (.not. validate_common(start_price, mu, sigma, strike, ctrl, result)) return
      n = size(start_price)
      call initialize_result(result, n, ctrl)
      call rng_seed(rng, ctrl%seed)
      dt = 1.0_dp / ctrl%days_per_year
      root_dt = sqrt(dt)

      do i = 1, n
         do j = 1, ctrl%monte_carlo
            terminal = start_price(i)
            if (ctrl%legacy_mode) then
               w = 0.0_dp
               do k = 0, ctrl%day
                  w = w + rng_normal(rng)
                  exponent = (mu(i) - 0.5_dp * sigma(i)**2) * &
                     (real(k - 1, dp) * dt) + sigma(i) * w
                  terminal = start_price(i) * safe_exp(exponent, result)
                  if (ctrl%store_paths) result%paths(i, j, k) = terminal
               end do
            else
               terminal = start_price(i)
               if (ctrl%store_paths) result%paths(i, j, 0) = terminal
               do k = 1, ctrl%day
                  exponent = (mu(i) - 0.5_dp * sigma(i)**2) * dt + &
                     sigma(i) * root_dt * rng_normal(rng)
                  terminal = terminal * safe_exp(exponent, result)
                  if (ctrl%store_paths) result%paths(i, j, k) = terminal
               end do
            end if
            result%terminal_price(i, j) = terminal
         end do
      end do
      call calculate_prices(strike, ctrl, result)
      call finish_warning(result)
   end subroutine normal_bs

   subroutine jdm_bs(start_price, mu, sigma, lambda, strike, result, control)
      real(dp), intent(in) :: start_price(:), mu(:), sigma(:), lambda
      real(dp), intent(in) :: strike(:)
      type(jdmbs_result), intent(out) :: result
      type(jdmbs_control), intent(in), optional :: control
      type(jdmbs_control) :: ctrl
      type(rng_state) :: rng
      integer :: i, j, k, n
      real(dp) :: dt, root_dt, horizon, rate, next_jump, multiplier
      real(dp) :: w, exponent, terminal, threshold

      ctrl = jdmbs_control()
      if (present(control)) ctrl = control
      if (.not. validate_common(start_price, mu, sigma, strike, ctrl, result)) return
      if (.not. ieee_is_finite(lambda)) then
         call fail(result, jdmbs_nonfinite_input, 'lambda must be finite')
         return
      end if
      if (lambda < 0.0_dp) then
         call fail(result, jdmbs_invalid_argument, 'lambda must be nonnegative')
         return
      end if

      n = size(start_price)
      call initialize_result(result, n, ctrl)
      call rng_seed(rng, ctrl%seed)
      dt = 1.0_dp / ctrl%days_per_year
      root_dt = sqrt(dt)
      horizon = real(ctrl%day, dp) * dt
      if (lambda > 0.0_dp) then
         rate = lambda / horizon
      else
         rate = 0.0_dp
      end if

      do i = 1, n
         do j = 1, ctrl%monte_carlo
            terminal = start_price(i)
            next_jump = next_arrival(rng, rate)
            multiplier = 1.0_dp
            if (ctrl%legacy_mode) then
               w = 0.0_dp
               do k = 0, ctrl%day
                  w = w + rng_normal(rng)
                  threshold = real(k - 1, dp) * dt
                  if (next_jump <= threshold) then
                     multiplier = multiplier * (2.0_dp * rng_uniform(rng))
                     result%jump_events = result%jump_events + 1
                     next_jump = next_arrival_after(rng, rate, next_jump)
                  end if
                  exponent = (mu(i) - 0.5_dp * sigma(i)**2) * threshold + sigma(i) * w
                  terminal = start_price(i) * safe_exp(exponent, result) * multiplier
                  if (ctrl%store_paths) result%paths(i, j, k) = terminal
               end do
            else
               terminal = start_price(i)
               if (ctrl%store_paths) result%paths(i, j, 0) = terminal
               do k = 1, ctrl%day
                  threshold = real(k, dp) * dt
                  do while (next_jump <= threshold)
                     terminal = terminal * (2.0_dp * rng_uniform(rng))
                     result%jump_events = result%jump_events + 1
                     next_jump = next_arrival_after(rng, rate, next_jump)
                  end do
                  exponent = (mu(i) - 0.5_dp * sigma(i)**2) * dt + &
                     sigma(i) * root_dt * rng_normal(rng)
                  terminal = terminal * safe_exp(exponent, result)
                  if (ctrl%store_paths) result%paths(i, j, k) = terminal
               end do
            end if
            result%terminal_price(i, j) = terminal
         end do
      end do
      call calculate_prices(strike, ctrl, result)
      call finish_warning(result)
   end subroutine jdm_bs

   subroutine jdm_new_bs(correlation_matrix, start_price, mu, sigma, lambda, &
      strike, result, control)
      real(dp), intent(in) :: correlation_matrix(:, :)
      real(dp), intent(in) :: start_price(:), mu(:), sigma(:), lambda
      real(dp), intent(in) :: strike(:)
      type(jdmbs_result), intent(out) :: result
      type(jdmbs_control), intent(in), optional :: control
      type(jdmbs_control) :: ctrl
      type(rng_state) :: rng
      integer :: i, j, k, n, event_index
      real(dp) :: dt, root_dt, horizon, rate
      real(dp) :: w, exponent, terminal, threshold, multiplier
      real(dp), allocatable :: event_time(:), event_size(:)
      integer, allocatable :: event_source(:)

      ctrl = jdmbs_control()
      if (present(control)) ctrl = control
      if (.not. validate_common(start_price, mu, sigma, strike, ctrl, result)) return
      n = size(start_price)
      if (size(correlation_matrix, 1) /= n .or. size(correlation_matrix, 2) /= n) then
         call fail(result, jdmbs_invalid_argument, &
            'correlation_matrix must be square with one row per asset')
         return
      end if
      if (.not. all(ieee_is_finite(correlation_matrix)) .or. .not. ieee_is_finite(lambda)) then
         call fail(result, jdmbs_nonfinite_input, 'correlation matrix and lambda must be finite')
         return
      end if
      if (lambda < 0.0_dp .or. any(abs(correlation_matrix) > 1.0_dp)) then
         call fail(result, jdmbs_invalid_argument, &
            'lambda must be nonnegative and transmission entries must lie in [-1,1]')
         return
      end if

      call initialize_result(result, n, ctrl)
      call rng_seed(rng, ctrl%seed)
      dt = 1.0_dp / ctrl%days_per_year
      root_dt = sqrt(dt)
      horizon = real(ctrl%day, dp) * dt
      if (lambda > 0.0_dp) then
         rate = lambda / horizon
      else
         rate = 0.0_dp
      end if

      do j = 1, ctrl%monte_carlo
         call generate_shared_events(rng, rate, horizon, n, event_time, &
            event_source, event_size)
         result%jump_events = result%jump_events + size(event_time)

         if (ctrl%legacy_mode) then
            do i = 1, n
               terminal = start_price(i)
               w = 0.0_dp
               multiplier = 1.0_dp
               event_index = 1
               do k = 0, ctrl%day
                  w = w + rng_normal(rng)
                  threshold = real(k, dp) * dt
                  if (event_index <= size(event_time)) then
                     if (event_time(event_index) <= threshold) then
                        multiplier = multiplier * (1.0_dp + &
                           event_size(event_index) * &
                           correlation_matrix(event_source(event_index), i))
                        event_index = event_index + 1
                     end if
                  end if
                  exponent = (mu(i) - 0.5_dp * sigma(i)**2) * threshold + sigma(i) * w
                  terminal = start_price(i) * safe_exp(exponent, result) * multiplier
                  if (ctrl%store_paths) result%paths(i, j, k) = terminal
               end do
               result%terminal_price(i, j) = terminal
            end do
         else
            result%terminal_price(:, j) = start_price
            if (ctrl%store_paths) result%paths(:, j, 0) = start_price
            event_index = 1
            do k = 1, ctrl%day
               threshold = real(k, dp) * dt
               do while (event_index <= size(event_time))
                  if (event_time(event_index) > threshold) exit
                  do i = 1, n
                     result%terminal_price(i, j) = result%terminal_price(i, j) * &
                        (1.0_dp + event_size(event_index) * &
                        correlation_matrix(event_source(event_index), i))
                  end do
                  event_index = event_index + 1
               end do
               do i = 1, n
                  exponent = (mu(i) - 0.5_dp * sigma(i)**2) * dt + &
                     sigma(i) * root_dt * rng_normal(rng)
                  result%terminal_price(i, j) = result%terminal_price(i, j) * &
                     safe_exp(exponent, result)
                  if (ctrl%store_paths) then
                     result%paths(i, j, k) = result%terminal_price(i, j)
                  end if
               end do
            end do
         end if
      end do
      call calculate_prices(strike, ctrl, result)
      call finish_warning(result)
   end subroutine jdm_new_bs

   logical function validate_common(start_price, mu, sigma, strike, ctrl, result)
      real(dp), intent(in) :: start_price(:), mu(:), sigma(:), strike(:)
      type(jdmbs_control), intent(in) :: ctrl
      type(jdmbs_result), intent(out) :: result
      integer :: n

      validate_common = .false.
      result%status = jdmbs_success
      result%message = 'success'
      n = size(start_price)
      if (n < 1 .or. size(mu) /= n .or. size(sigma) /= n .or. size(strike) /= n) then
         call fail(result, jdmbs_invalid_argument, &
            'start_price, mu, sigma, and strike must have equal nonzero size')
         return
      end if
      if (ctrl%day <= 0 .or. ctrl%monte_carlo <= 0 .or. &
          ctrl%days_per_year <= 0.0_dp) then
         call fail(result, jdmbs_invalid_argument, &
            'day, monte_carlo, and days_per_year must be positive')
         return
      end if
      if (.not. all(ieee_is_finite(start_price)) .or. &
          .not. all(ieee_is_finite(mu)) .or. &
          .not. all(ieee_is_finite(sigma)) .or. &
          .not. all(ieee_is_finite(strike)) .or. &
          .not. ieee_is_finite(ctrl%discount_rate) .or. &
          .not. ieee_is_finite(ctrl%days_per_year)) then
         call fail(result, jdmbs_nonfinite_input, 'all numerical inputs must be finite')
         return
      end if
      if (any(start_price <= 0.0_dp) .or. any(sigma < 0.0_dp) .or. &
          any(strike < 0.0_dp)) then
         call fail(result, jdmbs_invalid_argument, &
            'start prices must be positive; volatilities and strikes nonnegative')
         return
      end if
      validate_common = .true.
   end function validate_common

   subroutine initialize_result(result, n, ctrl)
      type(jdmbs_result), intent(inout) :: result
      integer, intent(in) :: n
      type(jdmbs_control), intent(in) :: ctrl

      allocate(result%call_price(n), result%put_price(n))
      allocate(result%call_se(n), result%put_se(n))
      allocate(result%terminal_price(n, ctrl%monte_carlo))
      result%call_price = 0.0_dp
      result%put_price = 0.0_dp
      result%call_se = 0.0_dp
      result%put_se = 0.0_dp
      result%terminal_price = 0.0_dp
      if (ctrl%store_paths) then
         allocate(result%paths(n, ctrl%monte_carlo, 0:ctrl%day))
         result%paths = 0.0_dp
      end if
   end subroutine initialize_result

   subroutine calculate_prices(strike, ctrl, result)
      real(dp), intent(in) :: strike(:)
      type(jdmbs_control), intent(in) :: ctrl
      type(jdmbs_result), intent(inout) :: result
      integer :: i, m
      real(dp) :: discount
      real(dp), allocatable :: payoff(:)

      m = ctrl%monte_carlo
      discount = exp(-ctrl%discount_rate * real(ctrl%day, dp) / ctrl%days_per_year)
      allocate(payoff(m))
      do i = 1, size(strike)
         payoff = max(result%terminal_price(i, :) - strike(i), 0.0_dp)
         result%call_price(i) = discount * sum(payoff) / real(m, dp)
         result%call_se(i) = discount * standard_error(payoff)
         payoff = max(strike(i) - result%terminal_price(i, :), 0.0_dp)
         result%put_price(i) = discount * sum(payoff) / real(m, dp)
         result%put_se(i) = discount * standard_error(payoff)
      end do
   end subroutine calculate_prices

   pure function standard_error(x) result(se)
      real(dp), intent(in) :: x(:)
      real(dp) :: se, mean_x

      if (size(x) <= 1) then
         se = 0.0_dp
      else
         mean_x = sum(x) / real(size(x), dp)
         se = sqrt(sum((x - mean_x)**2) / real(size(x) - 1, dp)) / &
            sqrt(real(size(x), dp))
      end if
   end function standard_error

   subroutine generate_shared_events(rng, rate, horizon, n_assets, times, &
      sources, sizes)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: rate, horizon
      integer, intent(in) :: n_assets
      real(dp), allocatable, intent(out) :: times(:), sizes(:)
      integer, allocatable, intent(out) :: sources(:)
      real(dp), allocatable :: work_times(:), work_sizes(:), grown_real(:)
      integer, allocatable :: work_sources(:), grown_integer(:)
      real(dp) :: current
      integer :: count, capacity

      capacity = 16
      allocate(work_times(capacity), work_sizes(capacity), work_sources(capacity))
      count = 0
      current = next_arrival(rng, rate)
      do while (current <= horizon)
         count = count + 1
         if (count > capacity) then
            allocate(grown_real(2 * capacity))
            grown_real(1:capacity) = work_times
            call move_alloc(grown_real, work_times)
            allocate(grown_real(2 * capacity))
            grown_real(1:capacity) = work_sizes
            call move_alloc(grown_real, work_sizes)
            allocate(grown_integer(2 * capacity))
            grown_integer(1:capacity) = work_sources
            call move_alloc(grown_integer, work_sources)
            capacity = 2 * capacity
         end if
         work_times(count) = current
         work_sources(count) = rng_integer(rng, 1, n_assets)
         work_sizes(count) = 2.0_dp * rng_uniform(rng) - 1.0_dp
         current = next_arrival_after(rng, rate, current)
      end do
      allocate(times(count), sizes(count), sources(count))
      if (count > 0) then
         times = work_times(1:count)
         sizes = work_sizes(1:count)
         sources = work_sources(1:count)
      end if
   end subroutine generate_shared_events

   function next_arrival(rng, rate) result(value)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: rate
      real(dp) :: value

      if (rate <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = -log(max(1.0_dp - rng_uniform(rng), tiny(1.0_dp))) / rate
      end if
   end function next_arrival

   function next_arrival_after(rng, rate, current) result(value)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: rate, current
      real(dp) :: value

      if (rate <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = current - log(max(1.0_dp - rng_uniform(rng), tiny(1.0_dp))) / rate
      end if
   end function next_arrival_after

   function safe_exp(x, result) result(value)
      real(dp), intent(in) :: x
      type(jdmbs_result), intent(inout) :: result
      real(dp) :: value
      real(dp), parameter :: upper = log(huge(1.0_dp)) - 2.0_dp
      real(dp), parameter :: lower = log(tiny(1.0_dp)) + 2.0_dp

      if (x > upper) then
         value = exp(upper)
         result%clipped_exponentials = result%clipped_exponentials + 1
      else if (x < lower) then
         value = 0.0_dp
         result%clipped_exponentials = result%clipped_exponentials + 1
      else
         value = exp(x)
      end if
   end function safe_exp

   subroutine finish_warning(result)
      type(jdmbs_result), intent(inout) :: result

      if (result%clipped_exponentials > 0) then
         result%status = jdmbs_numerical_warning
         write (result%message, '(a,i0,a)') 'success with ', &
            result%clipped_exponentials, ' clipped exponentials'
      end if
   end subroutine finish_warning

   subroutine fail(result, status, message)
      type(jdmbs_result), intent(out) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message

      result%status = status
      result%message = message
   end subroutine fail

end module jdmbs_model
