! SmithWilsonYieldCurve modern Fortran translation
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-only
!
! The algorithms and public behavior in this module are derived from
! SmithWilsonYieldCurve 1.1.1 by Phil Joubert, distributed under GPL-3.

module smith_wilson
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use smith_wilson_kinds, only : dp
   use smith_wilson_linalg, only : solve_linear_system
   implicit none
   private

   public :: dp

   integer, parameter, public :: sw_success = 0
   integer, parameter, public :: sw_invalid_argument = 1
   integer, parameter, public :: sw_dimension_error = 2
   integer, parameter, public :: sw_singular_system = 3
   integer, parameter, public :: sw_unknown_instrument = 4

   integer, parameter, public :: sw_libor = 1
   integer, parameter, public :: sw_swap = 2
   integer, parameter, public :: sw_bond = 3

   type, public :: market_instrument
      integer :: instrument_type = 0
      real(dp) :: tenor = 0.0_dp
      real(dp) :: frequency = 0.0_dp
      real(dp) :: rate = 0.0_dp
      real(dp) :: price = 1.0_dp
   end type market_instrument

   type, public :: smith_wilson_curve
      real(dp) :: ufr = 0.0_dp
      real(dp) :: alpha = 0.0_dp
      real(dp), allocatable :: times(:)
      real(dp), allocatable :: cashflows(:, :)
      real(dp), allocatable :: market_values(:)
      real(dp), allocatable :: xi(:)
      logical :: fitted = .false.
   contains
      procedure, private :: curve_discount_scalar
      procedure, private :: curve_discount_vector
      generic, public :: discount => curve_discount_scalar, curve_discount_vector
      procedure, private :: curve_compound_kernel_scalar
      procedure, private :: curve_compound_kernel_vector
      generic, public :: compound_kernel => curve_compound_kernel_scalar, curve_compound_kernel_vector
      procedure, private :: curve_continuous_spot_scalar
      procedure, private :: curve_continuous_spot_vector
      generic, public :: continuous_spot => curve_continuous_spot_scalar, curve_continuous_spot_vector
      procedure, public :: repriced_values => curve_repriced_values
   end type smith_wilson_curve

   public :: wilson_function
   public :: create_kernel_matrix
   public :: fit_kernel_weights
   public :: fit_smith_wilson_curve
   public :: fit_smith_wilson_curve_to_instruments
   public :: create_time_vector
   public :: create_cashflow_matrix
   public :: create_market_value_vector
   public :: get_instrument_times
   public :: get_instrument_cashflows
   public :: make_market_instrument

contains

   pure elemental function wilson_function(t, u, ufr, alpha) result(w)
      real(dp), intent(in) :: t, u, ufr, alpha
      real(dp) :: w
      real(dp) :: minimum_time, maximum_time, hyperbolic_term

      minimum_time = min(t, u)
      maximum_time = max(t, u)

      ! This is algebraically equal to
      ! exp(-ufr*(t+u)) * (alpha*min(t,u) - sinh(alpha*min(t,u)) /
      ! exp(alpha*max(t,u))), but avoids overflow in sinh for large terms.
      hyperbolic_term = 0.5_dp * (exp(-alpha * (maximum_time - minimum_time)) - &
                                      exp(-alpha * (maximum_time + minimum_time)))
      w = exp(-ufr * (t + u)) * (alpha * minimum_time - hyperbolic_term)
   end function wilson_function


   subroutine create_kernel_matrix(times, kernel_matrix, ufr, alpha, info, message)
      real(dp), intent(in) :: times(:)
      real(dp), allocatable, intent(out) :: kernel_matrix(:, :)
      real(dp), intent(in) :: ufr, alpha
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer :: i, j, local_info, n

      local_info = sw_success
      call set_message(message, 'ok')

      n = size(times)
      if (n < 1 .or. any(times < 0.0_dp)) then
         local_info = sw_invalid_argument
         call set_message(message, 'times must be nonnegative and nonempty')
         allocate(kernel_matrix(0, 0))
         call set_info(info, local_info)
         return
      end if
      if (alpha <= 0.0_dp) then
         local_info = sw_invalid_argument
         call set_message(message, 'alpha must be greater than zero')
         allocate(kernel_matrix(0, 0))
         call set_info(info, local_info)
         return
      end if

      allocate(kernel_matrix(n, n))
      do j = 1, n
         do i = 1, j
            kernel_matrix(i, j) = wilson_function(times(i), times(j), ufr, alpha)
            kernel_matrix(j, i) = kernel_matrix(i, j)
         end do
      end do

      if (ufr < 0.0_dp) call set_message(message, 'warning: ufr is negative')
      call set_info(info, local_info)
   end subroutine create_kernel_matrix


   subroutine fit_kernel_weights(cashflow_matrix, kernel_matrix, market_values, &
                                 base_zero_values, xi, info, message)
      real(dp), intent(in) :: cashflow_matrix(:, :)
      real(dp), intent(in) :: kernel_matrix(:, :)
      real(dp), intent(in) :: market_values(:)
      real(dp), intent(in) :: base_zero_values(:)
      real(dp), allocatable, intent(out) :: xi(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: system_matrix(:, :), rhs(:), solution(:)
      integer :: local_info, n_instruments, n_times, solve_info

      local_info = sw_success
      call set_message(message, 'ok')

      n_instruments = size(cashflow_matrix, 1)
      n_times = size(cashflow_matrix, 2)
      if (n_instruments < 1 .or. n_times < 1 .or. &
          size(kernel_matrix, 1) /= n_times .or. size(kernel_matrix, 2) /= n_times .or. &
          size(market_values) /= n_instruments .or. size(base_zero_values) /= n_times) then
         local_info = sw_dimension_error
         call set_message(message, 'inconsistent cashflow, kernel, market-value, or base-vector dimensions')
         allocate(xi(0))
         call set_info(info, local_info)
         return
      end if

      allocate(system_matrix(n_instruments, n_instruments), rhs(n_instruments), &
               solution(n_instruments))
      system_matrix = matmul(matmul(cashflow_matrix, kernel_matrix), transpose(cashflow_matrix))
      rhs = market_values - matmul(cashflow_matrix, base_zero_values)

      call solve_linear_system(system_matrix, rhs, solution, solve_info)
      if (solve_info /= 0) then
         local_info = sw_singular_system
         call set_message(message, 'the Smith-Wilson calibration system is singular or numerically rank deficient')
         allocate(xi(0))
         call set_info(info, local_info)
         return
      end if

      allocate(xi(n_instruments))
      xi = solution
      call set_info(info, local_info)
   end subroutine fit_kernel_weights


   subroutine fit_smith_wilson_curve(times, cashflow_matrix, market_values, ufr, alpha, &
                                     curve, info, message)
      real(dp), intent(in) :: times(:)
      real(dp), intent(in) :: cashflow_matrix(:, :)
      real(dp), intent(in) :: market_values(:)
      real(dp), intent(in) :: ufr, alpha
      type(smith_wilson_curve), intent(out) :: curve
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: kernel_matrix(:, :), base_zero_values(:), weights(:)
      integer :: local_info
      character(len=256) :: local_message

      call clear_curve(curve)
      local_info = sw_success
      local_message = 'ok'

      if (size(times) < 1 .or. any(times < 0.0_dp)) then
         local_info = sw_invalid_argument
         local_message = 'times must be nonnegative and nonempty'
      else if (alpha <= 0.0_dp) then
         local_info = sw_invalid_argument
         local_message = 'alpha must be greater than zero'
      else if (size(cashflow_matrix, 2) /= size(times) .or. &
               size(cashflow_matrix, 1) /= size(market_values) .or. &
               size(market_values) < 1) then
         local_info = sw_dimension_error
         local_message = 'inconsistent times, cashflow-matrix, and market-value dimensions'
      else if (.not. all_finite_vector(times) .or. .not. all_finite_matrix(cashflow_matrix) .or. &
               .not. all_finite_vector(market_values) .or. .not. is_finite(ufr) .or. &
               .not. is_finite(alpha)) then
         local_info = sw_invalid_argument
         local_message = 'all numerical inputs must be finite'
      end if

      if (local_info /= sw_success) then
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      call create_kernel_matrix(times, kernel_matrix, ufr, alpha, local_info, local_message)
      if (local_info /= sw_success) then
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      allocate(base_zero_values(size(times)))
      base_zero_values = exp(-ufr * times)
      call fit_kernel_weights(cashflow_matrix, kernel_matrix, market_values, &
                              base_zero_values, weights, local_info, local_message)
      if (local_info /= sw_success) then
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      curve%ufr = ufr
      curve%alpha = alpha
      curve%times = times
      curve%cashflows = cashflow_matrix
      curve%market_values = market_values
      curve%xi = weights
      curve%fitted = .true.

      if (ufr < 0.0_dp) local_message = 'warning: ufr is negative'
      call set_info(info, sw_success)
      call set_message(message, trim(local_message))
   end subroutine fit_smith_wilson_curve


   subroutine fit_smith_wilson_curve_to_instruments(instruments, ufr, alpha, curve, info, message)
      type(market_instrument), intent(in) :: instruments(:)
      real(dp), intent(in) :: ufr, alpha
      type(smith_wilson_curve), intent(out) :: curve
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: times(:), cashflow_matrix(:, :), market_values(:)
      integer :: local_info
      character(len=256) :: local_message

      call create_time_vector(instruments, times, local_info, local_message)
      if (local_info /= sw_success) then
         call clear_curve(curve)
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      call create_cashflow_matrix(instruments, times, cashflow_matrix, local_info, local_message)
      if (local_info /= sw_success) then
         call clear_curve(curve)
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      call create_market_value_vector(instruments, market_values, local_info, local_message)
      if (local_info /= sw_success) then
         call clear_curve(curve)
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      call fit_smith_wilson_curve(times, cashflow_matrix, market_values, ufr, alpha, &
                                  curve, local_info, local_message)
      call set_info(info, local_info)
      call set_message(message, trim(local_message))
   end subroutine fit_smith_wilson_curve_to_instruments


   subroutine make_market_instrument(type_name, tenor, rate, instrument, info, message, frequency, price)
      character(len=*), intent(in) :: type_name
      real(dp), intent(in) :: tenor, rate
      type(market_instrument), intent(out) :: instrument
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message
      real(dp), intent(in), optional :: frequency, price

      character(len=:), allocatable :: normalized_name
      integer :: local_info

      instrument = market_instrument()
      instrument%tenor = tenor
      instrument%rate = rate
      if (present(frequency)) instrument%frequency = frequency
      if (present(price)) instrument%price = price

      normalized_name = uppercase(trim(adjustl(type_name)))
      select case (normalized_name)
      case ('LIBOR')
         instrument%instrument_type = sw_libor
      case ('SWAP')
         instrument%instrument_type = sw_swap
      case ('BOND')
         instrument%instrument_type = sw_bond
      case default
         local_info = sw_unknown_instrument
         call set_info(info, local_info)
         call set_message(message, 'unknown instrument type: '//trim(type_name))
         return
      end select

      call validate_instrument(instrument, local_info, message)
      call set_info(info, local_info)
   end subroutine make_market_instrument


   subroutine create_time_vector(instruments, times, info, message)
      type(market_instrument), intent(in) :: instruments(:)
      real(dp), allocatable, intent(out) :: times(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: schedule(:), all_times(:), unique_times(:)
      integer :: i, local_info, n_all, n_unique, position
      character(len=256) :: local_message

      call set_message(message, 'ok')
      if (size(instruments) < 1) then
         allocate(times(0))
         call set_info(info, sw_invalid_argument)
         call set_message(message, 'at least one instrument is required')
         return
      end if

      n_all = 0
      do i = 1, size(instruments)
         call get_instrument_times(instruments(i), schedule, local_info, local_message)
         if (local_info /= sw_success) then
            allocate(times(0))
            call set_info(info, local_info)
            call set_message(message, 'instrument '//integer_string(i)//': '//trim(local_message))
            return
         end if
         n_all = n_all + size(schedule)
      end do

      allocate(all_times(n_all))
      position = 1
      do i = 1, size(instruments)
         call get_instrument_times(instruments(i), schedule, local_info, local_message)
         all_times(position:position + size(schedule) - 1) = schedule
         position = position + size(schedule)
      end do

      call sort_real(all_times)
      allocate(unique_times(n_all))
      n_unique = 0
      do i = 1, n_all
         if (n_unique == 0) then
            n_unique = 1
            unique_times(n_unique) = all_times(i)
         else if (.not. nearly_equal(all_times(i), unique_times(n_unique))) then
            n_unique = n_unique + 1
            unique_times(n_unique) = all_times(i)
         end if
      end do

      allocate(times(n_unique))
      times = unique_times(1:n_unique)
      call set_info(info, sw_success)
   end subroutine create_time_vector


   subroutine create_cashflow_matrix(instruments, times, cashflow_matrix, info, message)
      type(market_instrument), intent(in) :: instruments(:)
      real(dp), intent(in) :: times(:)
      real(dp), allocatable, intent(out) :: cashflow_matrix(:, :)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: schedule_times(:), schedule_cashflows(:)
      integer :: i, j, local_info, column
      character(len=256) :: local_message

      call set_message(message, 'ok')
      if (size(instruments) < 1 .or. size(times) < 1) then
         allocate(cashflow_matrix(0, 0))
         call set_info(info, sw_invalid_argument)
         call set_message(message, 'instruments and times must be nonempty')
         return
      end if

      allocate(cashflow_matrix(size(instruments), size(times)))
      cashflow_matrix = 0.0_dp

      do i = 1, size(instruments)
         call get_instrument_cashflows(instruments(i), schedule_times, schedule_cashflows, &
                                       local_info, local_message)
         if (local_info /= sw_success) then
            deallocate(cashflow_matrix)
            allocate(cashflow_matrix(0, 0))
            call set_info(info, local_info)
            call set_message(message, 'instrument '//integer_string(i)//': '//trim(local_message))
            return
         end if

         do j = 1, size(schedule_times)
            column = find_matching_time(schedule_times(j), times)
            if (column == 0) then
               deallocate(cashflow_matrix)
               allocate(cashflow_matrix(0, 0))
               call set_info(info, sw_dimension_error)
               call set_message(message, 'the supplied time vector omits an instrument cashflow time')
               return
            end if
            cashflow_matrix(i, column) = cashflow_matrix(i, column) + schedule_cashflows(j)
         end do
      end do

      call set_info(info, sw_success)
   end subroutine create_cashflow_matrix


   subroutine create_market_value_vector(instruments, market_values, info, message)
      type(market_instrument), intent(in) :: instruments(:)
      real(dp), allocatable, intent(out) :: market_values(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer :: i, local_info
      character(len=256) :: local_message

      call set_message(message, 'ok')
      if (size(instruments) < 1) then
         allocate(market_values(0))
         call set_info(info, sw_invalid_argument)
         call set_message(message, 'at least one instrument is required')
         return
      end if

      allocate(market_values(size(instruments)))
      do i = 1, size(instruments)
         call validate_instrument(instruments(i), local_info, local_message)
         if (local_info /= sw_success) then
            deallocate(market_values)
            allocate(market_values(0))
            call set_info(info, local_info)
            call set_message(message, 'instrument '//integer_string(i)//': '//trim(local_message))
            return
         end if
         if (instruments(i)%instrument_type == sw_bond) then
            market_values(i) = instruments(i)%price
         else
            market_values(i) = 1.0_dp
         end if
      end do

      call set_info(info, sw_success)
   end subroutine create_market_value_vector


   subroutine get_instrument_times(instrument, times, info, message)
      type(market_instrument), intent(in) :: instrument
      real(dp), allocatable, intent(out) :: times(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp) :: dcf, tx, tolerance
      integer :: i, local_info, n_payments
      character(len=256) :: local_message

      call validate_instrument(instrument, local_info, local_message)
      if (local_info /= sw_success) then
         allocate(times(0))
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      select case (instrument%instrument_type)
      case (sw_libor)
         allocate(times(1))
         times(1) = instrument%tenor

      case (sw_swap)
         n_payments = nint(instrument%tenor * instrument%frequency)
         allocate(times(n_payments))
         do i = 1, n_payments
            times(i) = real(i, dp) / instrument%frequency
         end do

      case (sw_bond)
         dcf = 1.0_dp / instrument%frequency
         tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, instrument%tenor)
         n_payments = 0
         tx = instrument%tenor
         do while (tx > tolerance)
            n_payments = n_payments + 1
            tx = tx - dcf
         end do
         if (n_payments < 1) then
            allocate(times(0))
            call set_info(info, sw_invalid_argument)
            call set_message(message, 'no bond cashflows were calculated; check tenor and frequency')
            return
         end if
         allocate(times(n_payments))
         tx = instrument%tenor
         do i = n_payments, 1, -1
            times(i) = tx
            tx = tx - dcf
         end do

      case default
         allocate(times(0))
         call set_info(info, sw_unknown_instrument)
         call set_message(message, 'unknown instrument type')
         return
      end select

      call set_info(info, sw_success)
      call set_message(message, 'ok')
   end subroutine get_instrument_times


   subroutine get_instrument_cashflows(instrument, times, cashflows, info, message)
      type(market_instrument), intent(in) :: instrument
      real(dp), allocatable, intent(out) :: times(:)
      real(dp), allocatable, intent(out) :: cashflows(:)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      integer :: local_info
      character(len=256) :: local_message

      call get_instrument_times(instrument, times, local_info, local_message)
      if (local_info /= sw_success) then
         allocate(cashflows(0))
         call set_info(info, local_info)
         call set_message(message, trim(local_message))
         return
      end if

      allocate(cashflows(size(times)))
      select case (instrument%instrument_type)
      case (sw_libor)
         cashflows(1) = 1.0_dp + instrument%rate * instrument%tenor
      case (sw_swap, sw_bond)
         cashflows = instrument%rate / instrument%frequency
         cashflows(size(cashflows)) = cashflows(size(cashflows)) + 1.0_dp
      case default
         cashflows = 0.0_dp
         call set_info(info, sw_unknown_instrument)
         call set_message(message, 'unknown instrument type')
         return
      end select

      call set_info(info, sw_success)
      call set_message(message, 'ok')
   end subroutine get_instrument_cashflows


   function curve_discount_scalar(self, t) result(price)
      class(smith_wilson_curve), intent(in) :: self
      real(dp), intent(in) :: t
      real(dp) :: price
      real(dp), allocatable :: kernel(:)

      if (.not. self%fitted .or. t < 0.0_dp .or. .not. is_finite(t)) then
         price = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      kernel = self%compound_kernel(t)
      price = exp(-self%ufr * t) + dot_product(self%xi, kernel)
   end function curve_discount_scalar


   function curve_discount_vector(self, t) result(price)
      class(smith_wilson_curve), intent(in) :: self
      real(dp), intent(in) :: t(:)
      real(dp), allocatable :: price(:)
      real(dp), allocatable :: kernel(:, :)
      integer :: i

      allocate(price(size(t)))
      if (.not. self%fitted .or. any(t < 0.0_dp) .or. .not. all_finite_vector(t)) then
         price = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      kernel = self%compound_kernel(t)
      do i = 1, size(t)
         price(i) = exp(-self%ufr * t(i)) + dot_product(self%xi, kernel(:, i))
      end do
   end function curve_discount_vector


   function curve_compound_kernel_scalar(self, t) result(kernel)
      class(smith_wilson_curve), intent(in) :: self
      real(dp), intent(in) :: t
      real(dp), allocatable :: kernel(:)
      real(dp), allocatable :: values(:)
      integer :: j

      if (.not. self%fitted) then
         allocate(kernel(0))
         return
      end if

      allocate(values(size(self%times)))
      do j = 1, size(self%times)
         values(j) = wilson_function(t, self%times(j), self%ufr, self%alpha)
      end do
      kernel = matmul(self%cashflows, values)
   end function curve_compound_kernel_scalar


   function curve_compound_kernel_vector(self, t) result(kernel)
      class(smith_wilson_curve), intent(in) :: self
      real(dp), intent(in) :: t(:)
      real(dp), allocatable :: kernel(:, :)
      real(dp), allocatable :: values(:, :)
      integer :: i, j

      if (.not. self%fitted) then
         allocate(kernel(0, 0))
         return
      end if

      allocate(values(size(self%times), size(t)))
      do i = 1, size(t)
         do j = 1, size(self%times)
            values(j, i) = wilson_function(t(i), self%times(j), self%ufr, self%alpha)
         end do
      end do
      kernel = matmul(self%cashflows, values)
   end function curve_compound_kernel_vector


   function curve_continuous_spot_scalar(self, t) result(rate)
      class(smith_wilson_curve), intent(in) :: self
      real(dp), intent(in) :: t
      real(dp) :: rate, price

      if (t <= 0.0_dp) then
         rate = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      price = self%discount(t)
      if (.not. is_finite(price) .or. price <= 0.0_dp) then
         rate = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         rate = -log(price) / t
      end if
   end function curve_continuous_spot_scalar


   function curve_continuous_spot_vector(self, t) result(rate)
      class(smith_wilson_curve), intent(in) :: self
      real(dp), intent(in) :: t(:)
      real(dp), allocatable :: rate(:), price(:)
      integer :: i

      allocate(rate(size(t)))
      price = self%discount(t)
      do i = 1, size(t)
         if (t(i) <= 0.0_dp .or. .not. is_finite(price(i)) .or. price(i) <= 0.0_dp) then
            rate(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            rate(i) = -log(price(i)) / t(i)
         end if
      end do
   end function curve_continuous_spot_vector


   function curve_repriced_values(self) result(values)
      class(smith_wilson_curve), intent(in) :: self
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: zero_prices(:)

      if (.not. self%fitted) then
         allocate(values(0))
         return
      end if
      zero_prices = self%discount(self%times)
      values = matmul(self%cashflows, zero_prices)
   end function curve_repriced_values


   subroutine validate_instrument(instrument, info, message)
      type(market_instrument), intent(in) :: instrument
      integer, intent(out) :: info
      character(len=*), intent(out), optional :: message

      real(dp) :: payment_count

      info = sw_success
      call set_message(message, 'ok')

      if (.not. is_finite(instrument%tenor) .or. instrument%tenor <= 0.0_dp .or. &
          .not. is_finite(instrument%rate) .or. .not. is_finite(instrument%price)) then
         info = sw_invalid_argument
         call set_message(message, 'tenor must be positive and all instrument values must be finite')
         return
      end if

      select case (instrument%instrument_type)
      case (sw_libor)
         continue
      case (sw_swap)
         if (.not. is_finite(instrument%frequency) .or. instrument%frequency <= 0.0_dp) then
            info = sw_invalid_argument
            call set_message(message, 'swap frequency must be greater than zero')
            return
         end if
         payment_count = instrument%tenor * instrument%frequency
         if (payment_count < 1.0_dp .or. &
             abs(payment_count - real(nint(payment_count), dp)) > &
             100.0_dp * epsilon(1.0_dp) * max(1.0_dp, payment_count)) then
            info = sw_invalid_argument
            call set_message(message, 'swap tenor times frequency must be a positive integer')
            return
         end if
      case (sw_bond)
         if (.not. is_finite(instrument%frequency) .or. instrument%frequency <= 0.0_dp) then
            info = sw_invalid_argument
            call set_message(message, 'bond frequency must be greater than zero')
            return
         end if
      case default
         info = sw_unknown_instrument
         call set_message(message, 'unknown instrument type')
      end select
   end subroutine validate_instrument


   subroutine clear_curve(curve)
      type(smith_wilson_curve), intent(out) :: curve

      curve%ufr = 0.0_dp
      curve%alpha = 0.0_dp
      curve%fitted = .false.
   end subroutine clear_curve


   pure logical function is_finite(x) result(finite)
      real(dp), intent(in) :: x

      finite = ieee_is_finite(x)
   end function is_finite


   pure logical function all_finite_vector(x) result(finite)
      real(dp), intent(in) :: x(:)
      integer :: i

      finite = .true.
      do i = 1, size(x)
         if (.not. is_finite(x(i))) then
            finite = .false.
            return
         end if
      end do
   end function all_finite_vector


   pure logical function all_finite_matrix(x) result(finite)
      real(dp), intent(in) :: x(:, :)
      integer :: i, j

      finite = .true.
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (.not. is_finite(x(i, j))) then
               finite = .false.
               return
            end if
         end do
      end do
   end function all_finite_matrix


   pure logical function nearly_equal(x, y) result(equal)
      real(dp), intent(in) :: x, y
      real(dp) :: scale

      scale = max(1.0_dp, abs(x), abs(y))
      equal = abs(x - y) <= 100.0_dp * epsilon(1.0_dp) * scale
   end function nearly_equal


   integer function find_matching_time(value, times) result(index)
      real(dp), intent(in) :: value, times(:)
      integer :: i

      index = 0
      do i = 1, size(times)
         if (nearly_equal(value, times(i))) then
            index = i
            return
         end if
      end do
   end function find_matching_time


   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j

      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real


   pure function uppercase(text) result(converted)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: converted
      integer :: code, i

      converted = text
      do i = 1, len(text)
         code = iachar(converted(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) then
            converted(i:i) = achar(code - iachar('a') + iachar('A'))
         end if
      end do
   end function uppercase


   function integer_string(value) result(text)
      integer, intent(in) :: value
      character(len=:), allocatable :: text
      character(len=32) :: buffer

      write(buffer, '(i0)') value
      text = trim(buffer)
   end function integer_string


   subroutine set_info(info, value)
      integer, intent(out), optional :: info
      integer, intent(in) :: value

      if (present(info)) info = value
   end subroutine set_info


   subroutine set_message(message, value)
      character(len=*), intent(out), optional :: message
      character(len=*), intent(in) :: value

      if (present(message)) message = value
   end subroutine set_message

end module smith_wilson
