! SPDX-License-Identifier: MIT
! SPDX-FileComment: Basic time-series adapters corresponding to selected R stats functions.
module r_stats_time_series
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_status, only: r_invalid_input, r_ok
   use r_stats_types, only: acf_result_t
   use r_time_series, only: r_arma_autocorrelation, r_arma_simulate, r_autocorrelation, &
                            r_autocovariance, r_convolution_filter, &
                            r_cross_correlation, r_cross_covariance, r_partial_autocorrelation, &
                            r_recursive_filter
   implicit none
   private

   public :: acf, arima_sim, arma_acf, ccf, filter, filter_linear, filter_recursive, pacf

contains

   pure function arma_acf(ar, ma, lag_max) result(values)
      !! Computes theoretical autocorrelations for a stationary ARMA model.
      real(dp), intent(in), optional :: ar(:) !! Autoregressive coefficients in increasing lag order.
      real(dp), intent(in), optional :: ma(:) !! Moving-average coefficients in increasing lag order.
      integer, intent(in), optional :: lag_max !! Largest lag; defaults to `max(p, q + 1)`.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: ar_values(:), ma_values(:), core_values(:)
      integer :: local_status, maximum_lag

      if (present(ar)) then
         ar_values = ar
      else
         allocate (ar_values(0))
      end if
      if (present(ma)) then
         ma_values = ma
      else
         allocate (ma_values(0))
      end if
      maximum_lag = max(size(ar_values), size(ma_values) + 1)
      if (present(lag_max)) maximum_lag = lag_max
      call r_arma_autocorrelation(ar_values, ma_values, maximum_lag, core_values, local_status)
      if (local_status == r_ok) values = core_values
   end function arma_acf

   pure function arima_sim(innovations, ar, ma, start_innovations, difference_order) &
      result(values)
      !! Simulates a possibly integrated nonseasonal ARIMA model from explicit innovations.
      real(dp), intent(in) :: innovations(:) !! Innovations retained in the returned series.
      real(dp), intent(in), optional :: ar(:) !! Autoregressive coefficients in increasing lag order.
      real(dp), intent(in), optional :: ma(:) !! Moving-average coefficients in increasing lag order.
      real(dp), intent(in), optional :: start_innovations(:) !! Earlier burn-in innovations.
      integer, intent(in), optional :: difference_order !! Nonnegative number of differences; defaults to zero.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: integrated(:)
      integer :: d, i, local_status, pass

      d = optval(difference_order, 0)
      if (d < 0) then
         allocate (values(0))
         return
      end if
      call r_arma_simulate(innovations, values, ar, ma, start_innovations, local_status)
      if (local_status /= r_ok) then
         if (allocated(values)) deallocate (values)
         allocate (values(0))
         return
      end if
      do pass = 1, d
         allocate (integrated(size(values) + 1))
         integrated(1) = 0.0_dp
         do i = 1, size(values)
            integrated(i + 1) = integrated(i) + values(i)
         end do
         call move_alloc(integrated, values)
      end do
   end function arima_sim

   pure function filter(x, coefficients, sides, method, initial) result(filtered)
      !! Applies an R-compatible convolution or recursive filter to a univariate series.
      real(dp), intent(in) :: x(:) !! Input time series.
      real(dp), intent(in) :: coefficients(:) !! Filter coefficients in increasing lag order.
      integer, intent(in), optional :: sides !! One for past values or two for a centered filter; defaults to two.
      character(len=*), intent(in), optional :: method !! `convolution` or `recursive`; defaults to convolution.
      real(dp), intent(in), optional :: initial(:) !! Recursive states `y(0), y(-1), ...`; defaults to zeros.
      real(dp), allocatable :: filtered(:)
      character(len=:), allocatable :: method_name

      method_name = "convolution"
      if (present(method)) method_name = trim(adjustl(method))
      select case (method_name)
      case ("convolution")
         filtered = filter_linear(x, coefficients, sides)
      case ("recursive")
         filtered = filter_recursive(x, coefficients, initial)
      case default
         allocate (filtered(0))
      end select
   end function filter

   pure function filter_linear(x, coefficients, sides) result(filtered)
      !! Applies a finite convolution filter and returns NaN where its window is incomplete.
      real(dp), intent(in) :: x(:) !! Input time series.
      real(dp), intent(in) :: coefficients(:) !! Filter coefficients in increasing lag order.
      integer, intent(in), optional :: sides !! One for past values or two for a centered filter; defaults to two.
      real(dp), allocatable :: filtered(:)
      call r_convolution_filter(x, coefficients, filtered, sides)
   end function filter_linear

   pure function filter_recursive(x, coefficients, initial) result(filtered)
      !! Applies an autoregressive filter using R's most-recent-first initial-state convention.
      real(dp), intent(in) :: x(:) !! Input time series.
      real(dp), intent(in) :: coefficients(:) !! Recursive coefficients in increasing lag order.
      real(dp), intent(in), optional :: initial(:) !! States `y(0), y(-1), ...`; defaults to zeros.
      real(dp), allocatable :: filtered(:)

      call r_recursive_filter(x, coefficients, filtered, initial)
   end function filter_recursive

   pure function acf(x, lag_max, acf_type, demean) result(fit)
      !! Computes univariate sample autocorrelations or biased autocovariances.
      real(dp), intent(in) :: x(:) !! Finite univariate time series.
      integer, intent(in), optional :: lag_max !! Largest nonnegative lag; defaults to R's logarithmic rule.
      character(len=*), intent(in), optional :: acf_type !! `correlation` or `covariance`; defaults to correlation.
      logical, intent(in), optional :: demean !! Subtract the sample mean when true; defaults to true.
      type(acf_result_t) :: fit
      real(dp), allocatable :: values(:)
      character(len=:), allocatable :: statistic_type
      integer :: i, lower, upper

      statistic_type = "correlation"
      if (present(acf_type)) statistic_type = trim(adjustl(acf_type))
      select case (statistic_type)
      case ("correlation")
         call r_autocorrelation(x, values, lag_max, demean, fit%status)
      case ("covariance")
         fit%covariance = .true.
         call r_autocovariance(x, values, lag_max, demean, fit%status)
      case default
         fit%status = r_invalid_input
         return
      end select
      if (fit%status /= r_ok) return
      lower = lbound(values, 1)
      upper = ubound(values, 1)
      allocate (fit%lag(size(values)), fit%value(size(values)))
      do i = lower, upper
         fit%lag(i - lower + 1) = real(i, dp)
         fit%value(i - lower + 1) = values(i)
      end do
      fit%n_used = size(x)
   end function acf

   pure function pacf(x, lag_max, demean) result(fit)
      !! Computes Yule-Walker partial autocorrelations for a univariate series.
      real(dp), intent(in) :: x(:) !! Finite univariate time series.
      integer, intent(in), optional :: lag_max !! Largest lag; defaults to R's logarithmic rule.
      logical, intent(in), optional :: demean !! Subtract the sample mean when true; defaults to true.
      type(acf_result_t) :: fit
      real(dp), allocatable :: values(:)
      integer :: i

      call r_partial_autocorrelation(x, values, lag_max, demean, fit%status)
      if (fit%status /= r_ok) return
      allocate (fit%lag(size(values)), fit%value(size(values)))
      do i = 1, size(values)
         fit%lag(i) = real(i, dp)
      end do
      fit%value = values
      fit%n_used = size(x)
   end function pacf

   pure function ccf(x, y, lag_max, acf_type, demean) result(fit)
      !! Computes cross-correlations or biased cross-covariances at positive and negative lags.
      real(dp), intent(in) :: x(:) !! First finite univariate time series.
      real(dp), intent(in) :: y(:) !! Second finite series with the same size as `x`.
      integer, intent(in), optional :: lag_max !! Largest absolute lag; defaults to R's logarithmic rule.
      character(len=*), intent(in), optional :: acf_type !! `correlation` or `covariance`; defaults to correlation.
      logical, intent(in), optional :: demean !! Subtract sample means when true; defaults to true.
      type(acf_result_t) :: fit
      real(dp), allocatable :: data(:, :), values(:, :, :)
      character(len=:), allocatable :: statistic_type
      integer :: i, lag_count, maximum_lag, n

      n = size(x)
      if (n /= size(y) .or. n == 0) then
         fit%status = r_invalid_input
         return
      end if
      maximum_lag = min(n - 1, max(1, int(10.0_dp*log10(real(n, dp)))))
      if (present(lag_max)) maximum_lag = lag_max
      if (maximum_lag < 0 .or. maximum_lag >= n) then
         fit%status = r_invalid_input
         return
      end if
      allocate (data(n, 2))
      data(:, 1) = x
      data(:, 2) = y
      statistic_type = "correlation"
      if (present(acf_type)) statistic_type = trim(adjustl(acf_type))
      select case (statistic_type)
      case ("correlation")
         call r_cross_correlation(data, values, maximum_lag, demean, fit%status)
      case ("covariance")
         fit%covariance = .true.
         call r_cross_covariance(data, values, maximum_lag, demean, fit%status)
      case default
         fit%status = r_invalid_input
         return
      end select
      if (fit%status /= r_ok) return
      lag_count = 2*maximum_lag + 1
      allocate (fit%lag(lag_count), fit%value(lag_count))
      do i = -maximum_lag, maximum_lag
         fit%lag(i + maximum_lag + 1) = real(i, dp)
         if (i < 0) then
            fit%value(i + maximum_lag + 1) = values(-i, 1, 2)
         else
            fit%value(i + maximum_lag + 1) = values(i, 2, 1)
         end if
      end do
      fit%n_used = n
   end function ccf

end module r_stats_time_series
