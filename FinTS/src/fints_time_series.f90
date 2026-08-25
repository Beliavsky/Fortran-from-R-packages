! SPDX-License-Identifier: GPL-2.0-or-later
module fints_time_series
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use fints_kinds, only : dp
   use fints_status, only : fints_ok, fints_invalid_input, fints_singular
   use fints_types, only : acf_result, cross_acf_result, test_result
   use fints_summary_mod, only : sample_mean
   use fints_special, only : chi_square_survival, rank_average
   use fints_linalg, only : least_squares
   use r_status, only : r_core_ok => r_ok, r_core_singular => r_singular
   use r_time_series, only : core_autocorrelation => r_autocorrelation, &
      core_autocovariance => r_autocovariance, core_cross_correlation => r_cross_correlation, &
      core_cross_covariance => r_cross_covariance
   implicit none
   private
   public :: acf, cross_acf, autocor_test, arch_test, pacf_from_acf

contains

   subroutine acf(x, result, lag_max, acf_type, demean)
      real(dp), intent(in) :: x(:)
      type(acf_result), intent(out) :: result
      integer, intent(in), optional :: lag_max
      character(len=*), intent(in), optional :: acf_type
      logical, intent(in), optional :: demean
      character(len=16) :: kind
      logical :: center
      integer :: n, m, k, stat, core_status
      real(dp), allocatable :: correlations(:), partial_values(:), core_values(:)

      result = acf_result()
      n = size(x)
      if (n < 2 .or. any(ieee_is_nan(x))) then
         result%status = fints_invalid_input
         return
      end if

      m = min(n - 1, max(1, int(10.0_dp * log10(real(n, dp)))))
      if (present(lag_max)) m = min(n - 1, max(0, lag_max))
      kind = 'correlation'
      if (present(acf_type)) kind = lowercase(adjustl(acf_type))
      center = .true.
      if (present(demean)) center = demean
      select case (trim(kind))
      case ('correlation', 'covariance')
         if (trim(kind) == 'correlation') then
            call core_autocorrelation(x, core_values, lag_max=m, demean=center, status=core_status)
         else
            call core_autocovariance(x, core_values, lag_max=m, demean=center, status=core_status)
         end if
         if (core_status /= r_core_ok) then
            if (core_status == r_core_singular) then
               result%status = fints_singular
            else
               result%status = fints_invalid_input
            end if
            return
         end if
         allocate(result%lag(m + 1), result%value(m + 1))
         do k = 0, m
            result%lag(k + 1) = real(k, dp)
            result%value(k + 1) = core_values(k)
         end do
      case ('partial')
         allocate(correlations(m + 1))
         call core_autocorrelation(x, core_values, lag_max=m, demean=center, status=core_status)
         if (core_status /= r_core_ok) then
            if (core_status == r_core_singular) then
               result%status = fints_singular
            else
               result%status = fints_invalid_input
            end if
            return
         end if
         do k = 0, m
            correlations(k + 1) = core_values(k)
         end do
         call pacf_from_acf(correlations, partial_values, stat)
         if (stat /= fints_ok) then
            result%status = stat
            return
         end if
         allocate(result%lag(m), result%value(m))
         do k = 1, m
            result%lag(k) = real(k, dp)
            result%value(k) = partial_values(k)
         end do
      case default
         result%status = fints_invalid_input
         return
      end select

      result%n_used = n
      result%acf_type = trim(kind)
      result%status = fints_ok
   end subroutine acf

   subroutine cross_acf(x, result, lag_max, acf_type, demean)
      real(dp), intent(in) :: x(:,:)
      type(cross_acf_result), intent(out) :: result
      integer, intent(in), optional :: lag_max
      character(len=*), intent(in), optional :: acf_type
      logical, intent(in), optional :: demean
      real(dp), allocatable :: core_values(:,:,:)
      character(len=16) :: kind
      logical :: center
      integer :: n, p, m, k, core_status

      result = cross_acf_result()
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1 .or. any(ieee_is_nan(x))) then
         result%status = fints_invalid_input
         return
      end if
      m = min(n - 1, max(1, int(10.0_dp * log10(real(n, dp)))))
      if (present(lag_max)) m = min(n - 1, max(0, lag_max))
      kind = 'correlation'
      if (present(acf_type)) kind = lowercase(adjustl(acf_type))
      if (trim(kind) /= 'correlation' .and. trim(kind) /= 'covariance') then
         result%status = fints_invalid_input
         return
      end if
      center = .true.
      if (present(demean)) center = demean

      if (trim(kind) == 'correlation') then
         call core_cross_correlation(x, core_values, lag_max=m, demean=center, status=core_status)
      else
         call core_cross_covariance(x, core_values, lag_max=m, demean=center, status=core_status)
      end if
      if (core_status /= r_core_ok) then
         if (core_status == r_core_singular) then
            result%status = fints_singular
         else
            result%status = fints_invalid_input
         end if
         return
      end if

      allocate(result%lag(m + 1), result%value(m + 1, p, p))
      do k = 0, m
         result%lag(k + 1) = real(k, dp)
         result%value(k + 1, :, :) = core_values(k, :, :)
      end do
      result%n_used = n
      result%acf_type = trim(kind)
      result%status = fints_ok
   end subroutine cross_acf

   subroutine pacf_from_acf(correlation, partial, status)
      real(dp), intent(in) :: correlation(:)
      real(dp), allocatable, intent(out) :: partial(:)
      integer, intent(out) :: status
      real(dp), allocatable :: phi(:,:), variance(:)
      real(dp) :: reflection
      integer :: m, k, j

      m = size(correlation) - 1
      allocate(partial(max(0, m)))
      if (m < 1 .or. abs(correlation(1)) <= tiny(1.0_dp)) then
         status = fints_invalid_input
         return
      end if
      allocate(phi(m, m), variance(0:m))
      phi = 0.0_dp
      variance(0) = correlation(1)
      do k = 1, m
         reflection = correlation(k + 1)
         if (k > 1) reflection = reflection - &
            dot_product(phi(k - 1, 1:k - 1), correlation(k:2:-1))
         if (abs(variance(k - 1)) <= tiny(1.0_dp)) then
            status = fints_singular
            return
         end if
         reflection = reflection / variance(k - 1)
         phi(k, k) = reflection
         do j = 1, k - 1
            phi(k, j) = phi(k - 1, j) - reflection * phi(k - 1, k - j)
         end do
         variance(k) = variance(k - 1) * (1.0_dp - reflection * reflection)
         partial(k) = reflection
      end do
      status = fints_ok
   end subroutine pacf_from_acf

   subroutine autocor_test(x, result, lag, test_type, degrees_freedom)
      real(dp), intent(in) :: x(:)
      type(test_result), intent(out) :: result
      integer, intent(in), optional :: lag
      character(len=*), intent(in), optional :: test_type
      real(dp), intent(in), optional :: degrees_freedom
      real(dp), allocatable :: data(:)
      type(acf_result) :: acf_values
      character(len=16) :: kind
      real(dp) :: statistic, df
      integer :: n, m, k

      result = test_result()
      n = size(x)
      if (n < 3 .or. any(ieee_is_nan(x))) then
         result%status = fints_invalid_input
         return
      end if
      m = ceiling(log(real(n, dp)))
      if (present(lag)) m = lag
      if (m < 1 .or. m >= n) then
         result%status = fints_invalid_input
         return
      end if
      kind = 'ljung-box'
      if (present(test_type)) kind = lowercase(adjustl(test_type))
      allocate(data(n))
      data = x
      if (trim(kind) == 'rank') then
         call rank_average(x, data)
         kind = 'ljung-box'
      end if
      if (trim(kind) /= 'ljung-box' .and. trim(kind) /= 'box-pierce') then
         result%status = fints_invalid_input
         return
      end if

      call acf(data, acf_values, lag_max=m, acf_type='correlation')
      if (acf_values%status /= fints_ok) then
         result%status = acf_values%status
         return
      end if
      statistic = 0.0_dp
      if (trim(kind) == 'box-pierce') then
         statistic = real(n, dp) * sum(acf_values%value(2:m + 1) ** 2)
         result%method = 'Box-Pierce test'
      else
         do k = 1, m
            statistic = statistic + acf_values%value(k + 1) ** 2 / real(n - k, dp)
         end do
         statistic = real(n * (n + 2), dp) * statistic
         result%method = 'Ljung-Box test'
      end if
      df = real(m, dp)
      if (present(degrees_freedom)) df = max(1.0_dp, degrees_freedom)
      result%statistic = statistic
      result%degrees_freedom = df
      result%p_value = chi_square_survival(statistic, df)
      result%lag = m
      result%n_observations = n
      result%status = fints_ok
   end subroutine autocor_test

   subroutine arch_test(x, result, lags, demean)
      real(dp), intent(in) :: x(:)
      type(test_result), intent(out) :: result
      integer, intent(in), optional :: lags
      logical, intent(in), optional :: demean
      real(dp), allocatable :: data(:), design(:,:), response(:), beta(:), residuals(:)
      real(dp) :: center, sse, sst, r_squared
      logical :: remove_mean
      integer :: n, m, rows, i, j, stat

      result = test_result()
      n = size(x)
      m = 12
      if (present(lags)) m = lags
      if (n <= m + 1 .or. m < 1 .or. any(ieee_is_nan(x))) then
         result%status = fints_invalid_input
         return
      end if
      remove_mean = .false.
      if (present(demean)) remove_mean = demean
      center = 0.0_dp
      if (remove_mean) center = sample_mean(x)
      allocate(data(n))
      data = (x - center) ** 2
      rows = n - m
      allocate(design(rows, m + 1), response(rows))
      design(:, 1) = 1.0_dp
      do i = 1, rows
         response(i) = data(m + i)
         do j = 1, m
            design(i, j + 1) = data(m + i - j)
         end do
      end do
      call least_squares(design, response, beta, residuals, sse, stat)
      if (stat /= fints_ok) then
         result%status = stat
         return
      end if
      sst = sum((response - sample_mean(response)) ** 2)
      r_squared = 0.0_dp
      if (sst > 0.0_dp) r_squared = max(0.0_dp, min(1.0_dp, 1.0_dp - sse / sst))
      result%statistic = r_squared * real(rows, dp)
      result%degrees_freedom = real(m, dp)
      result%p_value = chi_square_survival(result%statistic, result%degrees_freedom)
      result%lag = m
      result%n_observations = rows
      result%method = 'ARCH LM test'
      result%status = fints_ok
   end subroutine arch_test

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lowercase

end module fints_time_series
