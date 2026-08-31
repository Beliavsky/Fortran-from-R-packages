! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! MCMC convergence and autocorrelation diagnostics translated from mitml.
module mitml_convergence
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_linalg, only : solve_system
   use mitml_numeric, only : mean_real, sample_variance
   use mitml_types, only : MITML_ERR_ARGUMENT, MITML_ERR_LINALG, MITML_OK
   implicit none
   private

   public :: gelman_rubin
   public :: sd_proportion
   public :: reduced_acf
   public :: moving_average

contains

   subroutine gelman_rubin(x, n_chains, rhat, status)
      real(dp), intent(in) :: x(:, :) !! MCMC traces, shape n_parameter by concatenated iteration.
      integer, intent(in) :: n_chains !! Number of equal-length chains concatenated across the second dimension.
      real(dp), intent(out) :: rhat(:) !! Gelman-Rubin potential scale reduction factors, one per parameter.
      integer, intent(out) :: status !! MITML_OK on success or an argument status code.
      real(dp), allocatable :: means(:)
      real(dp), allocatable :: variances(:)
      real(dp) :: b_div_n
      real(dp) :: covariance_v_m
      real(dp) :: covariance_v_m2
      real(dp) :: df
      real(dp) :: muhat
      real(dp) :: sighat2
      real(dp) :: var_vhat
      real(dp) :: vhat
      real(dp) :: w
      integer :: chain
      integer :: first
      integer :: i
      integer :: last
      integer :: n
      integer :: used

      status = MITML_OK
      if (size(rhat) /= size(x, 1) .or. n_chains < 2) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      n = size(x, 2) / n_chains
      if (n < 2) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      used = n * n_chains
      allocate(means(n_chains), variances(n_chains))

      do i = 1, size(x, 1)
         if (all(ieee_is_nan(x(i, 1:used)))) then
            rhat(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            cycle
         end if
         if (any(ieee_is_nan(x(i, 1:used)))) then
            rhat(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            cycle
         end if
         do chain = 1, n_chains
            first = (chain - 1) * n + 1
            last = chain * n
            means(chain) = mean_real(x(i, first:last))
            variances(chain) = sample_variance(x(i, first:last))
         end do
         b_div_n = sample_variance(means)
         w = mean_real(variances)
         muhat = mean_real(x(i, 1:used))
         sighat2 = real(n - 1, dp) / real(n, dp) * w + b_div_n
         vhat = sighat2 + b_div_n / real(n_chains, dp)

         if (b_div_n <= 0.0_dp .and. all(variances <= 0.0_dp)) then
            rhat(i) = 1.0_dp
            cycle
         end if
         if (w <= 0.0_dp) then
            rhat(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            cycle
         end if

         covariance_v_m2 = sample_covariance_pair(variances, means * means)
         covariance_v_m = sample_covariance_pair(variances, means)
         var_vhat = (real(n - 1, dp) / real(n, dp))**2 * sample_variance(variances) / real(n_chains, dp)
         var_vhat = var_vhat + (real(n_chains + 1, dp) / real(n_chains * n, dp))**2 * &
            2.0_dp / real(n_chains - 1, dp) * (b_div_n * real(n, dp))**2
         var_vhat = var_vhat + 2.0_dp * real((n_chains + 1) * (n - 1), dp) / &
            real(n_chains * n * n, dp) * real(n, dp) / real(n_chains, dp) * &
            (covariance_v_m2 - 2.0_dp * muhat * covariance_v_m)
         if (var_vhat <= 0.0_dp) then
            rhat(i) = sqrt(vhat / w)
         else
            df = 2.0_dp * vhat * vhat / var_vhat
            if (df <= 2.0_dp) then
               rhat(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
               rhat(i) = sqrt((vhat / w) * df / (df - 2.0_dp))
            end if
         end if
      end do
   end subroutine gelman_rubin

   subroutine sd_proportion(x, sdp, n_eff, status, max_order)
      real(dp), intent(in) :: x(:, :) !! MCMC traces, shape n_parameter by n_iteration.
      real(dp), intent(out) :: sdp(:) !! Square root of variance-of-mean divided by marginal chain variance.
      real(dp), intent(out) :: n_eff(:) !! Approximate effective sample sizes implied by the spectral density at zero.
      integer, intent(out) :: status !! MITML_OK on success or an argument/linear-algebra status code.
      integer, intent(in), optional :: max_order !! Maximum Yule-Walker AR order; default follows R stats::ar scale.
      real(dp) :: marginal_variance
      real(dp) :: spectral_zero
      integer :: i
      integer :: n
      integer :: order_limit
      integer :: local_status

      status = MITML_OK
      n = size(x, 2)
      if (size(sdp) /= size(x, 1) .or. size(n_eff) /= size(x, 1) .or. n < 3) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      order_limit = min(n - 1, max(1, int(10.0_dp * log10(real(n, dp)))))
      if (present(max_order)) then
         if (max_order < 0) then
            status = MITML_ERR_ARGUMENT
            return
         end if
         order_limit = min(n - 1, max_order)
      end if

      do i = 1, size(x, 1)
         if (all(ieee_is_nan(x(i, :))) .or. any(ieee_is_nan(x(i, :)))) then
            sdp(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            n_eff(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            cycle
         end if
         marginal_variance = sample_variance(x(i, :))
         if (marginal_variance <= 0.0_dp) then
            sdp(i) = 0.0_dp
            n_eff(i) = real(n, dp)
            cycle
         end if
         call ar_spectral_zero(x(i, :), order_limit, spectral_zero, local_status)
         if (local_status /= MITML_OK .or. spectral_zero <= 0.0_dp) then
            status = MITML_ERR_LINALG
            return
         end if
         sdp(i) = sqrt((spectral_zero / real(n, dp)) / marginal_variance)
         n_eff(i) = marginal_variance * real(n, dp) / spectral_zero
      end do
   end subroutine sd_proportion

   subroutine reduced_acf(x, lag, value, status, smooth, smoothing_sd)
      real(dp), intent(in) :: x(:) !! Scalar time series whose lag autocorrelation is summarized.
      integer, intent(in) :: lag !! Central nonnegative lag to evaluate.
      real(dp), intent(out) :: value !! Autocorrelation or normal-kernel-smoothed autocorrelation around lag.
      integer, intent(out) :: status !! MITML_OK on success or an argument status code.
      integer, intent(in), optional :: smooth !! Number of neighboring lags on each side; default zero.
      real(dp), intent(in), optional :: smoothing_sd !! Normal-kernel standard deviation in lag units; default 0.5.
      real(dp) :: ac
      real(dp) :: centered(size(x))
      real(dp) :: denominator
      real(dp) :: sd0
      real(dp) :: weight
      real(dp) :: weight_sum
      integer :: j
      integer :: ll
      integer :: radius

      status = MITML_OK
      value = 0.0_dp
      radius = 0
      if (present(smooth)) radius = smooth
      sd0 = 0.5_dp
      if (present(smoothing_sd)) sd0 = smoothing_sd
      if (size(x) < 2 .or. lag < 0 .or. radius < 0 .or. sd0 <= 0.0_dp) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      if (lag - radius < 0 .or. lag + radius >= size(x)) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      if (any(ieee_is_nan(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      centered = x - mean_real(x)
      denominator = sum(centered * centered)
      if (denominator <= 0.0_dp) then
         value = 0.0_dp
         return
      end if

      weight_sum = 0.0_dp
      do j = -radius, radius
         ll = lag + j
         if (ll == 0) then
            ac = 1.0_dp
         else
            ac = sum(centered(1:size(x) - ll) * centered(1 + ll:size(x))) / denominator
         end if
         weight = exp(-0.5_dp * (real(j, dp) / sd0)**2)
         value = value + weight * ac
         weight_sum = weight_sum + weight
      end do
      value = value / weight_sum
   end subroutine reduced_acf

   subroutine moving_average(x, half_window, y, status, fill_edges)
      real(dp), intent(in) :: x(:) !! Time series to smooth with a centered odd-width moving window.
      integer, intent(in) :: half_window !! Half-width B, giving an interior window width of 2*B+1.
      real(dp), intent(out) :: y(:) !! Smoothed series, the same length as x.
      integer, intent(out) :: status !! MITML_OK on success or an argument status code.
      logical, intent(in), optional :: fill_edges !! Fill edges with growing odd windows; default true as in mitml.
      logical :: fill
      integer :: i
      integer :: n
      integer :: radius

      status = MITML_OK
      n = size(x)
      if (size(y) /= n .or. half_window < 0 .or. 2 * half_window + 1 > n) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      fill = .true.
      if (present(fill_edges)) fill = fill_edges
      y = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = half_window + 1, n - half_window
         y(i) = mean_real(x(i - half_window:i + half_window))
      end do
      if (.not. fill) return
      do i = 1, half_window
         radius = i - 1
         y(i) = mean_real(x(1:2 * radius + 1))
         y(n - i + 1) = mean_real(x(n - 2 * radius:n))
      end do
   end subroutine moving_average

   pure real(dp) function sample_covariance_pair(x, y) result(value)
      real(dp), intent(in) :: x(:) !! First vector in a sample covariance calculation.
      real(dp), intent(in) :: y(:) !! Second vector in a sample covariance calculation of the same length.
      real(dp) :: mean_x
      real(dp) :: mean_y

      mean_x = sum(x) / real(size(x), dp)
      mean_y = sum(y) / real(size(y), dp)
      value = sum((x - mean_x) * (y - mean_y)) / real(size(x) - 1, dp)
   end function sample_covariance_pair

   subroutine ar_spectral_zero(x, max_order, spectral_zero, status)
      real(dp), intent(in) :: x(:) !! Complete scalar trace fitted by candidate Yule-Walker autoregressions.
      integer, intent(in) :: max_order !! Largest AR order considered by the AIC search.
      real(dp), intent(out) :: spectral_zero !! Estimated spectral density at frequency zero without the 2*pi factor.
      integer, intent(out) :: status !! MITML_OK or MITML_ERR_LINALG when every positive-order solve fails.
      real(dp), allocatable :: acov(:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: phi(:)
      real(dp), allocatable :: toeplitz(:, :)
      real(dp) :: aic
      real(dp) :: best_aic
      real(dp) :: best_innovation
      real(dp) :: best_sum_phi
      real(dp) :: centered(size(x))
      real(dp) :: innovation
      real(dp) :: sum_phi
      integer :: i
      integer :: info
      integer :: j
      integer :: n
      integer :: order

      status = MITML_OK
      n = size(x)
      centered = x - mean_real(x)
      allocate(acov(0:max_order))
      do i = 0, max_order
         acov(i) = sum(centered(1:n - i) * centered(1 + i:n)) / real(n, dp)
      end do
      if (acov(0) <= 0.0_dp) then
         spectral_zero = 0.0_dp
         return
      end if
      best_innovation = acov(0)
      best_sum_phi = 0.0_dp
      best_aic = real(n, dp) * log(best_innovation) + 2.0_dp

      do order = 1, max_order
         allocate(toeplitz(order, order), rhs(order), phi(order))
         do i = 1, order
            rhs(i) = acov(i)
            do j = 1, order
               toeplitz(i, j) = acov(abs(i - j))
            end do
         end do
         call solve_system(toeplitz, rhs, phi, info)
         if (info == 0) then
            innovation = acov(0) - dot_product(phi, rhs)
            if (innovation > tiny(1.0_dp)) then
               aic = real(n, dp) * log(innovation) + 2.0_dp * real(order + 1, dp)
               if (aic < best_aic) then
                  best_aic = aic
                  best_innovation = innovation
                  best_sum_phi = sum(phi)
               end if
            end if
         end if
         deallocate(toeplitz, rhs, phi)
      end do
      sum_phi = 1.0_dp - best_sum_phi
      if (abs(sum_phi) <= sqrt(epsilon(1.0_dp))) then
         status = MITML_ERR_LINALG
         spectral_zero = 0.0_dp
      else
         spectral_zero = best_innovation / (sum_phi * sum_phi)
      end if
   end subroutine ar_spectral_zero

end module mitml_convergence
