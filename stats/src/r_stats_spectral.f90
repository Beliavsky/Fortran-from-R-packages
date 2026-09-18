! SPDX-License-Identifier: MIT
! SPDX-FileComment: Spectral-analysis functions corresponding to selected R stats functions.
module r_stats_spectral
   use r_fourier, only: r_fft, r_next_smooth_length
   use r_kinds, only: dp
   use r_missing, only: r_is_finite
   use r_optional, only: optval
   use r_spectral, only: r_circular_smooth, r_modified_daniell_kernel
   use r_spectral, only: r_spectral_kernel_bandwidth, r_spectral_kernel_df
   use r_status, only: r_invalid_input, r_ok
   use r_stats_types, only: multivariate_spectrum_result_t, spectrum_result_t
   implicit none
   private

   public :: spec_pgram, spectrum

   interface spec_pgram
      module procedure spec_pgram_matrix, spec_pgram_vector
   end interface spec_pgram

   interface spectrum
      module procedure spectrum_matrix, spectrum_vector
   end interface spectrum

contains

   pure function spec_pgram_vector(x, spans, taper, pad, fast, demean, detrend, &
                                   frequency) result(estimate)
      !! Computes a raw or modified-Daniell-smoothed periodogram using R's conventions.
      real(dp), intent(in) :: x(:) !! Finite univariate time series.
      integer, intent(in), optional :: spans(:) !! Odd smoothing spans applied successively.
      real(dp), intent(in), optional :: taper !! Fraction tapered at each end; defaults to 0.1.
      real(dp), intent(in), optional :: pad !! Fractional zero padding; defaults to zero.
      logical, intent(in), optional :: fast !! Pad to a 2/3/5-smooth transform length when true.
      logical, intent(in), optional :: demean !! Remove the sample mean; defaults to false.
      logical, intent(in), optional :: detrend !! Remove a linear trend when true; defaults to true.
      real(dp), intent(in), optional :: frequency !! Sampling frequency; defaults to one.
      type(spectrum_result_t) :: estimate
      complex(dp), allocatable :: transform(:)
      real(dp), allocatable :: periodogram(:), smoothed_periodogram(:)
      real(dp) :: bandwidth_scale, sampling_frequency
      real(dp) :: taper_fraction, taper_second_moment, taper_fourth_moment
      integer :: i, kernel_status, n_initial, n_padded, n_spectrum
      logical :: do_demean, do_detrend, do_fast

      n_initial = size(x)
      if (n_initial < 2 .or. any(.not. r_is_finite(x))) then
         estimate%status = r_invalid_input
         return
      end if
      taper_fraction = optval(taper, 0.1_dp)
      sampling_frequency = optval(frequency, 1.0_dp)
      if (taper_fraction < 0.0_dp .or. taper_fraction > 0.5_dp .or. &
          sampling_frequency <= 0.0_dp .or. optval(pad, 0.0_dp) < 0.0_dp) then
         estimate%status = r_invalid_input
         return
      end if
      do_fast = optval(fast, .true.)
      do_demean = optval(demean, .false.)
      do_detrend = optval(detrend, .true.)
      n_padded = n_initial + nint(real(n_initial, dp)*optval(pad, 0.0_dp))
      if (do_fast) n_padded = r_next_smooth_length(n_padded)
      allocate (transform(n_padded))
      call prepare_transform(x, taper_fraction, do_demean, do_detrend, transform)

      allocate (periodogram(n_padded))
      periodogram = abs(transform)**2/(real(n_initial, dp)*sampling_frequency)
      if (present(spans)) then
         call r_modified_daniell_kernel(spans, estimate%kernel, kernel_status)
         if (kernel_status /= r_ok) then
            estimate%status = kernel_status
            return
         end if
      else
         allocate (estimate%kernel(1), source=1.0_dp)
      end if
      if (size(estimate%kernel) > 1) then
         periodogram(1) = 0.5_dp*(periodogram(2) + periodogram(n_padded))
         call r_circular_smooth(periodogram, estimate%kernel, smoothed_periodogram, kernel_status)
         if (kernel_status /= r_ok) then
            estimate%status = kernel_status
            return
         end if
         call move_alloc(smoothed_periodogram, periodogram)
         estimate%smoothed = .true.
      end if

      taper_second_moment = 1.0_dp - (5.0_dp/4.0_dp)*taper_fraction
      taper_fourth_moment = 1.0_dp - (93.0_dp/64.0_dp)*taper_fraction
      n_spectrum = n_padded/2
      allocate (estimate%frequency(n_spectrum), estimate%spectrum(n_spectrum))
      do i = 1, n_spectrum
         estimate%frequency(i) = real(i, dp)*sampling_frequency/real(n_padded, dp)
         estimate%spectrum(i) = periodogram(i + 1)/taper_second_moment
      end do
      estimate%degrees_freedom = r_spectral_kernel_df(estimate%kernel)* &
         taper_second_moment**2/taper_fourth_moment*real(n_initial, dp)/real(n_padded, dp)
      bandwidth_scale = r_spectral_kernel_bandwidth(estimate%kernel)
      estimate%bandwidth = bandwidth_scale*sampling_frequency/real(n_padded, dp)
      estimate%taper = taper_fraction
      estimate%n_used = n_padded
      estimate%n_original = n_initial
      estimate%demeaned = do_demean
      estimate%detrended = do_detrend
   end function spec_pgram_vector

   pure function spectrum_vector(x, spans, taper, pad, fast, demean, detrend, &
                                 frequency) result(estimate)
      !! Provides the default periodogram path corresponding to R's `spectrum` generic.
      real(dp), intent(in) :: x(:) !! Finite univariate time series.
      integer, intent(in), optional :: spans(:) !! Odd smoothing spans applied successively.
      real(dp), intent(in), optional :: taper !! Fraction tapered at each end; defaults to 0.1.
      real(dp), intent(in), optional :: pad !! Fractional zero padding; defaults to zero.
      logical, intent(in), optional :: fast !! Pad to a 2/3/5-smooth transform length when true.
      logical, intent(in), optional :: demean !! Remove the sample mean; defaults to false.
      logical, intent(in), optional :: detrend !! Remove a linear trend when true; defaults to true.
      real(dp), intent(in), optional :: frequency !! Sampling frequency; defaults to one.
      type(spectrum_result_t) :: estimate

      estimate = spec_pgram_vector(x, spans, taper, pad, fast, demean, detrend, frequency)
   end function spectrum_vector

   pure function spec_pgram_matrix(x, spans, taper, pad, fast, demean, detrend, &
                                   frequency) result(estimate)
      !! Computes raw or modified-Daniell-smoothed multivariate periodograms.
      real(dp), intent(in) :: x(:, :) !! Finite time series arranged as `(observation, series)`.
      integer, intent(in), optional :: spans(:) !! Odd smoothing spans applied successively.
      real(dp), intent(in), optional :: taper !! Fraction tapered at each end; defaults to 0.1.
      real(dp), intent(in), optional :: pad !! Fractional zero padding; defaults to zero.
      logical, intent(in), optional :: fast !! Pad to a 2/3/5-smooth transform length when true.
      logical, intent(in), optional :: demean !! Remove each sample mean; defaults to false.
      logical, intent(in), optional :: detrend !! Remove each linear trend; defaults to true.
      real(dp), intent(in), optional :: frequency !! Sampling frequency; defaults to one.
      type(multivariate_spectrum_result_t) :: estimate
      complex(dp), allocatable :: cross_periodogram(:, :, :), smoothed_cross(:), transforms(:, :)
      real(dp) :: sampling_frequency, taper_fraction, taper_second_moment, taper_fourth_moment
      integer :: i, j, k, kernel_status, n_initial, n_padded, n_series, n_spectrum
      logical :: do_demean, do_detrend, do_fast

      n_initial = size(x, 1)
      n_series = size(x, 2)
      if (n_initial < 2 .or. n_series < 2 .or. any(.not. r_is_finite(x))) then
         estimate%status = r_invalid_input
         return
      end if
      taper_fraction = optval(taper, 0.1_dp)
      sampling_frequency = optval(frequency, 1.0_dp)
      if (taper_fraction < 0.0_dp .or. taper_fraction > 0.5_dp .or. &
          sampling_frequency <= 0.0_dp .or. optval(pad, 0.0_dp) < 0.0_dp) then
         estimate%status = r_invalid_input
         return
      end if
      do_fast = optval(fast, .true.)
      do_demean = optval(demean, .false.)
      do_detrend = optval(detrend, .true.)
      n_padded = n_initial + nint(real(n_initial, dp)*optval(pad, 0.0_dp))
      if (do_fast) n_padded = r_next_smooth_length(n_padded)
      allocate (transforms(n_padded, n_series))
      do j = 1, n_series
         call prepare_transform(x(:, j), taper_fraction, do_demean, do_detrend, transforms(:, j))
      end do

      if (present(spans)) then
         call r_modified_daniell_kernel(spans, estimate%kernel, kernel_status)
         if (kernel_status /= r_ok) then
            estimate%status = kernel_status
            return
         end if
      else
         allocate (estimate%kernel(1), source=1.0_dp)
      end if
      allocate (cross_periodogram(n_padded, n_series, n_series))
      do j = 1, n_series
         do k = 1, n_series
            cross_periodogram(:, j, k) = transforms(:, j)*conjg(transforms(:, k))/ &
                                         (real(n_initial, dp)*sampling_frequency)
            if (size(estimate%kernel) > 1) then
               cross_periodogram(1, j, k) = 0.5_dp* &
                  (cross_periodogram(2, j, k) + cross_periodogram(n_padded, j, k))
               call r_circular_smooth(cross_periodogram(:, j, k), estimate%kernel, &
                                      smoothed_cross, kernel_status)
               if (kernel_status /= r_ok) then
                  estimate%status = kernel_status
                  return
               end if
               cross_periodogram(:, j, k) = smoothed_cross
            end if
         end do
      end do

      taper_second_moment = 1.0_dp - (5.0_dp/4.0_dp)*taper_fraction
      taper_fourth_moment = 1.0_dp - (93.0_dp/64.0_dp)*taper_fraction
      n_spectrum = n_padded/2
      allocate (estimate%frequency(n_spectrum), estimate%spectrum(n_spectrum, n_series))
      allocate (estimate%cross_spectrum(n_spectrum, n_series, n_series))
      allocate (estimate%coherence(n_spectrum, n_series, n_series))
      allocate (estimate%phase(n_spectrum, n_series, n_series))
      do i = 1, n_spectrum
         estimate%frequency(i) = real(i, dp)*sampling_frequency/real(n_padded, dp)
      end do
      do j = 1, n_series
         estimate%spectrum(:, j) = real(cross_periodogram(2:n_spectrum + 1, j, j), dp)/ &
                                   taper_second_moment
      end do
      do j = 1, n_series
         do k = 1, n_series
            estimate%cross_spectrum(:, j, k) = cross_periodogram(2:n_spectrum + 1, j, k)/ &
                                                taper_second_moment
            estimate%coherence(:, j, k) = abs(cross_periodogram(2:n_spectrum + 1, j, k))**2/ &
                                           (estimate%spectrum(:, j)*estimate%spectrum(:, k)* &
                                            taper_second_moment**2)
            estimate%phase(:, j, k) = atan2(aimag(cross_periodogram(2:n_spectrum + 1, j, k)), &
                                             real(cross_periodogram(2:n_spectrum + 1, j, k), dp))
         end do
      end do
      estimate%degrees_freedom = r_spectral_kernel_df(estimate%kernel)* &
         taper_second_moment**2/taper_fourth_moment*real(n_initial, dp)/real(n_padded, dp)
      estimate%bandwidth = r_spectral_kernel_bandwidth(estimate%kernel)* &
                           sampling_frequency/real(n_padded, dp)
      estimate%taper = taper_fraction
      estimate%n_used = n_padded
      estimate%n_original = n_initial
      estimate%demeaned = do_demean
      estimate%detrended = do_detrend
      estimate%smoothed = size(estimate%kernel) > 1
   end function spec_pgram_matrix

   pure function spectrum_matrix(x, spans, taper, pad, fast, demean, detrend, &
                                 frequency) result(estimate)
      !! Provides the multivariate default periodogram path corresponding to R's `spectrum` generic.
      real(dp), intent(in) :: x(:, :) !! Finite time series arranged as `(observation, series)`.
      integer, intent(in), optional :: spans(:) !! Odd smoothing spans applied successively.
      real(dp), intent(in), optional :: taper !! Fraction tapered at each end; defaults to 0.1.
      real(dp), intent(in), optional :: pad !! Fractional zero padding; defaults to zero.
      logical, intent(in), optional :: fast !! Pad to a 2/3/5-smooth transform length when true.
      logical, intent(in), optional :: demean !! Remove each sample mean; defaults to false.
      logical, intent(in), optional :: detrend !! Remove each linear trend; defaults to true.
      real(dp), intent(in), optional :: frequency !! Sampling frequency; defaults to one.
      type(multivariate_spectrum_result_t) :: estimate

      estimate = spec_pgram_matrix(x, spans, taper, pad, fast, demean, detrend, frequency)
   end function spectrum_matrix

   pure subroutine prepare_transform(x, taper_fraction, demean, detrend, transform)
      !! Detrends, tapers, zero-pads, and Fourier-transforms one finite series.
      real(dp), intent(in) :: x(:) !! Input observations before preprocessing.
      real(dp), intent(in) :: taper_fraction !! Fraction tapered at each end.
      logical, intent(in) :: demean !! Whether to remove the sample mean.
      logical, intent(in) :: detrend !! Whether to remove a fitted linear trend and intercept.
      complex(dp), intent(out) :: transform(:) !! Transform sized to the padded length.
      real(dp), allocatable :: work(:)
      real(dp) :: slope, sum_time_squared, time_midpoint
      integer :: i, n

      n = size(x)
      allocate (work(size(transform)), source=0.0_dp)
      work(1:n) = x
      if (detrend) then
         time_midpoint = 0.5_dp*real(n + 1, dp)
         sum_time_squared = real(n, dp)*(real(n, dp)**2 - 1.0_dp)/12.0_dp
         slope = 0.0_dp
         do i = 1, n
            slope = slope + work(i)*(real(i, dp) - time_midpoint)
         end do
         slope = slope/sum_time_squared
         work(1:n) = work(1:n) - sum(work(1:n))/real(n, dp)
         do i = 1, n
            work(i) = work(i) - slope*(real(i, dp) - time_midpoint)
         end do
      else if (demean) then
         work(1:n) = work(1:n) - sum(work(1:n))/real(n, dp)
      end if
      call taper_series(work(1:n), taper_fraction)
      transform = cmplx(work, 0.0_dp, dp)
      call r_fft(transform)
   end subroutine prepare_transform

   pure subroutine taper_series(x, fraction)
      !! Applies R's split-cosine-bell taper in place.
      real(dp), intent(inout) :: x(:) !! Series values to taper.
      real(dp), intent(in) :: fraction !! Fraction tapered at each end.
      real(dp) :: weight
      integer :: i, taper_count

      taper_count = int(floor(real(size(x), dp)*fraction))
      do i = 1, taper_count
         weight = 0.5_dp*(1.0_dp - cos(acos(-1.0_dp)*real(2*i - 1, dp)/real(2*taper_count, dp)))
         x(i) = weight*x(i)
         x(size(x) - i + 1) = weight*x(size(x) - i + 1)
      end do
   end subroutine taper_series

end module r_stats_spectral
