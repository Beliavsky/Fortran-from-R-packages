! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_covariance
   use r_kinds, only : dp
   implicit none
   private

   public :: ftsa_facf, ftsa_long_run_covariance

contains

   pure subroutine ftsa_facf(fun_data, lags, rho, status)
      real(dp), intent(in) :: fun_data(:, :) !! Time-by-grid functional observations, matching upstream `facf` orientation.
      integer, intent(in) :: lags(:) !! Integer lags; absolute values must be smaller than the number of time points.
      real(dp), allocatable, intent(out) :: rho(:) !! Functional autocorrelations; lag zero is returned as zero.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid lag or empty data.
      real(dp), allocatable :: centered(:, :), gamma0(:, :), gammal(:, :), means(:)
      real(dp) :: denom
      integer :: h, i, n, p, t

      n = size(fun_data, 1)
      p = size(fun_data, 2)
      if (n < 2 .or. p < 1 .or. any(abs(lags) >= n)) then
         allocate(rho(0))
         status = 1
         return
      end if
      allocate(centered(n, p), means(p), gamma0(p, p), gammal(p, p), rho(size(lags)))
      means = sum(fun_data, dim=1) / real(n, dp)
      do t = 1, n
         centered(t, :) = fun_data(t, :) - means
      end do
      gamma0 = matmul(transpose(centered), centered) / real(n, dp)
      denom = sum([(gamma0(i, i), i=1, p)])
      if (denom <= tiny(1.0_dp)) then
         rho = 0.0_dp
         status = 1
         return
      end if
      do i = 1, size(lags)
         h = lags(i)
         if (h == 0) then
            rho(i) = 0.0_dp
         else
            call lag_covariance_rows(centered, h, gammal)
            rho(i) = sqrt(sum(gammal * gammal)) / denom
         end if
      end do
      status = 0
   end subroutine ftsa_facf

   pure subroutine ftsa_long_run_covariance(dat, covariance, bandwidth, c0, h_window, status)
      real(dp), intent(in) :: dat(:, :) !! Grid-by-time functional observations.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Bartlett long-run covariance estimate.
      real(dp), intent(out) :: bandwidth !! Plug-in Bartlett bandwidth estimated from a flat-top pilot.
      real(dp), intent(in), optional :: c0 !! Adaptive pilot threshold multiplier; upstream default is 3.
      integer, intent(in), optional :: h_window !! Consecutive small-correlation window; upstream default is 3.
      integer, intent(out) :: status !! Zero on success; nonzero for insufficient dimensions or degeneracy.
      real(dp), allocatable :: centered(:, :), c0mat(:, :), c2mat(:, :), gamma(:, :), rho(:)
      real(dp) :: c0_local, const, denom, first_part, kernel_square_int, second_part, w, pilot_band
      integer :: h_local, i, lag, n, p, selected

      p = size(dat, 1)
      n = size(dat, 2)
      c0_local = 3.0_dp
      if (present(c0)) c0_local = c0
      h_local = 3
      if (present(h_window)) h_local = h_window
      if (p < 1 .or. n < max(4, h_local + 1) .or. c0_local <= 0.0_dp .or. h_local < 1) then
         allocate(covariance(0, 0))
         bandwidth = 0.0_dp
         status = 1
         return
      end if
      allocate(centered(p, n), rho(n - 1), gamma(p, p), c0mat(p, p), c2mat(p, p))
      do i = 1, p
         centered(i, :) = dat(i, :) - sum(dat(i, :)) / real(n, dp)
      end do
      call lag_covariance_columns(centered, 0, gamma)
      denom = sum([(gamma(i, i), i=1, p)])**2
      if (denom <= tiny(1.0_dp)) then
         allocate(covariance(p, p))
         covariance = 0.0_dp
         bandwidth = 0.0_dp
         status = 1
         return
      end if
      do lag = 1, n - 1
         call lag_covariance_columns(centered, lag, gamma)
         rho(lag) = sqrt(sum(gamma * gamma) / denom)
      end do
      const = c0_local * sqrt(log10(real(n, dp)) / real(n, dp))
      selected = adaptive_cutoff(rho, const, h_local)
      ! Upstream computes the adaptive cutoff but its returned value is not used in the final plug-in call.
      pilot_band = real(n, dp)**0.2_dp + 0.0_dp * real(max(selected, 0), dp)
      c0mat = 0.0_dp
      c2mat = 0.0_dp
      call lag_covariance_columns(centered, 0, gamma)
      c0mat = gamma
      do lag = 1, n - 1
         w = flat_top(real(lag, dp) / pilot_band)
         if (w <= 0.0_dp) cycle
         call lag_covariance_columns(centered, lag, gamma)
         c0mat = c0mat + w * (gamma + transpose(gamma))
         c2mat = c2mat + w * real(lag, dp) * (gamma + transpose(gamma))
      end do
      first_part = (2.0_dp * sum(c2mat * c2mat))**(1.0_dp / 3.0_dp)
      kernel_square_int = 2.0_dp / 3.0_dp
      second_part = ((sum(c0mat * c0mat) + sum([(c0mat(i, i), i=1, p)])**2) * &
         kernel_square_int)**(-1.0_dp / 3.0_dp)
      bandwidth = first_part * second_part * real(n, dp)**(1.0_dp / 3.0_dp)
      if (.not. (bandwidth > 0.0_dp)) bandwidth = max(1.0_dp, pilot_band)
      allocate(covariance(p, p))
      call lag_covariance_columns(centered, 0, gamma)
      covariance = gamma
      do lag = 1, n - 1
         w = bartlett(real(lag, dp) / bandwidth)
         if (w <= 0.0_dp) cycle
         call lag_covariance_columns(centered, lag, gamma)
         covariance = covariance + w * (gamma + transpose(gamma))
      end do
      status = 0
   end subroutine ftsa_long_run_covariance

   pure subroutine lag_covariance_rows(centered, lag, covariance)
      real(dp), intent(in) :: centered(:, :) !! Time-by-grid centered observations.
      integer, intent(in) :: lag !! Signed lag with magnitude less than the number of rows.
      real(dp), intent(out) :: covariance(:, :) !! Lag covariance divided by the full time length.
      integer :: h, n, p, t

      n = size(centered, 1)
      p = size(centered, 2)
      h = abs(lag)
      covariance = 0.0_dp
      if (lag >= 0) then
         do t = 1, n - h
            covariance = covariance + spread(centered(t, :), 2, p) * spread(centered(t + h, :), 1, p)
         end do
      else
         do t = 1, n - h
            covariance = covariance + spread(centered(t + h, :), 2, p) * spread(centered(t, :), 1, p)
         end do
      end if
      covariance = covariance / real(n, dp)
   end subroutine lag_covariance_rows

   pure subroutine lag_covariance_columns(centered, lag, covariance)
      real(dp), intent(in) :: centered(:, :) !! Grid-by-time centered observations.
      integer, intent(in) :: lag !! Nonnegative lag smaller than the number of columns.
      real(dp), intent(out) :: covariance(:, :) !! Lag covariance divided by the full time length.
      integer :: n, p, t

      p = size(centered, 1)
      n = size(centered, 2)
      covariance = 0.0_dp
      do t = 1, n - lag
         covariance = covariance + spread(centered(:, t), 2, p) * spread(centered(:, t + lag), 1, p)
      end do
      covariance = covariance / real(n, dp)
   end subroutine lag_covariance_columns

   pure integer function adaptive_cutoff(rho, threshold, width) result(cutoff)
      real(dp), intent(in) :: rho(:) !! Positive-lag functional autocorrelations.
      real(dp), intent(in) :: threshold !! Small-correlation threshold.
      integer, intent(in) :: width !! Number of consecutive values required below threshold.
      integer :: i, last

      cutoff = -1
      last = size(rho) - width + 1
      do i = 1, last
         if (all(rho(i:i + width - 1) <= threshold)) then
            cutoff = i - 1
            return
         end if
      end do
   end function adaptive_cutoff

   pure real(dp) function flat_top(x) result(value)
      real(dp), intent(in) :: x !! Nonnegative normalized lag.

      if (x < 0.5_dp) then
         value = 1.0_dp
      else if (x <= 1.0_dp) then
         value = 2.0_dp - 2.0_dp * x
      else
         value = 0.0_dp
      end if
   end function flat_top

   pure real(dp) function bartlett(x) result(value)
      real(dp), intent(in) :: x !! Nonnegative normalized lag.

      value = max(0.0_dp, 1.0_dp - abs(x))
   end function bartlett

end module ftsa_covariance
