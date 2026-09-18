! SPDX-License-Identifier: MIT
! Direct-convolution translation of kde1d-cpp bandwidth and KDE utilities.
! Upstream copyright: Thomas Nagler and Thibault Vatter.
! The direct convolution evaluates the same binned Gaussian kernels without Eigen FFT.
module kde1d_bandwidth
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp, r_pi
   use r_quantiles, only : r_median, r_qrule_hf7, r_quantile_type7
   use r_quantiles, only : r_weighted_quantile_survey
   implicit none
   private

   public :: gaussian_derivative
   public :: kde_drv_binned
   public :: linbin
   public :: plugin_bandwidth
   public :: weighted_quantile_cpp

contains

   pure subroutine linbin(x, lower, upper, num_bins, weights, counts)
      real(dp), intent(in) :: x(:) !! Observations to linearly bin on [lower,upper].
      real(dp), intent(in) :: lower !! Lower edge of the regular binning grid.
      real(dp), intent(in) :: upper !! Upper edge of the regular binning grid.
      integer, intent(in) :: num_bins !! Number of equal-width cells; output has num_bins+1 knots.
      real(dp), intent(in) :: weights(:) !! Observation weights, conforming with x.
      real(dp), intent(out) :: counts(num_bins + 1) !! Linearly interpolated weighted knot counts.

      real(dp) :: delta
      real(dp) :: lxi
      real(dp) :: rem
      integer :: i
      integer :: li

      counts = 0.0_dp
      if (num_bins < 1 .or. upper <= lower) return
      delta = (upper - lower) / real(num_bins, dp)
      do i = 1, size(x)
         lxi = (x(i) - lower) / delta
         li = int(floor(lxi))
         rem = lxi - real(li, dp)
         if (li >= 0 .and. li < num_bins) then
            counts(li + 1) = counts(li + 1) + (1.0_dp - rem) * weights(i)
            counts(li + 2) = counts(li + 2) + rem * weights(i)
         else if (x(i) == upper) then
            counts(num_bins + 1) = counts(num_bins + 1) + weights(i)
         end if
      end do
   end subroutine linbin

   pure real(dp) function gaussian_derivative(x, drv) result(value)
      real(dp), intent(in) :: x !! Standard-normal evaluation point.
      integer, intent(in) :: drv !! Nonnegative derivative order.

      real(dp) :: h0
      real(dp) :: h1
      real(dp) :: h2
      integer :: j

      h0 = 1.0_dp
      if (drv == 0) then
         value = exp(-0.5_dp * x * x) / sqrt(2.0_dp * r_pi)
         return
      end if
      h1 = x
      if (drv == 1) then
         value = -h1 * exp(-0.5_dp * x * x) / sqrt(2.0_dp * r_pi)
         return
      end if
      do j = 2, drv
         h2 = x * h1 - real(j - 1, dp) * h0
         h0 = h1
         h1 = h2
      end do
      value = h1 * exp(-0.5_dp * x * x) / sqrt(2.0_dp * r_pi)
      if (mod(drv, 2) == 1) value = -value
   end function gaussian_derivative

   pure subroutine kde_drv_binned(x, bandwidth, lower, upper, num_bins, weights, drv, values, bin_counts)
      real(dp), intent(in) :: x(:) !! Observations used to construct the binned KDE.
      real(dp), intent(in) :: bandwidth !! Positive Gaussian-kernel bandwidth.
      real(dp), intent(in) :: lower !! Lower knot of the regular evaluation grid.
      real(dp), intent(in) :: upper !! Upper knot of the regular evaluation grid.
      integer, intent(in) :: num_bins !! Number of grid cells; values has num_bins+1 entries.
      real(dp), intent(in) :: weights(:) !! Observation weights; size zero requests equal weights.
      integer, intent(in) :: drv !! Gaussian KDE derivative order.
      real(dp), intent(out) :: values(num_bins + 1) !! Binned derivative estimate on regular grid knots.
      real(dp), intent(out), optional :: bin_counts(num_bins + 1) !! Normalized weighted linear-bin counts.

      real(dp) :: counts(num_bins + 1)
      real(dp) :: delta
      real(dp) :: kernel
      real(dp) :: norm
      real(dp) :: tau
      real(dp) :: w(size(x))
      real(dp) :: z
      integer :: j
      integer :: k

      values = 0.0_dp
      if (present(bin_counts)) bin_counts = 0.0_dp
      if (size(x) == 0 .or. num_bins < 1 .or. bandwidth <= 0.0_dp .or. upper <= lower) return

      if (size(weights) == size(x)) then
         if (sum(weights) <= 0.0_dp) return
         w = weights / (sum(weights) / real(size(weights), dp))
      else
         w = 1.0_dp
      end if
      call linbin(x, lower, upper, num_bins, w, counts)
      if (present(bin_counts)) bin_counts = counts
      norm = sum(counts)
      if (norm <= 0.0_dp) return

      delta = (upper - lower) / real(num_bins, dp)
      tau = 4.0_dp + real(drv, dp)
      do j = 0, num_bins
         do k = 0, num_bins
            z = real(j - k, dp) * delta / bandwidth
            if (abs(z) <= tau) then
               kernel = gaussian_derivative(z, drv) / bandwidth**(drv + 1)
               values(j + 1) = values(j + 1) + counts(k + 1) * kernel
            end if
         end do
      end do
      values = values / norm
   end subroutine kde_drv_binned

   pure real(dp) function weighted_quantile_cpp(x, probability, weights) result(value)
      real(dp), intent(in) :: x(:) !! Data values for the kde1d-cpp weighted quantile convention.
      real(dp), intent(in) :: probability !! Quantile probability in [0,1].
      real(dp), intent(in) :: weights(:) !! Nonnegative weights conforming with x.

      value = r_weighted_quantile_survey(x, weights, probability, r_qrule_hf7)
   end function weighted_quantile_cpp

   pure real(dp) function plugin_bandwidth(x, degree, weights) result(bandwidth)
      real(dp), intent(in) :: x(:) !! Transformed observations used for plug-in bandwidth selection.
      integer, intent(in) :: degree !! Local-polynomial degree, restricted to 0, 1, or 2.
      real(dp), intent(in) :: weights(:) !! Optional normalized case weights; size zero means equal weights.

      real(dp) :: fallback
      real(dp) :: ibias2
      real(dp) :: ivar
      real(dp) :: n_eff
      real(dp) :: scale
      real(dp) :: w(size(x))
      integer :: bandwidth_power
      logical :: ok

      if (size(x) < 2 .or. degree < 0 .or. degree > 2) then
         bandwidth = 1.0_dp
         return
      end if
      if (size(weights) == size(x)) then
         w = weights * real(size(x), dp) / sum(weights)
      else
         w = 1.0_dp
      end if
      n_eff = sum(w)**2 / sum(w * w)
      scale = scale_estimate(x, w)
      bandwidth_power = merge(4, 8, degree < 2)
      fallback = 4.0_dp * 1.06_dp * scale * n_eff**(-1.0_dp / real(bandwidth_power + 1, dp))

      if (maxval(x) <= minval(x) .or. .not. (scale > 0.0_dp)) then
         bandwidth = fallback
         return
      end if

      call integrated_bias2(x, w, degree, scale, ibias2, ok)
      ivar = merge(1.0_dp, 27.0_dp / 16.0_dp, degree < 2) * 0.5_dp / sqrt(r_pi)
      if (.not. ok .or. .not. (ibias2 > 0.0_dp) .or. .not. ieee_is_finite(ibias2)) then
         bandwidth = fallback
      else
         bandwidth = (ivar / (real(bandwidth_power, dp) * n_eff * ibias2)) &
            ** (1.0_dp / real(bandwidth_power + 1, dp))
         if (.not. ieee_is_finite(bandwidth) .or. .not. (bandwidth > 0.0_dp)) bandwidth = fallback
      end if
   end function plugin_bandwidth

   pure real(dp) function scale_estimate(x, weights) result(scale)
      real(dp), intent(in) :: x(:) !! Observations whose robust spread is estimated.
      real(dp), intent(in) :: weights(:) !! Weights normalized to sum to the observation count.

      real(dp) :: mean_x
      real(dp) :: q25
      real(dp) :: q75
      real(dp) :: sd_x

      mean_x = sum(x * weights) / real(size(x), dp)
      if (size(x) > 1) then
         sd_x = sqrt(sum((x - mean_x)**2 * weights) / real(size(x) - 1, dp))
      else
         sd_x = 0.0_dp
      end if
      if (maxval(weights) == minval(weights)) then
         q25 = r_quantile_type7(x, 0.25_dp)
         q75 = r_quantile_type7(x, 0.75_dp)
      else
         q25 = weighted_quantile_cpp(x, 0.25_dp, weights)
         q75 = weighted_quantile_cpp(x, 0.75_dp, weights)
      end if
      scale = min((q75 - q25) / 1.349_dp, sd_x)
      if (scale == 0.0_dp) scale = merge(sd_x, 1.0_dp, sd_x > 0.0_dp)
   end function scale_estimate

   pure subroutine integrated_bias2(x, weights, degree, scale, ibias2, ok)
      real(dp), intent(in) :: x(:) !! Observations used by the plug-in functional estimator.
      real(dp), intent(in) :: weights(:) !! Weights normalized to mean one.
      integer, intent(in) :: degree !! Local-polynomial degree, 0 through 2.
      real(dp), intent(in) :: scale !! Robust scale estimate used by normal-reference pilot bandwidths.
      real(dp), intent(out) :: ibias2 !! Estimated integrated squared-bias coefficient.
      logical, intent(out) :: ok !! True when all pilot calculations remained finite and positive.

      integer, parameter :: nb = 400
      real(dp) :: arg(nb + 1)
      real(dp) :: bin_counts(nb + 1)
      real(dp) :: f0(nb + 1)
      real(dp) :: f1(nb + 1)
      real(dp) :: f2(nb + 1)
      real(dp) :: f4(nb + 1)
      real(dp) :: h
      real(dp) :: tiny_f
      integer :: i

      ok = .false.
      ibias2 = 0.0_dp
      if (maxval(x) <= minval(x)) return

      if (degree == 0 .or. degree == 1) then
         call bkfe_bandwidth(x, weights, 4, scale, h, ok)
      else
         call bkfe_bandwidth(x, weights, 8, scale, h, ok)
      end if
      if (.not. ok) return

      call kde_drv_binned(x, h, minval(x), maxval(x), nb, weights, 0, f0, bin_counts)
      tiny_f = epsilon(1.0_dp) * max(1.0_dp, maxval(abs(f0)))
      select case (degree)
      case (0)
         call kde_drv_binned(x, h, minval(x), maxval(x), nb, weights, 4, f4)
         arg = 0.25_dp * f4
      case (1)
         call kde_drv_binned(x, h, minval(x), maxval(x), nb, weights, 1, f1)
         call kde_drv_binned(x, h, minval(x), maxval(x), nb, weights, 2, f2)
         do i = 1, nb + 1
            if (f0(i) <= tiny_f) return
            arg(i) = (0.5_dp * f2(i) + f1(i)**2 / f0(i))**2 / f0(i)
         end do
      case (2)
         call kde_drv_binned(x, h, minval(x), maxval(x), nb, weights, 1, f1)
         call kde_drv_binned(x, h, minval(x), maxval(x), nb, weights, 2, f2)
         call kde_drv_binned(x, h, minval(x), maxval(x), nb, weights, 4, f4)
         do i = 1, nb + 1
            if (f0(i) <= tiny_f) return
            arg(i) = f4(i) - 3.0_dp * f2(i)**2 / f0(i) + 2.0_dp * f1(i)**4 / f0(i)**3
            arg(i) = (0.125_dp * arg(i))**2 / f0(i)
         end do
      case default
         return
      end select

      if (sum(bin_counts) <= 0.0_dp) return
      ibias2 = sum(bin_counts * arg) / sum(bin_counts)
      ok = ieee_is_finite(ibias2) .and. ibias2 > 0.0_dp
   end subroutine integrated_bias2

   pure subroutine bkfe_bandwidth(x, weights, drv, scale, bandwidth, ok)
      real(dp), intent(in) :: x(:) !! Observations used by the kernel-functional bandwidth selector.
      real(dp), intent(in) :: weights(:) !! Weights normalized to mean one.
      integer, intent(in) :: drv !! Even derivative order of the target kernel functional.
      real(dp), intent(in) :: scale !! Positive normal-reference scale estimate.
      real(dp), intent(out) :: bandwidth !! Selected pilot bandwidth when ok is true.
      logical, intent(out) :: ok !! True when both normal-reference and plug-in stages are valid.

      integer, parameter :: nb = 400
      real(dp) :: bin_counts(nb + 1)
      real(dp) :: deriv(nb + 1)
      real(dp) :: h0
      real(dp) :: kr
      real(dp) :: n_eff
      real(dp) :: psi
      real(dp) :: ratio
      integer :: r

      ok = .false.
      bandwidth = 0.0_dp
      if (mod(drv, 2) /= 0 .or. .not. (scale > 0.0_dp)) return
      n_eff = sum(weights)**2 / sum(weights * weights)

      r = drv + 4
      psi = merge(1.0_dp, -1.0_dp, mod(r / 2, 2) == 0)
      psi = psi * factorial_real(r)
      psi = psi / ((2.0_dp * scale)**(r + 1) * factorial_real(r / 2) * sqrt(r_pi))
      kr = gaussian_derivative(0.0_dp, r - 2)
      ratio = -2.0_dp * kr / (psi * n_eff)
      if (.not. (ratio > 0.0_dp) .or. .not. ieee_is_finite(ratio)) return
      h0 = ratio**(1.0_dp / real(r + 1, dp))

      call kde_drv_binned(x, h0, minval(x), maxval(x), nb, weights, drv + 2, deriv, bin_counts)
      if (sum(bin_counts) <= 0.0_dp) return
      psi = sum(bin_counts * deriv) / sum(bin_counts)
      r = r - 2
      kr = gaussian_derivative(0.0_dp, r - 2)
      ratio = -2.0_dp * kr / (psi * n_eff)
      if (.not. (ratio > 0.0_dp) .or. .not. ieee_is_finite(ratio)) return
      bandwidth = ratio**(1.0_dp / real(r + 1, dp))
      ok = ieee_is_finite(bandwidth) .and. bandwidth > 0.0_dp
   end subroutine bkfe_bandwidth

   pure real(dp) function factorial_real(n) result(value)
      integer, intent(in) :: n !! Nonnegative integer whose factorial is returned as real(dp).

      integer :: i

      value = 1.0_dp
      do i = 2, n
         value = value * real(i, dp)
      end do
   end function factorial_real



end module kde1d_bandwidth
