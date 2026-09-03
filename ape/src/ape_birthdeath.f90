! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Standard birth-death likelihood translated from ape R/birthdeath.R.
! Upstream copyright and provenance are documented in NOTICE.md.
module ape_birthdeath
   use r_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ape_types, only : phylo_tree
   use ape_tree_algorithms, only : branching_times
   implicit none
   private

   type, public :: birthdeath_result
      integer :: n_tip = 0
      real(dp) :: death_birth_ratio = 0.0_dp
      real(dp) :: net_diversification = 0.0_dp
      real(dp) :: death_birth_ratio_se = 0.0_dp
      real(dp) :: net_diversification_se = 0.0_dp
      real(dp) :: deviance = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: confidence_interval(2, 2) = 0.0_dp
   end type birthdeath_result

   public :: birthdeath_fit
   public :: birthdeath_from_times
   public :: birthdeath_deviance

contains

   pure subroutine birthdeath_fit(tree, result, info)
      !! Fits ape's standard two-parameter birth-death model to a rooted tree.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths and root numbered `n_tip+1`.
      type(birthdeath_result), intent(out) :: result !! Fitted ratio/net-rate, SEs, deviance, and fixed-profile 95% intervals.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid tree/times or numerical failure.
      real(dp), allocatable :: times(:)
      integer :: status

      result = birthdeath_result()
      if (.not. tree%valid() .or. .not. tree%has_lengths()) then
         info = 1
         return
      end if
      if (tree%n_tip < 3 .or. tree%n_node /= tree%n_tip - 1) then
         info = 2
         return
      end if
      if (tree%root() /= tree%n_tip + 1) then
         info = 3
         return
      end if
      call branching_times(tree, times, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call birthdeath_from_times(times, result, info)
   end subroutine birthdeath_fit

   pure subroutine birthdeath_from_times(times, result, info)
      !! Fits the standard birth-death model from ape-style internal-node branching times.
      real(dp), intent(in) :: times(:) !! Internal-node ages with root age first; length is `n_tip-1`.
      type(birthdeath_result), intent(out) :: result !! Fitted parameters, SEs, likelihood, and fixed-profile intervals.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid ages or failed optimization/Hessian.
      real(dp) :: a_best
      real(dp) :: a_grid
      real(dp) :: a_left
      real(dp) :: a_right
      real(dp) :: dev_best
      real(dp) :: dev_grid
      real(dp) :: dev_refined
      real(dp) :: r_best
      real(dp) :: r_grid
      real(dp) :: r_refined
      real(dp) :: step
      integer :: best_index
      integer :: i
      integer, parameter :: grid_count = 101

      result = birthdeath_result()
      info = 0
      if (size(times) < 2) then
         info = 1
         return
      end if
      if (.not. all(ieee_is_finite(times)) .or. any(times < 0.0_dp)) then
         info = 2
         return
      end if
      if (times(1) <= 0.0_dp) then
         info = 3
         return
      end if

      dev_best = huge(1.0_dp)
      a_best = 0.0_dp
      r_best = 0.0_dp
      best_index = 1
      step = 0.99_dp / real(grid_count - 1, dp)
      do i = 1, grid_count
         a_grid = real(i - 1, dp) * step
         call optimize_net_rate(times, a_grid, r_grid, dev_grid)
         if (dev_grid < dev_best) then
            dev_best = dev_grid
            a_best = a_grid
            r_best = r_grid
            best_index = i
         end if
      end do

      if (best_index > 1) then
         a_left = real(best_index - 2, dp) * step
         a_right = min(1.0_dp - 1.0e-12_dp, real(best_index, dp) * step)
         call refine_ratio(times, a_left, a_right, a_best, r_refined, dev_refined)
         if (dev_refined < dev_best) then
            r_best = r_refined
            dev_best = dev_refined
         end if
      else
         a_best = 0.0_dp
         call optimize_net_rate(times, a_best, r_best, dev_best)
      end if

      if (.not. ieee_is_finite(dev_best) .or. r_best <= 0.0_dp) then
         info = 4
         return
      end if
      if (a_best < 5.0e-8_dp) a_best = 0.0_dp
      result%n_tip = size(times) + 1
      result%death_birth_ratio = a_best
      result%net_diversification = r_best
      result%deviance = dev_best
      result%log_likelihood = -0.5_dp * dev_best
      call birthdeath_standard_errors(times, a_best, r_best, result%death_birth_ratio_se, &
         result%net_diversification_se, info)
      if (info /= 0) return
      call birthdeath_profile_intervals(times, a_best, r_best, dev_best, result%confidence_interval)
   end subroutine birthdeath_from_times

   pure real(dp) function birthdeath_deviance(times, death_birth_ratio, net_diversification) result(value)
      !! Evaluates the exact deviance used by ape `birthdeath` for fixed parameters.
      real(dp), intent(in) :: times(:) !! Internal-node ages with the root first.
      real(dp), intent(in) :: death_birth_ratio !! Extinction/speciation ratio `d/b`; values below one are admissible.
      real(dp), intent(in) :: net_diversification !! Net diversification rate `b-d`; must be positive.
      real(dp) :: log_likelihood
      real(dp) :: z
      integer :: i
      integer :: n_tip

      value = huge(1.0_dp)
      if (size(times) < 2) return
      if (.not. all(ieee_is_finite(times))) return
      if (.not. ieee_is_finite(death_birth_ratio) .or. .not. ieee_is_finite(net_diversification)) return
      if (death_birth_ratio >= 1.0_dp .or. net_diversification <= 0.0_dp) return
      if (1.0_dp - death_birth_ratio <= 0.0_dp) return
      n_tip = size(times) + 1
      log_likelihood = log_gamma(real(n_tip, dp))
      log_likelihood = log_likelihood + real(n_tip - 2, dp) * log(net_diversification)
      log_likelihood = log_likelihood + net_diversification * sum(times(2:))
      log_likelihood = log_likelihood + real(n_tip, dp) * log(1.0_dp - death_birth_ratio)
      do i = 1, size(times)
         z = net_diversification * times(i)
         if (1.0_dp - death_birth_ratio * exp(-z) <= 0.0_dp) return
         log_likelihood = log_likelihood - 2.0_dp * (z + log(1.0_dp - death_birth_ratio * exp(-z)))
      end do
      value = -2.0_dp * log_likelihood
   end function birthdeath_deviance

   pure subroutine optimize_net_rate(times, ratio, rate, deviance)
      !! Minimizes the birth-death deviance over positive net diversification for fixed `d/b`.
      real(dp), intent(in) :: times(:) !! Internal-node ages with root first.
      real(dp), intent(in) :: ratio !! Fixed death/birth ratio below one.
      real(dp), intent(out) :: rate !! Positive minimizing net diversification rate.
      real(dp), intent(out) :: deviance !! Deviance at `rate`.
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      real(dp) :: log_scale
      real(dp), parameter :: golden = 0.6180339887498948482_dp
      integer :: iter

      log_scale = -log(max(sum(times) / real(size(times), dp), sqrt(tiny(1.0_dp))))
      a = log_scale - 20.0_dp
      b = log_scale + 20.0_dp
      c = b - golden * (b - a)
      d = a + golden * (b - a)
      fc = birthdeath_deviance(times, ratio, exp(c))
      fd = birthdeath_deviance(times, ratio, exp(d))
      do iter = 1, 160
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - golden * (b - a)
            fc = birthdeath_deviance(times, ratio, exp(c))
         else
            a = c
            c = d
            fc = fd
            d = a + golden * (b - a)
            fd = birthdeath_deviance(times, ratio, exp(d))
         end if
      end do
      if (fc <= fd) then
         rate = exp(c)
         deviance = fc
      else
         rate = exp(d)
         deviance = fd
      end if
   end subroutine optimize_net_rate

   pure subroutine refine_ratio(times, lower, upper, ratio, rate, deviance)
      !! Refines the best interior `d/b` ratio while reoptimizing the net diversification rate.
      real(dp), intent(in) :: times(:) !! Internal-node ages with root first.
      real(dp), intent(in) :: lower !! Lower ratio bound from the coarse search.
      real(dp), intent(in) :: upper !! Upper ratio bound from the coarse search; less than one.
      real(dp), intent(out) :: ratio !! Refined death/birth ratio.
      real(dp), intent(out) :: rate !! Net diversification rate optimized at `ratio`.
      real(dp), intent(out) :: deviance !! Deviance at the refined optimum.
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      real(dp) :: rc
      real(dp) :: rd
      real(dp), parameter :: golden = 0.6180339887498948482_dp
      integer :: iter

      a = lower
      b = upper
      c = b - golden * (b - a)
      d = a + golden * (b - a)
      call optimize_net_rate(times, c, rc, fc)
      call optimize_net_rate(times, d, rd, fd)
      do iter = 1, 120
         if (fc <= fd) then
            b = d
            d = c
            rd = rc
            fd = fc
            c = b - golden * (b - a)
            call optimize_net_rate(times, c, rc, fc)
         else
            a = c
            c = d
            rc = rd
            fc = fd
            d = a + golden * (b - a)
            call optimize_net_rate(times, d, rd, fd)
         end if
      end do
      if (fc <= fd) then
         ratio = c
         rate = rc
         deviance = fc
      else
         ratio = d
         rate = rd
         deviance = fd
      end if
   end subroutine refine_ratio

   pure subroutine birthdeath_standard_errors(times, ratio, rate, ratio_se, rate_se, info)
      !! Computes ape-compatible local standard errors from the Hessian of the deviance.
      real(dp), intent(in) :: times(:) !! Internal-node ages with root first.
      real(dp), intent(in) :: ratio !! Fitted death/birth ratio.
      real(dp), intent(in) :: rate !! Fitted net diversification rate.
      real(dp), intent(out) :: ratio_se !! Standard error for `d/b`, zero at the boundary fit.
      real(dp), intent(out) :: rate_se !! Standard error for `b-d`.
      integer, intent(out) :: info !! Zero on success; nonzero if the numerical Hessian is not positive definite.
      real(dp) :: det
      real(dp) :: f0
      real(dp) :: faa
      real(dp) :: far
      real(dp) :: frr
      real(dp) :: ha
      real(dp) :: hr

      info = 0
      ratio_se = 0.0_dp
      rate_se = 0.0_dp
      f0 = birthdeath_deviance(times, ratio, rate)
      hr = max(1.0e-6_dp, 1.0e-4_dp * rate)
      if (rate - hr <= 0.0_dp) hr = 0.25_dp * rate
      frr = (birthdeath_deviance(times, ratio, rate + hr) - 2.0_dp * f0 &
         + birthdeath_deviance(times, ratio, rate - hr)) / (hr * hr)
      if (ratio <= 0.0_dp) then
         if (frr <= 0.0_dp .or. .not. ieee_is_finite(frr)) then
            info = 1
            return
         end if
         rate_se = sqrt(1.0_dp / frr)
         return
      end if

      ha = max(1.0e-6_dp, 1.0e-4_dp * max(1.0_dp, abs(ratio)))
      ha = min(ha, 0.25_dp * ratio, 0.25_dp * (1.0_dp - ratio))
      if (ha <= 0.0_dp) then
         info = 2
         return
      end if
      faa = (birthdeath_deviance(times, ratio + ha, rate) - 2.0_dp * f0 &
         + birthdeath_deviance(times, ratio - ha, rate)) / (ha * ha)
      far = (birthdeath_deviance(times, ratio + ha, rate + hr) &
         - birthdeath_deviance(times, ratio + ha, rate - hr) &
         - birthdeath_deviance(times, ratio - ha, rate + hr) &
         + birthdeath_deviance(times, ratio - ha, rate - hr)) / (4.0_dp * ha * hr)
      det = faa * frr - far * far
      if (faa <= 0.0_dp .or. frr <= 0.0_dp .or. det <= 0.0_dp &
         .or. .not. ieee_is_finite(det)) then
         info = 3
         return
      end if
      ratio_se = sqrt(frr / det)
      rate_se = sqrt(faa / det)
   end subroutine birthdeath_standard_errors

   pure subroutine birthdeath_profile_intervals(times, ratio, rate, deviance, interval)
      !! Reproduces ape's fixed-other-parameter 3.84-deviance confidence-interval stepping rule.
      real(dp), intent(in) :: times(:) !! Internal-node ages with root first.
      real(dp), intent(in) :: ratio !! Fitted death/birth ratio.
      real(dp), intent(in) :: rate !! Fitted net diversification rate.
      real(dp), intent(in) :: deviance !! Minimum deviance used as the interval reference.
      real(dp), intent(out) :: interval(2, 2) !! Rows `(d/b,b-d)`, columns lower and upper limits.

      interval(1, 1) = profile_bound(times, ratio, rate, deviance, 1, -1)
      interval(2, 1) = profile_bound(times, ratio, rate, deviance, 2, -1)
      interval(1, 2) = profile_bound(times, ratio, rate, deviance, 1, 1)
      interval(2, 2) = profile_bound(times, ratio, rate, deviance, 2, 1)
   end subroutine birthdeath_profile_intervals

   pure real(dp) function profile_bound(times, ratio, rate, deviance, which_parameter, direction) result(bound)
      !! Applies the decimal-step search used by ape for one fixed-parameter confidence bound.
      real(dp), intent(in) :: times(:) !! Internal-node ages with root first.
      real(dp), intent(in) :: ratio !! Fitted death/birth ratio.
      real(dp), intent(in) :: rate !! Fitted net diversification rate.
      real(dp), intent(in) :: deviance !! Minimum deviance.
      integer, intent(in) :: which_parameter !! One for `d/b`, two for `b-d`.
      integer, intent(in) :: direction !! Minus one for lower, plus one for upper.
      real(dp) :: increment
      real(dp) :: p
      real(dp) :: trial
      integer :: guard

      increment = 0.1_dp
      if (which_parameter == 1) then
         p = ratio + real(direction, dp) * increment
      else
         p = rate + real(direction, dp) * increment
      end if
      do while (increment > 1.0e-9_dp)
         guard = 0
         do
            if (which_parameter == 1) then
               trial = birthdeath_deviance(times, p, rate)
            else
               trial = birthdeath_deviance(times, ratio, p)
            end if
            if (trial >= deviance + 3.84_dp) exit
            p = p + real(direction, dp) * increment
            guard = guard + 1
            if (guard > 100000) exit
         end do
         p = p - real(direction, dp) * increment
         increment = increment / 10.0_dp
      end do
      bound = p
   end function profile_bound

end module ape_birthdeath
