! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Extended birth-death likelihood translated from ape R/birthdeath.R `bd.ext`.
! Upstream copyright and provenance are documented in NOTICE.md.
module ape_birthdeath_extended
   use r_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ape_types, only : phylo_tree
   use ape_tree_algorithms, only : branching_times
   implicit none
   private

   type, public :: birthdeath_extended_result
      integer :: n_tip = 0
      logical :: conditional = .true.
      real(dp) :: death_birth_ratio = 0.0_dp
      real(dp) :: net_diversification = 0.0_dp
      real(dp) :: death_birth_ratio_se = 0.0_dp
      real(dp) :: net_diversification_se = 0.0_dp
      real(dp) :: deviance = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
   end type birthdeath_extended_result

   public :: birthdeath_extended_fit
   public :: birthdeath_extended_from_data
   public :: birthdeath_extended_deviance

contains

   pure subroutine birthdeath_extended_fit(tree, species_count, result, info, conditional)
      !! Fits ape `bd.ext` to a rooted tree and per-tip taxonomic richness counts.
      type(phylo_tree), intent(in) :: tree !! Rooted binary tree with positive branch lengths and root numbered `n_tip+1`.
      integer, intent(in) :: species_count(:) !! Taxonomic richness represented by each tip; one positive count per tip.
      type(birthdeath_extended_result), intent(out) :: result !! Fitted ratio/rate, Hessian SEs, deviance, and log-likelihood.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid tree, richness counts, or numerical failure.
      logical, intent(in), optional :: conditional !! Use ape's conditional likelihood when true; default is true.
      real(dp), allocatable :: terminal_length(:)
      real(dp), allocatable :: times(:)
      logical :: use_conditional
      integer :: child
      integer :: e
      integer :: status

      result = birthdeath_extended_result()
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
      if (size(species_count) /= tree%n_tip .or. any(species_count < 1)) then
         info = 4
         return
      end if
      if (.not. all(ieee_is_finite(tree%edge_length)) .or. any(tree%edge_length <= 0.0_dp)) then
         info = 5
         return
      end if
      allocate(terminal_length(tree%n_tip))
      terminal_length = -1.0_dp
      do e = 1, tree%nedge()
         child = tree%edge(e, 2)
         if (child <= tree%n_tip) terminal_length(child) = tree%edge_length(e)
      end do
      if (any(terminal_length <= 0.0_dp)) then
         info = 6
         return
      end if
      call branching_times(tree, times, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      call birthdeath_extended_from_data(times, terminal_length, species_count, result, info, use_conditional)
   end subroutine birthdeath_extended_fit

   pure subroutine birthdeath_extended_from_data(times, terminal_length, species_count, result, info, conditional)
      !! Fits ape `bd.ext` from internal-node ages, terminal lengths, and taxonomic richness counts.
      real(dp), intent(in) :: times(:) !! Internal-node branching ages; length must be one less than the number of tips.
      real(dp), intent(in) :: terminal_length(:) !! Positive terminal-branch lengths, one value per tip.
      integer, intent(in) :: species_count(:) !! Positive taxonomic richness counts, one value per terminal branch.
      type(birthdeath_extended_result), intent(out) :: result !! Fitted ratio/rate, Hessian SEs, deviance, and log-likelihood.
      integer, intent(out) :: info !! Zero on success; nonzero for inconsistent inputs or failed optimization/Hessian.
      logical, intent(in), optional :: conditional !! Use ape's conditional likelihood when true; default is true.
      real(dp) :: dev_best
      real(dp) :: dev_grid
      real(dp) :: dev_refined
      real(dp) :: ratio_best
      real(dp) :: ratio_grid
      real(dp) :: ratio_left
      real(dp) :: ratio_right
      real(dp) :: rate_best
      real(dp) :: rate_grid
      real(dp) :: rate_refined
      real(dp) :: step
      logical :: use_conditional
      integer :: best_index
      integer :: i
      integer, parameter :: grid_count = 101

      result = birthdeath_extended_result()
      info = 0
      use_conditional = .true.
      if (present(conditional)) use_conditional = conditional
      if (size(species_count) < 3) then
         info = 1
         return
      end if
      if (size(terminal_length) /= size(species_count) .or. size(times) /= size(species_count) - 1) then
         info = 2
         return
      end if
      if (any(species_count < 1)) then
         info = 3
         return
      end if
      if (.not. all(ieee_is_finite(times)) .or. any(times < 0.0_dp)) then
         info = 4
         return
      end if
      if (.not. all(ieee_is_finite(terminal_length)) .or. any(terminal_length <= 0.0_dp)) then
         info = 5
         return
      end if

      dev_best = huge(1.0_dp)
      ratio_best = 0.0_dp
      rate_best = 0.0_dp
      best_index = 1
      step = 0.99_dp / real(grid_count - 1, dp)
      do i = 1, grid_count
         ratio_grid = real(i - 1, dp) * step
         call optimize_extended_rate(times, terminal_length, species_count, ratio_grid, use_conditional, &
            rate_grid, dev_grid)
         if (dev_grid < dev_best) then
            dev_best = dev_grid
            ratio_best = ratio_grid
            rate_best = rate_grid
            best_index = i
         end if
      end do

      if (best_index > 1) then
         ratio_left = real(best_index - 2, dp) * step
         ratio_right = min(1.0_dp - 1.0e-12_dp, real(best_index, dp) * step)
         call refine_extended_ratio(times, terminal_length, species_count, ratio_left, ratio_right, &
            use_conditional, ratio_best, rate_refined, dev_refined)
         if (dev_refined < dev_best) then
            rate_best = rate_refined
            dev_best = dev_refined
         end if
      else
         ratio_best = 0.0_dp
         call optimize_extended_rate(times, terminal_length, species_count, ratio_best, use_conditional, &
            rate_best, dev_best)
      end if

      if (.not. ieee_is_finite(dev_best) .or. rate_best <= 0.0_dp) then
         info = 6
         return
      end if
      if (ratio_best < 5.0e-8_dp) ratio_best = 0.0_dp
      result%n_tip = size(species_count)
      result%conditional = use_conditional
      result%death_birth_ratio = ratio_best
      result%net_diversification = rate_best
      result%deviance = dev_best
      result%log_likelihood = -0.5_dp * dev_best
      call extended_standard_errors(times, terminal_length, species_count, ratio_best, rate_best, &
         use_conditional, result%death_birth_ratio_se, result%net_diversification_se, info)
   end subroutine birthdeath_extended_from_data

   pure real(dp) function birthdeath_extended_deviance(times, terminal_length, species_count, &
      death_birth_ratio, net_diversification, conditional) result(value)
      !! Evaluates the exact conditional or unconditional deviance used by ape `bd.ext`.
      real(dp), intent(in) :: times(:) !! Internal-node branching ages; one fewer value than the number of tips.
      real(dp), intent(in) :: terminal_length(:) !! Positive terminal-branch lengths, one per tip.
      integer, intent(in) :: species_count(:) !! Positive taxonomic richness counts, one per tip.
      real(dp), intent(in) :: death_birth_ratio !! Extinction/speciation ratio `d/b`, constrained to `[0,1)`.
      real(dp), intent(in) :: net_diversification !! Positive net diversification rate `b-d`.
      logical, intent(in) :: conditional !! Selects ape's conditional likelihood when true.
      real(dp) :: exp_negative
      real(dp) :: log_denom
      real(dp) :: log_likelihood
      real(dp) :: log_one_minus_exp
      real(dp) :: z
      integer :: i
      integer :: n_tip

      value = huge(1.0_dp)
      n_tip = size(species_count)
      if (n_tip < 3) return
      if (size(times) /= n_tip - 1 .or. size(terminal_length) /= n_tip) return
      if (any(species_count < 1)) return
      if (.not. all(ieee_is_finite(times)) .or. any(times < 0.0_dp)) return
      if (.not. all(ieee_is_finite(terminal_length)) .or. any(terminal_length <= 0.0_dp)) return
      if (.not. ieee_is_finite(death_birth_ratio) .or. .not. ieee_is_finite(net_diversification)) return
      if (death_birth_ratio < 0.0_dp .or. death_birth_ratio >= 1.0_dp) return
      if (net_diversification <= 0.0_dp) return

      log_likelihood = log_gamma(real(n_tip, dp))
      log_likelihood = log_likelihood + real(n_tip - 2, dp) * log(net_diversification)
      if (conditional) then
         log_likelihood = log_likelihood + real(n_tip, dp) * log(1.0_dp - death_birth_ratio)
      else
         log_likelihood = log_likelihood + real(3 * n_tip, dp) * log(1.0_dp - death_birth_ratio)
      end if
      log_likelihood = log_likelihood + 2.0_dp * net_diversification * sum(times)
      do i = 1, size(times)
         z = net_diversification * times(i)
         exp_negative = exp(-z)
         if (1.0_dp - death_birth_ratio * exp_negative <= 0.0_dp) return
         log_denom = z + log(1.0_dp - death_birth_ratio * exp_negative)
         log_likelihood = log_likelihood - 2.0_dp * log_denom
      end do

      do i = 1, n_tip
         z = net_diversification * terminal_length(i)
         exp_negative = exp(-z)
         if (exp_negative >= 1.0_dp) return
         if (1.0_dp - death_birth_ratio * exp_negative <= 0.0_dp) return
         log_one_minus_exp = log(1.0_dp - exp_negative)
         log_denom = log(1.0_dp - death_birth_ratio * exp_negative)
         if (conditional) then
            log_likelihood = log_likelihood + log(1.0_dp - death_birth_ratio) - z - log_denom
            log_likelihood = log_likelihood + real(species_count(i) - 1, dp) * (log_one_minus_exp - log_denom)
         else
            log_likelihood = log_likelihood + z
            log_likelihood = log_likelihood + real(species_count(i) - 1, dp) * (z + log_one_minus_exp)
            log_likelihood = log_likelihood - real(species_count(i) + 1, dp) * (z + log_denom)
         end if
      end do
      value = -2.0_dp * log_likelihood
   end function birthdeath_extended_deviance

   pure subroutine optimize_extended_rate(times, terminal_length, species_count, ratio, conditional, rate, deviance)
      !! Minimizes the extended birth-death deviance over positive net rate for a fixed ratio.
      real(dp), intent(in) :: times(:) !! Internal-node branching ages.
      real(dp), intent(in) :: terminal_length(:) !! Positive terminal-branch lengths.
      integer, intent(in) :: species_count(:) !! Positive per-tip taxonomic richness counts.
      real(dp), intent(in) :: ratio !! Fixed death/birth ratio in `[0,1)`.
      logical, intent(in) :: conditional !! Selects the conditional or unconditional ape likelihood.
      real(dp), intent(out) :: rate !! Positive minimizing net diversification rate.
      real(dp), intent(out) :: deviance !! Deviance at the minimizing rate.
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      real(dp) :: log_scale
      real(dp) :: scale_sum
      real(dp), parameter :: golden = 0.6180339887498948482_dp
      integer :: iter

      scale_sum = sum(times) + sum(terminal_length)
      log_scale = -log(max(scale_sum / real(size(times) + size(terminal_length), dp), sqrt(tiny(1.0_dp))))
      a = log_scale - 20.0_dp
      b = log_scale + 20.0_dp
      c = b - golden * (b - a)
      d = a + golden * (b - a)
      fc = birthdeath_extended_deviance(times, terminal_length, species_count, ratio, exp(c), conditional)
      fd = birthdeath_extended_deviance(times, terminal_length, species_count, ratio, exp(d), conditional)
      do iter = 1, 160
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - golden * (b - a)
            fc = birthdeath_extended_deviance(times, terminal_length, species_count, ratio, exp(c), conditional)
         else
            a = c
            c = d
            fc = fd
            d = a + golden * (b - a)
            fd = birthdeath_extended_deviance(times, terminal_length, species_count, ratio, exp(d), conditional)
         end if
      end do
      if (fc <= fd) then
         rate = exp(c)
         deviance = fc
      else
         rate = exp(d)
         deviance = fd
      end if
   end subroutine optimize_extended_rate

   pure subroutine refine_extended_ratio(times, terminal_length, species_count, lower, upper, conditional, &
      ratio, rate, deviance)
      !! Refines the extended-model ratio while reoptimizing the net diversification rate.
      real(dp), intent(in) :: times(:) !! Internal-node branching ages.
      real(dp), intent(in) :: terminal_length(:) !! Positive terminal-branch lengths.
      integer, intent(in) :: species_count(:) !! Positive per-tip taxonomic richness counts.
      real(dp), intent(in) :: lower !! Lower ratio bound from the coarse profile search.
      real(dp), intent(in) :: upper !! Upper ratio bound from the coarse profile search; strictly below one.
      logical, intent(in) :: conditional !! Selects the conditional or unconditional ape likelihood.
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
      call optimize_extended_rate(times, terminal_length, species_count, c, conditional, rc, fc)
      call optimize_extended_rate(times, terminal_length, species_count, d, conditional, rd, fd)
      do iter = 1, 120
         if (fc <= fd) then
            b = d
            d = c
            rd = rc
            fd = fc
            c = b - golden * (b - a)
            call optimize_extended_rate(times, terminal_length, species_count, c, conditional, rc, fc)
         else
            a = c
            c = d
            rc = rd
            fc = fd
            d = a + golden * (b - a)
            call optimize_extended_rate(times, terminal_length, species_count, d, conditional, rd, fd)
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
   end subroutine refine_extended_ratio

   pure subroutine extended_standard_errors(times, terminal_length, species_count, ratio, rate, conditional, &
      ratio_se, rate_se, info)
      !! Computes local Hessian standard errors for the extended birth-death fit.
      real(dp), intent(in) :: times(:) !! Internal-node branching ages.
      real(dp), intent(in) :: terminal_length(:) !! Positive terminal-branch lengths.
      integer, intent(in) :: species_count(:) !! Positive per-tip taxonomic richness counts.
      real(dp), intent(in) :: ratio !! Fitted death/birth ratio.
      real(dp), intent(in) :: rate !! Fitted net diversification rate.
      logical, intent(in) :: conditional !! Selects the conditional or unconditional ape likelihood.
      real(dp), intent(out) :: ratio_se !! Standard error for `d/b`; zero when the optimum is on the zero boundary.
      real(dp), intent(out) :: rate_se !! Standard error for `b-d`.
      integer, intent(out) :: info !! Zero on success; nonzero when the numerical Hessian is not positive definite.
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
      f0 = birthdeath_extended_deviance(times, terminal_length, species_count, ratio, rate, conditional)
      hr = max(1.0e-6_dp, 1.0e-4_dp * rate)
      if (rate - hr <= 0.0_dp) hr = 0.25_dp * rate
      frr = (birthdeath_extended_deviance(times, terminal_length, species_count, ratio, rate + hr, conditional) &
         - 2.0_dp * f0 &
         + birthdeath_extended_deviance(times, terminal_length, species_count, ratio, rate - hr, conditional)) &
         / (hr * hr)
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
      faa = (birthdeath_extended_deviance(times, terminal_length, species_count, ratio + ha, rate, conditional) &
         - 2.0_dp * f0 &
         + birthdeath_extended_deviance(times, terminal_length, species_count, ratio - ha, rate, conditional)) &
         / (ha * ha)
      far = (birthdeath_extended_deviance(times, terminal_length, species_count, ratio + ha, rate + hr, conditional) &
         - birthdeath_extended_deviance(times, terminal_length, species_count, ratio + ha, rate - hr, conditional) &
         - birthdeath_extended_deviance(times, terminal_length, species_count, ratio - ha, rate + hr, conditional) &
         + birthdeath_extended_deviance(times, terminal_length, species_count, ratio - ha, rate - hr, conditional)) &
         / (4.0_dp * ha * hr)
      det = faa * frr - far * far
      if (faa <= 0.0_dp .or. frr <= 0.0_dp .or. det <= 0.0_dp .or. .not. ieee_is_finite(det)) then
         info = 3
         return
      end if
      ratio_se = sqrt(frr / det)
      rate_se = sqrt(faa / det)
   end subroutine extended_standard_errors

end module ape_birthdeath_extended
