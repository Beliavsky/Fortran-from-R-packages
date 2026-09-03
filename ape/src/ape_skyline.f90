! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Deterministic skyline computations translated from ape R/collapsed.intervals.R
! and R/skyline.R. Plotting and R class plumbing are intentionally omitted.
module ape_skyline
   use ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use ape_types, only : phylo_tree
   use ape_statistics, only : coalescent_intervals
   implicit none
   private

   type, public :: skyline_result
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: interval_length(:)
      real(dp), allocatable :: population_size(:)
      integer :: parameter_count = 0
      real(dp) :: epsilon = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: log_likelihood_aicc = 0.0_dp
   end type skyline_result

   public :: collapsed_intervals
   public :: skyline_from_intervals
   public :: skyline_tree
   public :: find_skyline_epsilon

contains

   pure subroutine collapsed_intervals(widths, epsilon, group, group_count, info)
      !! Reproduces ape `collapsed.intervals` grouping of coalescent intervals from the tips backward.
      real(dp), intent(in) :: widths(:) !! Nonnegative coalescent interval widths ordered from present toward the root.
      real(dp), intent(in) :: epsilon !! Nonnegative minimum accumulated width for a collapsed interval.
      integer, allocatable, intent(out) :: group(:) !! One-based collapsed-interval index assigned to each original interval.
      integer, intent(out) :: group_count !! Number of collapsed intervals after the trailing merge rule.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid widths or epsilon.
      real(dp) :: accumulated
      integer :: i
      integer :: p

      info = 0
      group_count = 0
      allocate(group(size(widths)))
      group = 0
      if (.not. ieee_is_finite(epsilon) .or. .not. all(ieee_is_finite(widths)) &
         .or. epsilon < 0.0_dp .or. any(widths < 0.0_dp)) then
         info = 1
         return
      end if
      if (size(widths) == 0) return
      if (epsilon > sum(widths)) then
         info = 2
         return
      end if

      p = 1
      accumulated = 0.0_dp
      do i = 1, size(widths)
         group(i) = p
         accumulated = accumulated + widths(i)
         if (accumulated >= epsilon) then
            p = p + 1
            accumulated = 0.0_dp
         end if
      end do
      if (sum(widths, mask=group == p) < epsilon) then
         p = p - 1
         where (group == p + 1) group = p
      end if
      group_count = group(size(widths))
   end subroutine collapsed_intervals

   pure subroutine skyline_from_intervals(lineages, widths, epsilon, result, info, old_style)
      !! Computes ape skyline population-size estimates and log-likelihood from coalescent intervals.
      integer, intent(in) :: lineages(:) !! Number of lineages at the start of each coalescent interval.
      real(dp), intent(in) :: widths(:) !! Coalescent interval widths corresponding one-to-one with `lineages`.
      real(dp), intent(in) :: epsilon !! Collapse threshold passed to ape `collapsed.intervals`.
      type(skyline_result), intent(out) :: result !! Skyline times, sizes, likelihood, AICc, and parameter count.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid intervals or degenerate population sizes.
      logical, intent(in), optional :: old_style !! Use ape's historical skyline estimator when true; default false.
      integer, allocatable :: group(:)
      real(dp), allocatable :: b(:)
      real(dp), allocatable :: classic(:)
      real(dp), allocatable :: ng(:)
      integer, allocatable :: count_group(:)
      logical :: use_old_style
      integer :: g
      integer :: groups
      integer :: i
      integer :: s

      result = skyline_result()
      info = 0
      s = size(widths)
      if (size(lineages) /= s .or. s == 0) then
         info = 1
         return
      end if
      if (any(lineages < 2)) then
         info = 2
         return
      end if
      call collapsed_intervals(widths, epsilon, group, groups, info)
      if (info /= 0) return
      if (groups <= 0) then
         info = 3
         return
      end if

      use_old_style = .false.
      if (present(old_style)) use_old_style = old_style
      allocate(result%time(groups), result%interval_length(groups), result%population_size(groups))
      allocate(ng(groups), count_group(groups), b(s), classic(s))
      result%interval_length = 0.0_dp
      result%population_size = 0.0_dp
      ng = 0.0_dp
      count_group = 0
      do i = 1, s
         b(i) = 0.5_dp * real(lineages(i), dp) * real(lineages(i) - 1, dp)
         classic(i) = widths(i) * b(i)
         g = group(i)
         result%interval_length(g) = result%interval_length(g) + widths(i)
         count_group(g) = count_group(g) + 1
         if (use_old_style) then
            if (count_group(g) == 1) ng(g) = real(lineages(i), dp)
         else
            ng(g) = ng(g) + classic(i)
         end if
      end do

      do g = 1, groups
         if (use_old_style) then
            result%population_size(g) = result%interval_length(g) * &
               ng(g) * (ng(g) - real(count_group(g), dp)) / (2.0_dp * real(count_group(g), dp))
         else
            result%population_size(g) = ng(g) / real(count_group(g), dp)
         end if
         if (result%population_size(g) <= 0.0_dp) then
            info = 4
            return
         end if
      end do

      result%time = cumulative_sum(result%interval_length)
      result%log_likelihood = 0.0_dp
      do i = 1, s
         g = group(i)
         result%log_likelihood = result%log_likelihood + log(b(i) / result%population_size(g)) &
            - b(i) * widths(i) / result%population_size(g)
      end do
      result%parameter_count = groups
      result%epsilon = epsilon
      if (s - groups > 1) then
         result%log_likelihood_aicc = result%log_likelihood - real(groups, dp) &
            - real(groups * (groups + 1), dp) / real(s - groups - 1, dp)
      else
         result%log_likelihood_aicc = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end subroutine skyline_from_intervals

   subroutine skyline_tree(tree, epsilon, result, info, old_style)
      !! Computes the deterministic skyline pipeline directly from a rooted binary phylogenetic tree.
      type(phylo_tree), intent(in) :: tree !! Rooted binary tree accepted by ape-style coalescent interval extraction.
      real(dp), intent(in) :: epsilon !! Nonnegative interval-collapse threshold.
      type(skyline_result), intent(out) :: result !! Skyline estimate derived from the tree's coalescent intervals.
      integer, intent(out) :: info !! Zero on success or a tree/skyline validation status code.
      logical, intent(in), optional :: old_style !! Select the historical estimator when true.
      integer, allocatable :: lineages(:)
      real(dp), allocatable :: widths(:)
      real(dp) :: total_depth

      call coalescent_intervals(tree, lineages, widths, total_depth, info)
      if (info /= 0) return
      if (present(old_style)) then
         call skyline_from_intervals(lineages, widths, epsilon, result, info, old_style)
      else
         call skyline_from_intervals(lineages, widths, epsilon, result, info)
      end if
   end subroutine skyline_tree

   subroutine find_skyline_epsilon(lineages, widths, best_epsilon, best_result, info, grid, minimum_epsilon, old_style)
      !! Reproduces ape's deterministic grid search for the AICc-optimal skyline collapse threshold.
      integer, intent(in) :: lineages(:) !! Coalescent lineage counts corresponding to `widths`.
      real(dp), intent(in) :: widths(:) !! Nonnegative coalescent interval widths.
      real(dp), intent(out) :: best_epsilon !! Grid-search epsilon giving the largest eligible AICc log-likelihood.
      type(skyline_result), intent(out) :: best_result !! Skyline result evaluated at `best_epsilon`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid grid settings or skyline evaluation failure.
      integer, intent(in), optional :: grid !! Number of equal grid divisions across total depth; default 1000.
      real(dp), intent(in), optional :: minimum_epsilon !! Strict lower search limit; default `1e-6`.
      logical, intent(in), optional :: old_style !! Select the historical skyline estimator when true.
      type(skyline_result) :: candidate
      real(dp) :: best_aicc
      real(dp) :: delta
      real(dp) :: eps
      real(dp) :: minimum_value
      real(dp) :: total_depth
      integer :: grid_size
      integer :: status
      logical :: use_old_style

      info = 0
      best_result = skyline_result()
      if (size(widths) == 0 .or. size(lineages) /= size(widths)) then
         info = 1
         return
      end if
      grid_size = 1000
      if (present(grid)) grid_size = grid
      minimum_value = 1.0e-6_dp
      if (present(minimum_epsilon)) minimum_value = minimum_epsilon
      if (grid_size <= 0 .or. minimum_value < 0.0_dp) then
         info = 2
         return
      end if
      use_old_style = .false.
      if (present(old_style)) use_old_style = old_style
      total_depth = sum(widths)
      if (total_depth <= 0.0_dp) then
         info = 3
         return
      end if

      best_epsilon = total_depth
      call skyline_from_intervals(lineages, widths, best_epsilon, best_result, status, use_old_style)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      best_aicc = best_result%log_likelihood_aicc
      delta = total_depth / real(grid_size, dp)
      eps = total_depth - delta
      do while (eps > minimum_value)
         call skyline_from_intervals(lineages, widths, eps, candidate, status, use_old_style)
         if (status == 0 .and. .not. ieee_is_nan(candidate%log_likelihood_aicc)) then
            if ((ieee_is_nan(best_aicc) .or. candidate%log_likelihood_aicc > best_aicc) &
               .and. candidate%parameter_count < size(widths) - 1) then
               best_epsilon = eps
               best_aicc = candidate%log_likelihood_aicc
               best_result = candidate
            end if
         end if
         eps = eps - delta
      end do
   end subroutine find_skyline_epsilon

   pure function cumulative_sum(values) result(cumulative)
      !! Returns the cumulative sum of a real vector.
      real(dp), intent(in) :: values(:) !! Input values in their natural order.
      real(dp) :: cumulative(size(values))
      integer :: i

      if (size(values) == 0) return
      cumulative(1) = values(1)
      do i = 2, size(values)
         cumulative(i) = cumulative(i - 1) + values(i)
      end do
   end function cumulative_sum

end module ape_skyline
