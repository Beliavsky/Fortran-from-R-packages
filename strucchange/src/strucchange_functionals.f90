! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_functionals
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq, r_qchisq
   use strucchange_pvalues, only : efp_pvalue, fstats_pvalue
   implicit none
   private
   public :: cat_l2_bb_critical_value
   public :: cat_l2_bb_pvalue
   public :: cat_l2_bb_statistic
   public :: max_mosum_critical_value
   public :: max_mosum_pvalue
   public :: max_mosum_statistic
   public :: sup_lm_critical_value
   public :: sup_lm_pvalue
   public :: sup_lm_statistic
contains
   real(dp) function sup_lm_statistic(process, from, to) result(statistic)
      real(dp), intent(in) :: process(:, :)
      real(dp), intent(in) :: from
      real(dp), intent(in), optional :: to
      real(dp) :: to_use, t, value
      integer :: first_row, i, i1, i2, n

      statistic = -huge(1.0_dp)
      if (size(process, 1) < 2 .or. size(process, 2) < 1) return
      first_row = process_start_row(process)
      n = size(process, 1) - first_row + 1
      if (n < 2 .or. from <= 0.0_dp .or. from >= 1.0_dp) return
      to_use = 1.0_dp - from
      if (present(to)) to_use = to
      if (to_use <= from .or. to_use >= 1.0_dp) return

      i1 = max(1, floor(from * real(n, dp)))
      i2 = min(n, floor(to_use * real(n, dp)))
      do i = i1, i2
         t = real(i, dp) / real(n, dp)
         value = sum(process(first_row + i - 1, :) ** 2) / (t * (1.0_dp - t))
         statistic = max(statistic, value)
      end do
   end function sup_lm_statistic

   real(dp) function sup_lm_pvalue(statistic, nproc, from, to) result(p)
      real(dp), intent(in) :: statistic, from
      real(dp), intent(in), optional :: to
      integer, intent(in) :: nproc
      real(dp) :: lambda, to_use

      to_use = 1.0_dp - from
      if (present(to)) to_use = to
      if (from <= 0.0_dp .or. to_use <= from .or. to_use >= 1.0_dp) then
         p = 1.0_dp
         return
      end if
      lambda = ((1.0_dp - from) * to_use) / (from * (1.0_dp - to_use))
      p = fstats_pvalue(statistic, "supF", nproc, lambda)
   end function sup_lm_pvalue

   real(dp) function sup_lm_critical_value(alpha, nproc, from, to) result(value)
      real(dp), intent(in) :: alpha, from
      real(dp), intent(in), optional :: to
      integer, intent(in) :: nproc
      real(dp) :: low, high, midpoint, p, to_use
      integer :: iteration

      to_use = 1.0_dp - from
      if (present(to)) to_use = to
      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      low = 0.0_dp
      high = 1000.0_dp
      do iteration = 1, 160
         midpoint = 0.5_dp * (low + high)
         p = sup_lm_pvalue(midpoint, nproc, from, to_use)
         if (p > alpha) then
            low = midpoint
         else
            high = midpoint
         end if
         if (high - low <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, midpoint)) exit
      end do
      value = 0.5_dp * (low + high)
   end function sup_lm_critical_value

   real(dp) function max_mosum_statistic(process, width) result(statistic)
      real(dp), intent(in) :: process(:, :)
      real(dp), intent(in) :: width
      integer :: first_row, i, j, n, nh
      real(dp) :: difference

      statistic = -huge(1.0_dp)
      if (size(process, 1) < 2 .or. size(process, 2) < 1) return
      first_row = process_start_row(process)
      n = size(process, 1) - first_row + 1
      if (width < 1.0_dp) then
         nh = floor(real(n, dp) * width)
      else
         nh = nint(width)
      end if
      if (nh < 1 .or. nh > n) return

      do j = 1, size(process, 2)
         do i = 0, n - nh
            difference = process_value(process, first_row, i + nh, j) - &
               process_value(process, first_row, i, j)
            statistic = max(statistic, abs(difference))
         end do
      end do
   end function max_mosum_statistic

   real(dp) function max_mosum_pvalue(statistic, nproc, width) result(p)
      real(dp), intent(in) :: statistic, width
      integer, intent(in) :: nproc

      p = efp_pvalue(statistic, "Brownian bridge increments", "max", &
         nproc, h = width)
   end function max_mosum_pvalue

   real(dp) function max_mosum_critical_value(alpha, nproc, width) result(value)
      real(dp), intent(in) :: alpha, width
      integer, intent(in) :: nproc
      real(dp) :: low, high, midpoint, p
      integer :: iteration

      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      low = 0.0_dp
      high = 100.0_dp
      do iteration = 1, 160
         midpoint = 0.5_dp * (low + high)
         p = max_mosum_pvalue(midpoint, nproc, width)
         if (p > alpha) then
            low = midpoint
         else
            high = midpoint
         end if
         if (high - low <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, midpoint)) exit
      end do
      value = 0.5_dp * (low + high)
   end function max_mosum_critical_value

   real(dp) function cat_l2_bb_statistic(process, frequencies) result(statistic)
      real(dp), intent(in) :: process(:, :), frequencies(:)
      real(dp), allocatable :: frequency(:)
      real(dp) :: cumulative, delta, previous
      integer :: first_row, i, index, j, n

      statistic = -huge(1.0_dp)
      if (size(process, 1) < 2 .or. size(process, 2) < 1) return
      if (size(frequencies) < 2 .or. any(frequencies <= 0.0_dp)) return
      allocate(frequency(size(frequencies)))
      frequency = frequencies / sum(frequencies)
      first_row = process_start_row(process)
      n = size(process, 1) - first_row + 1
      statistic = 0.0_dp
      do j = 1, size(process, 2)
         cumulative = 0.0_dp
         previous = 0.0_dp
         do i = 1, size(frequency)
            cumulative = cumulative + frequency(i)
            index = min(n, max(1, nint(cumulative * real(n, dp))))
            delta = process(first_row + index - 1, j) - previous
            statistic = statistic + delta * delta / frequency(i)
            previous = process(first_row + index - 1, j)
         end do
      end do
   end function cat_l2_bb_statistic

   real(dp) function cat_l2_bb_pvalue(statistic, nproc, n_categories) result(p)
      real(dp), intent(in) :: statistic
      integer, intent(in) :: nproc, n_categories
      real(dp) :: degrees_freedom

      degrees_freedom = real((n_categories - 1) * nproc, dp)
      if (degrees_freedom <= 0.0_dp) then
         p = 1.0_dp
      else
         p = r_pchisq(statistic, degrees_freedom, lower_tail = .false.)
      end if
   end function cat_l2_bb_pvalue

   real(dp) function cat_l2_bb_critical_value(alpha, nproc, n_categories) result(value)
      real(dp), intent(in) :: alpha
      integer, intent(in) :: nproc, n_categories
      real(dp) :: degrees_freedom

      degrees_freedom = real((n_categories - 1) * nproc, dp)
      if (degrees_freedom <= 0.0_dp .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         value = huge(1.0_dp)
      else
         value = r_qchisq(alpha, degrees_freedom, lower_tail = .false.)
      end if
   end function cat_l2_bb_critical_value

   integer function process_start_row(process) result(first_row)
      real(dp), intent(in) :: process(:, :)

      first_row = 1
      if (size(process, 1) > 1) then
         if (maxval(abs(process(1, :))) <= 16.0_dp * epsilon(1.0_dp)) first_row = 2
      end if
   end function process_start_row

   real(dp) function process_value(process, first_row, index, column) result(value)
      real(dp), intent(in) :: process(:, :)
      integer, intent(in) :: first_row, index, column

      if (index == 0) then
         value = 0.0_dp
      else
         value = process(first_row + index - 1, column)
      end if
   end function process_value
end module strucchange_functionals
