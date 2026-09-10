! SPDX-License-Identifier: GPL-2.0-or-later
! Automatic component grouping translated from Rssa 1.1.
module rssa_autogroup
   use rssa_kinds, only : dp
   use rssa_metrics, only : wcor_ssa
   use rssa_reconstruction, only : elementary_series_ssa
   use rssa_types, only : grouping_result, rssa_invalid_input, rssa_success, ssa_result
   implicit none
   private
   public :: grouping_auto, grouping_auto_pgram_ssa, grouping_auto_wcor_ssa
contains
   function grouping_auto_wcor_ssa(object, indices, threshold) result(out)
      type(ssa_result), intent(in) :: object !! Decomposed SSA object whose elementary series are clustered.
      integer, intent(in), optional :: indices(:) !! Components to cluster; defaults to all stored components.
      real(dp), intent(in), optional :: threshold !! Absolute w-correlation merge threshold; default 0.5.
      type(grouping_result) :: out
      integer, allocatable :: idx(:), label(:)
      real(dp), allocatable :: cor(:, :)
      real(dp) :: cutoff
      integer :: i, j, n, next_label
      if (.not. allocated(object%sigma)) then
         out%info = rssa_invalid_input
         allocate(out%group(0), out%score(0))
         return
      end if
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i=1, size(idx))]
      end if
      n = size(idx)
      cutoff = 0.5_dp
      if (present(threshold)) cutoff = threshold
      cor = wcor_ssa(object, idx)
      allocate(label(n), out%score(n))
      label = 0
      out%score = 0.0_dp
      next_label = 0
      do i = 1, n
         if (label(i) /= 0) cycle
         next_label = next_label + 1
         label(i) = next_label
         do j = i + 1, n
            if (abs(cor(i, j)) >= cutoff) label(j) = next_label
         end do
         if (n > 1) out%score(i) = maxval(abs(cor(i, :)), mask=[(.true., j=1,n)])
      end do
      out%group = label
      out%n_groups = next_label
      out%info = rssa_success
   end function grouping_auto_wcor_ssa

   function grouping_auto_pgram_ssa(object, indices, n_groups) result(out)
      type(ssa_result), intent(in) :: object !! Decomposed SSA object whose components are grouped by dominant frequency.
      integer, intent(in), optional :: indices(:) !! Components to group; defaults to all stored components.
      integer, intent(in), optional :: n_groups !! Number of frequency bands; default two.
      type(grouping_result) :: out
      integer, allocatable :: idx(:)
      real(dp), allocatable :: component(:)
      real(dp) :: pi, angle, re, im, power, best_power, freq
      integer :: i, j, k, n, ng, best_k
      pi = acos(-1.0_dp)
      if (.not. allocated(object%sigma)) then
         out%info = rssa_invalid_input
         allocate(out%group(0), out%score(0))
         return
      end if
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i=1, size(idx))]
      end if
      ng = 2
      if (present(n_groups)) ng = max(1, n_groups)
      allocate(out%group(size(idx)), out%score(size(idx)))
      do i = 1, size(idx)
         component = elementary_series_ssa(object, idx(i))
         n = size(component)
         best_power = -1.0_dp
         best_k = 0
         do k = 0, n / 2
            re = 0.0_dp
            im = 0.0_dp
            do j = 1, n
               angle = 2.0_dp * pi * real(k * (j - 1), dp) / real(max(1, n), dp)
               re = re + component(j) * cos(angle)
               im = im - component(j) * sin(angle)
            end do
            power = re * re + im * im
            if (power > best_power) then
               best_power = power
               best_k = k
            end if
         end do
         freq = real(best_k, dp) / real(max(1, n), dp)
         out%score(i) = freq
         out%group(i) = min(ng, 1 + int(2.0_dp * freq * real(ng, dp)))
      end do
      out%n_groups = ng
      out%info = rssa_success
   end function grouping_auto_pgram_ssa

   function grouping_auto(object, method, indices, n_groups, threshold) result(out)
      type(ssa_result), intent(in) :: object !! SSA object to group.
      character(len=*), intent(in), optional :: method !! `wcor` or `pgram`; default `wcor`.
      integer, intent(in), optional :: indices(:) !! Components to group.
      integer, intent(in), optional :: n_groups !! Number of spectral groups for `pgram`.
      real(dp), intent(in), optional :: threshold !! W-correlation threshold for `wcor`.
      type(grouping_result) :: out
      character(len=16) :: selected
      selected = 'wcor'
      if (present(method)) selected = trim(adjustl(method))
      if (selected == 'pgram') then
         out = grouping_auto_pgram_ssa(object, indices, n_groups)
      else
         out = grouping_auto_wcor_ssa(object, indices, threshold)
      end if
   end function grouping_auto
end module rssa_autogroup
