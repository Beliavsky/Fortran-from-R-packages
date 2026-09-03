! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_api
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use cmprsk_status, only : cmprsk_success, cmprsk_invalid_argument, cmprsk_no_failure_of_interest
   use cmprsk_utils, only : stable_order_real
   use cmprsk_cuminc, only : cuminc_curve, gray_test_result, cumulative_incidence, gray_test, curve_timepoints
   implicit none
   private

   type, public :: cuminc_entry
      integer :: group_label = 0
      integer :: cause_label = 0
      type(cuminc_curve) :: curve
   end type cuminc_entry

   type, public :: cuminc_result
      type(cuminc_entry), allocatable :: estimates(:)
      integer, allocatable :: tested_cause(:)
      type(gray_test_result), allocatable :: tests(:)
      real(dp) :: rho = 0.0_dp
   end type cuminc_result

   public :: fit_cuminc
   public :: cuminc_timepoints

contains

   pure subroutine fit_cuminc(ftime, fstatus, result, status, group, strata, cencode, rho)
      real(dp), intent(in) :: ftime(:) !! Failure/censoring times in arbitrary row order.
      integer, intent(in) :: fstatus(:) !! Cause code for each row, with `cencode` denoting censoring.
      type(cuminc_result), intent(out) :: result !! Group-by-cause cumulative-incidence curves and optional Gray tests.
      integer, intent(out) :: status !! Success or a documented `cmprsk_status` code.
      integer, intent(in), optional :: group(:) !! Group labels; defaults to one group and may use arbitrary integers.
      integer, intent(in), optional :: strata(:) !! Gray-test strata; defaults to one stratum and accepts arbitrary integer labels.
      integer, intent(in), optional :: cencode !! Censoring code; defaults to `0`.
      real(dp), intent(in), optional :: rho !! Gray-test weight exponent; defaults to `0`.

      integer :: cencode_value
      integer :: cidx
      integer :: eidx
      integer :: gidx
      integer :: i
      integer :: n
      integer :: ng
      integer :: ns
      integer :: ncauses
      integer, allocatable :: causes(:)
      integer, allocatable :: failed_any(:)
      integer, allocatable :: failed_cause(:)
      integer, allocatable :: group_code(:)
      integer, allocatable :: group_label(:)
      integer, allocatable :: group_raw(:)
      integer, allocatable :: order(:)
      integer, allocatable :: status_gray(:)
      integer, allocatable :: status_sorted(:)
      integer, allocatable :: strata_code(:)
      integer, allocatable :: strata_raw(:)
      integer, allocatable :: subset_index(:)
      real(dp) :: rho_value
      real(dp), allocatable :: time_subset(:)
      real(dp), allocatable :: time_sorted(:)

      n = size(ftime)
      cencode_value = 0
      if (present(cencode)) cencode_value = cencode
      rho_value = 0.0_dp
      if (present(rho)) rho_value = rho
      if (n < 1 .or. size(fstatus) /= n .or. any(ieee_is_nan(ftime)) .or. any(ftime < 0.0_dp) .or. &
          ieee_is_nan(rho_value)) then
         status = cmprsk_invalid_argument
         return
      end if
      if (present(group)) then
         if (size(group) /= n) then
            status = cmprsk_invalid_argument
            return
         end if
      end if
      if (present(strata)) then
         if (size(strata) /= n) then
            status = cmprsk_invalid_argument
            return
         end if
      end if

      allocate(order(n), time_sorted(n), status_sorted(n), group_raw(n), strata_raw(n))
      call stable_order_real(ftime, order)
      do i = 1, n
         time_sorted(i) = ftime(order(i))
         status_sorted(i) = fstatus(order(i))
         if (present(group)) then
            group_raw(i) = group(order(i))
         else
            group_raw(i) = 1
         end if
         if (present(strata)) then
            strata_raw(i) = strata(order(i))
         else
            strata_raw(i) = 1
         end if
      end do
      call integer_levels(group_raw, group_label)
      ng = size(group_label)
      allocate(group_code(n), strata_code(n))
      call map_to_levels(group_raw, group_label, group_code)
      call integer_levels(strata_raw, causes)
      call map_to_levels(strata_raw, causes, strata_code)
      deallocate(causes)

      call noncensor_causes(status_sorted, cencode_value, causes)
      ncauses = size(causes)
      if (ncauses == 0) then
         status = cmprsk_no_failure_of_interest
         return
      end if
      allocate(result%estimates(ng*ncauses))
      result%rho = rho_value
      eidx = 0
      do cidx = 1, ncauses
         do gidx = 1, ng
            ns = count(group_code == gidx)
            allocate(subset_index(ns), time_subset(ns), failed_any(ns), failed_cause(ns))
            ns = 0
            do i = 1, n
               if (group_code(i) /= gidx) cycle
               ns = ns + 1
               subset_index(ns) = i
               time_subset(ns) = time_sorted(i)
               failed_any(ns) = merge(0, 1, status_sorted(i) == cencode_value)
               failed_cause(ns) = merge(1, 0, status_sorted(i) == causes(cidx))
            end do
            eidx = eidx + 1
            result%estimates(eidx)%group_label = group_label(gidx)
            result%estimates(eidx)%cause_label = causes(cidx)
            call cumulative_incidence(time_subset, failed_any, failed_cause, result%estimates(eidx)%curve, status)
            if (status /= cmprsk_success) return
            deallocate(subset_index, time_subset, failed_any, failed_cause)
         end do
      end do

      if (ng > 1) then
         allocate(result%tested_cause(ncauses), result%tests(ncauses), status_gray(n))
         result%tested_cause = causes
         do cidx = 1, ncauses
            do i = 1, n
               if (status_sorted(i) == cencode_value) then
                  status_gray(i) = 0
               else if (status_sorted(i) == causes(cidx)) then
                  status_gray(i) = 1
               else
                  status_gray(i) = 2
               end if
            end do
            call gray_test(time_sorted, status_gray, group_code, strata_code, rho_value, result%tests(cidx), status)
            if (status /= cmprsk_success) return
         end do
      else
         allocate(result%tested_cause(0), result%tests(0))
      end if
      status = cmprsk_success
   end subroutine fit_cuminc

   pure subroutine cuminc_timepoints(result, requested_time, estimate, variance, present_value)
      type(cuminc_result), intent(in) :: result !! Cumulative-incidence collection returned by `fit_cuminc`.
      real(dp), intent(in) :: requested_time(:) !! Requested evaluation times, preferably in ascending order.
      real(dp), allocatable, intent(out) :: estimate(:, :) !! Curve estimates by collection entry and requested time.
      real(dp), allocatable, intent(out) :: variance(:, :) !! Curve variances by collection entry and requested time.
      logical, allocatable, intent(out) :: present_value(:, :) !! True where a curve extends to the requested time.

      integer :: i
      integer :: n_entries
      integer :: n_times

      n_entries = size(result%estimates)
      n_times = size(requested_time)
      allocate(estimate(n_entries, n_times), variance(n_entries, n_times), present_value(n_entries, n_times))
      do i = 1, n_entries
         call curve_timepoints(result%estimates(i)%curve, requested_time, estimate(i, :), variance(i, :), present_value(i, :))
      end do
   end subroutine cuminc_timepoints

   pure subroutine integer_levels(values, levels)
      integer, intent(in) :: values(:) !! Integer labels from which sorted distinct levels are requested.
      integer, allocatable, intent(out) :: levels(:) !! Sorted distinct integer labels.

      integer :: i
      integer :: j
      integer :: count
      integer :: key
      integer, allocatable :: work(:)

      if (size(values) == 0) then
         allocate(levels(0))
         return
      end if
      allocate(work(size(values)))
      work = values
      do i = 2, size(work)
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      count = 1
      do i = 2, size(work)
         if (work(i) /= work(count)) then
            count = count + 1
            work(count) = work(i)
         end if
      end do
      allocate(levels(count))
      levels = work(1:count)
   end subroutine integer_levels

   pure subroutine map_to_levels(values, levels, code)
      integer, intent(in) :: values(:) !! Integer labels to map.
      integer, intent(in) :: levels(:) !! Distinct labels defining consecutive code order.
      integer, intent(out) :: code(:) !! One-based index in `levels` for each input label.

      integer :: i
      integer :: j

      code = 0
      do i = 1, size(values)
         do j = 1, size(levels)
            if (values(i) == levels(j)) then
               code(i) = j
               exit
            end if
         end do
      end do
   end subroutine map_to_levels

   pure subroutine noncensor_causes(status_code, cencode, causes)
      integer, intent(in) :: status_code(:) !! Observed cause codes.
      integer, intent(in) :: cencode !! Cause code denoting censoring and excluded from the result.
      integer, allocatable, intent(out) :: causes(:) !! Sorted distinct noncensoring cause codes.

      integer :: i
      integer :: k
      integer, allocatable :: work(:)
      integer, allocatable :: levels(:)

      allocate(work(count(status_code /= cencode)))
      k = 0
      do i = 1, size(status_code)
         if (status_code(i) == cencode) cycle
         k = k + 1
         work(k) = status_code(i)
      end do
      call integer_levels(work, levels)
      allocate(causes(size(levels)))
      causes = levels
   end subroutine noncensor_causes

end module cmprsk_api
