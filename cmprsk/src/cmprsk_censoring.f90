! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_censoring
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use cmprsk_status, only : cmprsk_success, cmprsk_invalid_argument
   implicit none
   private

   public :: censoring_survival_left

contains

   pure subroutine censoring_survival_left(ftime, censored, group, survival_left, status)
      real(dp), intent(in) :: ftime(:) !! Event/censoring times, sorted in nondecreasing order.
      logical, intent(in) :: censored(:) !! True for censoring observations and false otherwise.
      integer, intent(in) :: group(:) !! Censoring-group codes in `1:maxval(group)`.
      real(dp), intent(out) :: survival_left(:, :) !! Group-specific Kaplan-Meier censoring survival at each `ftime(i)-`.
      integer, intent(out) :: status !! `cmprsk_success` or `cmprsk_invalid_argument`.

      integer :: g
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: ng
      integer :: nrisk
      integer :: ncens
      real(dp) :: s
      real(dp) :: time_value

      n = size(ftime)
      if (size(censored) /= n .or. size(group) /= n) then
         status = cmprsk_invalid_argument
         return
      end if
      if (n == 0) then
         if (size(survival_left, 2) /= 0) then
            status = cmprsk_invalid_argument
         else
            survival_left = 1.0_dp
            status = cmprsk_success
         end if
         return
      end if
      if (any(ieee_is_nan(ftime)) .or. any(group < 1)) then
         status = cmprsk_invalid_argument
         return
      end if
      ng = maxval(group)
      if (size(survival_left, 1) /= ng .or. size(survival_left, 2) /= n) then
         status = cmprsk_invalid_argument
         return
      end if
      do i = 2, n
         if (ftime(i) < ftime(i - 1)) then
            status = cmprsk_invalid_argument
            return
         end if
      end do

      survival_left = 1.0_dp
      do g = 1, ng
         nrisk = count(group == g)
         s = 1.0_dp
         i = 1
         do while (i <= n)
            time_value = ftime(i)
            j = i
            do while (j <= n)
               if (ftime(j) /= time_value) exit
               j = j + 1
            end do
            do k = i, j - 1
               if (group(k) == g) survival_left(g, k) = s
            end do
            ncens = 0
            do k = i, j - 1
               if (group(k) == g .and. censored(k)) ncens = ncens + 1
            end do
            if (nrisk > 0 .and. ncens > 0) s = s * real(nrisk - ncens, dp) / real(nrisk, dp)
            do k = i, j - 1
               if (group(k) == g) nrisk = nrisk - 1
            end do
            i = j
         end do
      end do
      status = cmprsk_success
   end subroutine censoring_survival_left

end module cmprsk_censoring
