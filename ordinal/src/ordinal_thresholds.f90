! Based on ordinal/R/utils.R makeThresholds() and ordinal/R/clm.start.R.
! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_thresholds
   use ordinal_kinds, only : dp
   implicit none
   private
   integer, parameter, public :: threshold_flexible = 1
   integer, parameter, public :: threshold_symmetric = 2
   integer, parameter, public :: threshold_symmetric2 = 3
   integer, parameter, public :: threshold_equidistant = 4
   public :: threshold_parameter_count, threshold_jacobian, thresholds_from_alpha
   public :: threshold_start
contains
   pure integer function threshold_parameter_count(nclass, structure) result(nalpha)
      integer, intent(in) :: nclass !! Number of ordered response categories; must be at least two.
      integer, intent(in) :: structure !! Threshold structure identifier.
      integer :: ntheta
      ntheta = nclass - 1
      select case (structure)
      case (threshold_flexible)
         nalpha = ntheta
      case (threshold_symmetric)
         if (mod(ntheta, 2) == 1) then
            nalpha = (ntheta + 1)/2
         else
            nalpha = (ntheta + 2)/2
         end if
      case (threshold_symmetric2)
         if (mod(ntheta, 2) == 1) then
            nalpha = (ntheta - 1)/2
         else
            nalpha = ntheta/2
         end if
      case (threshold_equidistant)
         nalpha = 2
      case default
         nalpha = 0
      end select
   end function threshold_parameter_count

   pure subroutine threshold_jacobian(nclass, structure, tjac, status)
      integer, intent(in) :: nclass !! Number of ordered response categories; must be at least two.
      integer, intent(in) :: structure !! Threshold structure identifier.
      real(dp), allocatable, intent(out) :: tjac(:, :) !! Matrix mapping free alpha parameters to actual cut points.
      integer, intent(out) :: status !! Zero on success; nonzero for an invalid category count or structure.
      integer :: ntheta, nalpha, i, j, m
      status = 0
      ntheta = nclass - 1
      nalpha = threshold_parameter_count(nclass, structure)
      if (nclass < 2 .or. nalpha <= 0) then
         status = 1
         allocate(tjac(0, 0))
         return
      end if
      if (structure /= threshold_flexible .and. nclass < 3) then
         status = 2
         allocate(tjac(0, 0))
         return
      end if
      allocate(tjac(ntheta, nalpha))
      tjac = 0.0_dp
      select case (structure)
      case (threshold_flexible)
         do i = 1, ntheta
            tjac(i, i) = 1.0_dp
         end do
      case (threshold_equidistant)
         tjac(:, 1) = 1.0_dp
         do i = 1, ntheta
            tjac(i, 2) = real(i - 1, dp)
         end do
      case (threshold_symmetric2)
         m = nalpha
         if (mod(ntheta, 2) == 1) then
            do i = 1, m
               tjac(i, m - i + 1) = -1.0_dp
               tjac(m + 1 + i, i) = 1.0_dp
            end do
         else
            do i = 1, m
               tjac(i, m - i + 1) = -1.0_dp
               tjac(m + i, i) = 1.0_dp
            end do
         end if
      case (threshold_symmetric)
         if (mod(ntheta, 2) == 1) then
            m = nalpha
            tjac(:, 1) = 1.0_dp
            do j = 2, m
               do i = 1, m - j + 1
                  tjac(i, j) = -1.0_dp
               end do
               do i = m + j - 1, ntheta
                  tjac(i, j) = 1.0_dp
               end do
            end do
         else
            m = ntheta/2
            tjac(1:m, 1) = 1.0_dp
            tjac(m + 1:ntheta, 2) = 1.0_dp
            do j = 3, nalpha
               do i = 1, m - j + 2
                  tjac(i, j) = -1.0_dp
               end do
               do i = m + j - 1, ntheta
                  tjac(i, j) = 1.0_dp
               end do
            end do
         end if
      end select
   end subroutine threshold_jacobian

   pure subroutine thresholds_from_alpha(alpha, nclass, structure, theta, status)
      real(dp), intent(in) :: alpha(:) !! Free threshold parameters in the selected structure.
      integer, intent(in) :: nclass !! Number of ordered response categories.
      integer, intent(in) :: structure !! Threshold structure identifier.
      real(dp), allocatable, intent(out) :: theta(:) !! Actual ordered-category cut points before any nominal effects.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or structure.
      real(dp), allocatable :: tjac(:, :)
      call threshold_jacobian(nclass, structure, tjac, status)
      if (status /= 0) then
         allocate(theta(0))
         return
      end if
      if (size(alpha) /= size(tjac, 2)) then
         status = 3
         allocate(theta(0))
         return
      end if
      theta = matmul(tjac, alpha)
   end subroutine thresholds_from_alpha

   pure subroutine threshold_start(nclass, structure, alpha, status)
      integer, intent(in) :: nclass !! Number of ordered response categories.
      integer, intent(in) :: structure !! Threshold structure identifier.
      real(dp), allocatable, intent(out) :: alpha(:) !! Default starting threshold parameters matching ordinal's &
         !! logistic-quantile starts.
      integer, intent(out) :: status !! Zero on success; nonzero for an invalid threshold configuration.
      real(dp), allocatable :: theta(:)
      integer :: ntheta, nalpha, i, mid
      ntheta = nclass - 1
      nalpha = threshold_parameter_count(nclass, structure)
      if (nclass < 2 .or. nalpha <= 0) then
         status = 1
         allocate(alpha(0))
         return
      end if
      if (structure /= threshold_flexible .and. nclass < 3) then
         status = 2
         allocate(alpha(0))
         return
      end if
      allocate(theta(ntheta))
      do i = 1, ntheta
         theta(i) = log(real(i, dp)/real(ntheta + 1 - i, dp))
      end do
      allocate(alpha(nalpha))
      select case (structure)
      case (threshold_flexible)
         alpha = theta
      case (threshold_equidistant)
         alpha(1) = theta(1)
         alpha(2) = sum(theta(2:) - theta(:ntheta - 1))/real(ntheta - 1, dp)
      case (threshold_symmetric)
         if (mod(ntheta, 2) == 1) then
            mid = (ntheta + 1)/2
            alpha(1) = theta(mid)
            do i = 2, nalpha
               alpha(i) = theta(mid + i - 1) - theta(mid + i - 2)
            end do
         else
            mid = ntheta/2
            alpha(1) = theta(mid)
            alpha(2) = theta(mid + 1)
            do i = 3, nalpha
               alpha(i) = theta(mid + i - 1) - theta(mid + i - 2)
            end do
         end if
      case (threshold_symmetric2)
         if (mod(ntheta, 2) == 1) then
            mid = (ntheta + 1)/2
            do i = 1, nalpha
               alpha(i) = theta(mid + i)
            end do
         else
            mid = ntheta/2
            do i = 1, nalpha
               alpha(i) = theta(mid + i)
            end do
         end if
      end select
      status = 0
   end subroutine threshold_start
end module ordinal_thresholds
