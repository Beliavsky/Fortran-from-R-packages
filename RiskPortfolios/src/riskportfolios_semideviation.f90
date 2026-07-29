! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_semideviation
   use riskportfolios_kinds, only : dp
   use riskportfolios_stats, only : column_means
   implicit none
   private

   integer, parameter, public :: SEMIDEV_NAIVE = 1
   integer, parameter, public :: SEMIDEV_EWMA = 2
   public :: semideviation_estimation, naive_semideviation, ewma_semideviation

contains

   subroutine semideviation_estimation(rets, semidev, method, lambda, info)
      real(dp), intent(in) :: rets(:, :)
      real(dp), allocatable, intent(out) :: semidev(:)
      integer, intent(in), optional :: method
      real(dp), intent(in), optional :: lambda
      integer, intent(out), optional :: info
      integer :: m
      real(dp) :: lam

      m = SEMIDEV_NAIVE
      if (present(method)) m = method
      lam = 0.94_dp
      if (present(lambda)) lam = lambda
      allocate(semidev(size(rets, 2)))
      select case (m)
      case (SEMIDEV_NAIVE)
         semidev = naive_semideviation(rets)
         if (present(info)) info = 0
      case (SEMIDEV_EWMA)
         semidev = ewma_semideviation(rets, lam)
         if (present(info)) info = 0
      case default
         semidev = 0.0_dp
         if (present(info)) info = -1
      end select
   end subroutine semideviation_estimation

   pure function naive_semideviation(rets) result(semidev)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: semidev(size(rets, 2))
      semidev = ewma_semideviation(rets, 1.0_dp)
   end function naive_semideviation

   pure function ewma_semideviation(rets, lambda) result(semidev)
      real(dp), intent(in) :: rets(:, :)
      real(dp), intent(in) :: lambda
      real(dp) :: semidev(size(rets, 2))
      real(dp) :: mu(size(rets, 2)), weights(size(rets, 1))
      real(dp) :: denominator, d
      integer :: i, j, t

      t = size(rets, 1)
      mu = column_means(rets)
      semidev = 0.0_dp
      if (t == 0) return
      do i = 1, t
         weights(i) = lambda ** real(t - i + 1, dp)
      end do
      do j = 1, size(rets, 2)
         denominator = 0.0_dp
         d = 0.0_dp
         do i = 1, t
            if (rets(i, j) < mu(j)) then
               denominator = denominator + weights(i)
               d = d + weights(i) * (rets(i, j) - mu(j)) ** 2
            end if
         end do
         if (denominator > tiny(1.0_dp)) semidev(j) = sqrt(d / denominator)
      end do
   end function ewma_semideviation

end module riskportfolios_semideviation
