! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_mean
   use riskportfolios_kinds, only : dp
   use riskportfolios_stats, only : column_means, sample_covariance, standard_deviations
   use riskportfolios_linalg, only : solve_linear
   implicit none
   private

   integer, parameter, public :: MEAN_NAIVE = 1
   integer, parameter, public :: MEAN_EWMA = 2
   integer, parameter, public :: MEAN_BAYES_STEIN = 3
   integer, parameter, public :: MEAN_MARTELLINI = 4
   public :: mean_estimation, naive_mean, ewma_mean, bayes_stein_mean, martellini_mean

contains

   subroutine mean_estimation(rets, mu, method, lambda, info)
      real(dp), intent(in) :: rets(:, :)
      real(dp), allocatable, intent(out) :: mu(:)
      integer, intent(in), optional :: method
      real(dp), intent(in), optional :: lambda
      integer, intent(out), optional :: info
      integer :: m, stat
      real(dp) :: lam

      m = MEAN_NAIVE
      if (present(method)) m = method
      lam = 0.94_dp
      if (present(lambda)) lam = lambda
      stat = 0
      select case (m)
      case (MEAN_NAIVE)
         mu = naive_mean(rets)
      case (MEAN_EWMA)
         mu = ewma_mean(rets, lam)
      case (MEAN_BAYES_STEIN)
         call bayes_stein_mean(rets, mu, stat)
      case (MEAN_MARTELLINI)
         mu = martellini_mean(rets)
      case default
         allocate(mu(0))
         stat = -1
      end select
      if (present(info)) info = stat
   end subroutine mean_estimation

   pure function naive_mean(rets) result(mu)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: mu(size(rets, 2))
      mu = column_means(rets)
   end function naive_mean

   pure function ewma_mean(rets, lambda) result(mu)
      real(dp), intent(in) :: rets(:, :)
      real(dp), intent(in) :: lambda
      real(dp) :: mu(size(rets, 2))
      real(dp) :: weights(size(rets, 1)), denom
      integer :: i, t

      t = size(rets, 1)
      mu = 0.0_dp
      if (t == 0) return
      do i = 1, t
         weights(i) = lambda ** real(t - i + 1, dp)
      end do
      denom = sum(weights)
      if (denom <= tiny(1.0_dp)) return
      weights = weights / denom
      mu = matmul(transpose(rets), weights)
   end function ewma_mean

   subroutine bayes_stein_mean(rets, mu_bs, info)
      real(dp), intent(in) :: rets(:, :)
      real(dp), allocatable, intent(out) :: mu_bs(:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: inv_sigma_one(:), inv_sigma_dev(:)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp) :: mu(size(rets, 2)), one(size(rets, 2))
      real(dp) :: w_min(size(rets, 2)), mu_min, phi, denom
      integer :: n, t, stat

      n = size(rets, 2)
      t = size(rets, 1)
      allocate(mu_bs(n))
      if (n == 0 .or. t <= 1) then
         mu_bs = 0.0_dp
         if (present(info)) info = -1
         return
      end if
      mu = column_means(rets)
      sigma = sample_covariance(rets)
      one = 1.0_dp
      call solve_linear(sigma, one, inv_sigma_one, stat)
      if (stat /= 0) then
         mu_bs = mu
         if (present(info)) info = stat
         return
      end if
      denom = dot_product(one, inv_sigma_one)
      if (abs(denom) <= tiny(1.0_dp)) then
         mu_bs = mu
         if (present(info)) info = -2
         return
      end if
      w_min = inv_sigma_one / denom
      mu_min = dot_product(mu, w_min)
      call solve_linear(sigma, mu - mu_min, inv_sigma_dev, stat)
      if (stat /= 0) then
         mu_bs = mu
         if (present(info)) info = stat
         return
      end if
      phi = real(n + 2, dp) / (real(n + 2, dp) + real(t, dp) * &
         dot_product(mu - mu_min, inv_sigma_dev))
      phi = max(0.0_dp, min(1.0_dp, phi))
      mu_bs = (1.0_dp - phi) * mu + phi * mu_min
      if (present(info)) info = 0
   end subroutine bayes_stein_mean

   pure function martellini_mean(rets) result(mu)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: mu(size(rets, 2))
      mu = standard_deviations(rets)
   end function martellini_mean

end module riskportfolios_mean
