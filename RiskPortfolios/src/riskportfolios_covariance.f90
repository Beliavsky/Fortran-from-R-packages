! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_covariance
   use riskportfolios_kinds, only : dp
   use riskportfolios_stats, only : column_means, sample_covariance, &
      population_covariance, standard_deviations, covariance_to_correlation, &
      frobenius_norm, symmetrize
   use riskportfolios_linalg, only : solve_linear, symmetric_eigen
   implicit none
   private

   integer, parameter, public :: COV_NAIVE = 1
   integer, parameter, public :: COV_EWMA = 2
   integer, parameter, public :: COV_LEDOIT_WOLF = 3
   integer, parameter, public :: COV_FACTOR = 4
   integer, parameter, public :: COV_CONSTANT = 5
   integer, parameter, public :: COV_COR_SHRINKAGE = 6
   integer, parameter, public :: COV_ONE_PARAMETER = 7
   integer, parameter, public :: COV_DIAGONAL = 8
   integer, parameter, public :: COV_LARGE = 9
   integer, parameter, public :: COV_BAYES_STEIN = 10

   public :: covariance_estimation
   public :: naive_covariance, ewma_covariance, constant_covariance
   public :: factor_covariance, ledoit_wolf_covariance, large_covariance
   public :: correlation_shrinkage_covariance, diagonal_shrinkage_covariance
   public :: one_parameter_covariance, bayes_stein_covariance

contains

   subroutine covariance_estimation(rets, sigma, method, lambda, n_factors, info)
      real(dp), intent(in) :: rets(:, :)
      real(dp), allocatable, intent(out) :: sigma(:, :)
      integer, intent(in), optional :: method, n_factors
      real(dp), intent(in), optional :: lambda
      integer, intent(out), optional :: info
      integer :: m, k, stat, n
      real(dp) :: lam

      n = size(rets, 2)
      allocate(sigma(n, n))
      m = COV_NAIVE
      if (present(method)) m = method
      lam = 0.94_dp
      if (present(lambda)) lam = lambda
      k = 1
      if (present(n_factors)) k = n_factors
      stat = 0

      select case (m)
      case (COV_NAIVE)
         sigma = naive_covariance(rets)
      case (COV_EWMA)
         sigma = ewma_covariance(rets, lam)
      case (COV_LEDOIT_WOLF)
         sigma = ledoit_wolf_covariance(rets)
      case (COV_FACTOR)
         call factor_covariance(rets, k, sigma, stat)
      case (COV_CONSTANT)
         sigma = constant_covariance(rets)
      case (COV_COR_SHRINKAGE)
         sigma = correlation_shrinkage_covariance(rets)
      case (COV_ONE_PARAMETER)
         sigma = one_parameter_covariance(rets)
      case (COV_DIAGONAL)
         sigma = diagonal_shrinkage_covariance(rets)
      case (COV_LARGE)
         sigma = large_covariance(rets)
      case (COV_BAYES_STEIN)
         call bayes_stein_covariance(rets, sigma, stat)
      case default
         sigma = 0.0_dp
         stat = -1
      end select
      sigma = symmetrize(sigma)
      if (present(info)) info = stat
   end subroutine covariance_estimation

   pure function naive_covariance(rets) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      sigma = sample_covariance(rets)
   end function naive_covariance

   pure function ewma_covariance(rets, lambda) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp), intent(in) :: lambda
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp) :: mu(size(rets, 2)), r(size(rets, 2))
      real(dp) :: alpha, denominator
      integer :: i, t

      t = size(rets, 1)
      sigma = sample_covariance(rets)
      if (t == 0) return
      mu = column_means(rets)
      denominator = 1.0_dp - lambda ** real(t, dp)
      if (abs(denominator) <= tiny(1.0_dp)) then
         alpha = 1.0_dp / real(t, dp)
      else
         alpha = (1.0_dp - lambda) / denominator
      end if
      do i = 1, t
         r = rets(i, :) - mu
         sigma = alpha * outer_product(r, r) + lambda * sigma
      end do
      sigma = symmetrize(sigma)
   end function ewma_covariance

   pure function constant_covariance(rets) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp) :: rho(size(rets, 2), size(rets, 2))
      real(dp) :: sd(size(rets, 2)), rbar
      integer :: i, j, n, count

      n = size(rets, 2)
      sigma = sample_covariance(rets)
      rho = covariance_to_correlation(sigma)
      sd = standard_deviations(rets)
      rbar = 0.0_dp
      count = 0
      do j = 1, n
         do i = j + 1, n
            rbar = rbar + rho(i, j)
            count = count + 1
         end do
      end do
      if (count > 0) rbar = rbar / real(count, dp)
      do j = 1, n
         do i = 1, n
            if (i == j) then
               sigma(i, j) = sd(i) ** 2
            else
               sigma(i, j) = rbar * sd(i) * sd(j)
            end if
         end do
      end do
   end function constant_covariance

   subroutine factor_covariance(rets, n_factors, sigma, info)
      real(dp), intent(in) :: rets(:, :)
      integer, intent(in) :: n_factors
      real(dp), intent(out) :: sigma(:, :)
      integer, intent(out), optional :: info
      real(dp), allocatable :: values(:), vectors(:, :), loadings(:, :)
      real(dp) :: sample(size(rets, 2), size(rets, 2))
      real(dp) :: rho(size(rets, 2), size(rets, 2))
      real(dp) :: rebuilt(size(rets, 2), size(rets, 2))
      real(dp) :: sd(size(rets, 2)), uniqueness(size(rets, 2))
      integer :: i, j, f, idx, k, n, stat

      n = size(rets, 2)
      if (n == 0 .or. n_factors < 1) then
         sigma = 0.0_dp
         if (present(info)) info = -1
         return
      end if
      k = min(n_factors, max(1, n - 1))
      sample = sample_covariance(rets)
      rho = covariance_to_correlation(sample)
      sd = standard_deviations(rets)
      call symmetric_eigen(rho, values, vectors, stat)
      if (stat /= 0) then
         sigma = sample
         if (present(info)) info = stat
         return
      end if
      allocate(loadings(n, k))
      do f = 1, k
         idx = n - f + 1
         loadings(:, f) = vectors(:, idx) * sqrt(max(values(idx), 0.0_dp))
      end do
      uniqueness = 1.0_dp
      do i = 1, n
         uniqueness(i) = max(1.0e-8_dp, 1.0_dp - sum(loadings(i, :) ** 2))
      end do
      rebuilt = matmul(loadings, transpose(loadings))
      do i = 1, n
         rebuilt(i, i) = rebuilt(i, i) + uniqueness(i)
      end do
      do j = 1, n
         do i = 1, n
            sigma(i, j) = sd(i) * rebuilt(i, j) * sd(j)
         end do
      end do
      sigma = symmetrize(sigma)
      if (present(info)) info = 0
   end subroutine factor_covariance

   function ledoit_wolf_covariance(rets) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp), allocatable :: x(:, :), y(:, :), z(:, :)
      real(dp), allocatable :: sample(:, :), prior(:, :), phi_mat(:, :)
      real(dp), allocatable :: rho1(:, :), rho3(:, :), rho_mat(:, :)
      real(dp), allocatable :: mkt(:), cov_mkt(:)
      real(dp) :: var_mkt, phi, rho_value, gamma, kappa, shrinkage
      integer :: i, j, t, n

      t = size(rets, 1)
      n = size(rets, 2)
      if (t == 0 .or. n == 0) then
         sigma = 0.0_dp
         return
      end if
      allocate(x(t, n), y(t, n), z(t, n), sample(n, n), prior(n, n))
      allocate(phi_mat(n, n), rho1(n, n), rho3(n, n), rho_mat(n, n))
      allocate(mkt(t), cov_mkt(n))
      x = rets - spread(column_means(rets), 1, t)
      y = x * x
      mkt = sum(x, dim=2) / real(n, dp)
      sample = matmul(transpose(x), x) / real(t, dp)
      cov_mkt = matmul(transpose(x), mkt) / real(t, dp)
      var_mkt = dot_product(mkt, mkt) / real(t, dp)
      if (var_mkt <= tiny(1.0_dp)) then
         sigma = sample
         return
      end if
      prior = outer_product(cov_mkt, cov_mkt) / var_mkt
      do i = 1, n
         prior(i, i) = sample(i, i)
      end do
      do j = 1, n
         z(:, j) = x(:, j) * mkt
      end do
      phi_mat = matmul(transpose(y), y) / real(t, dp) - &
         2.0_dp * matmul(transpose(x), x) * sample / real(t, dp) + sample * sample
      phi = sum(phi_mat)
      rho1 = matmul(transpose(y), z) / real(t, dp)
      do j = 1, n
         rho1(:, j) = rho1(:, j) * cov_mkt(j) / var_mkt
      end do
      rho3 = matmul(transpose(z), z) / real(t, dp) * &
         outer_product(cov_mkt, cov_mkt) / (var_mkt * var_mkt)
      rho_mat = 2.0_dp * rho1 - rho3 - prior * sample
      do i = 1, n
         rho_mat(i, i) = phi_mat(i, i)
      end do
      rho_value = sum(rho_mat)
      gamma = frobenius_norm(sample - prior) ** 2
      if (gamma <= tiny(1.0_dp)) then
         shrinkage = 1.0_dp
      else
         kappa = (phi - rho_value) / gamma
         shrinkage = max(0.0_dp, min(1.0_dp, kappa / real(t, dp)))
      end if
      sigma = shrinkage * prior + (1.0_dp - shrinkage) * sample
      sigma = symmetrize(sigma)
   end function ledoit_wolf_covariance

   function large_covariance(rets) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp), allocatable :: x(:, :), y(:, :), z(:, :), mkt(:), cov_mkt(:)
      real(dp), allocatable :: sample(:, :), prior(:, :), v1(:, :), v3(:, :)
      real(dp) :: var_mkt, d, r2, phi_diag, phi_off1, phi_off3, phi_off
      real(dp) :: phi, shrinkage, temp
      integer :: i, j, t, n

      t = size(rets, 1)
      n = size(rets, 2)
      if (t == 0 .or. n == 0) then
         sigma = 0.0_dp
         return
      end if
      allocate(x(t, n), y(t, n), z(t, n), mkt(t), cov_mkt(n))
      allocate(sample(n, n), prior(n, n), v1(n, n), v3(n, n))
      x = rets - spread(column_means(rets), 1, t)
      y = x * x
      mkt = sum(x, dim=2) / real(n, dp)
      sample = matmul(transpose(x), x) / real(t, dp)
      cov_mkt = matmul(transpose(x), mkt) / real(t, dp)
      var_mkt = dot_product(mkt, mkt) / real(t, dp)
      if (var_mkt <= tiny(1.0_dp)) then
         sigma = sample
         return
      end if
      prior = outer_product(cov_mkt, cov_mkt) / var_mkt
      do i = 1, n
         prior(i, i) = sample(i, i)
      end do
      do j = 1, n
         z(:, j) = x(:, j) * mkt
      end do
      d = frobenius_norm(sample - prior) ** 2 / real(n, dp)
      r2 = sum(matmul(transpose(y), y)) / &
         (real(n, dp) * real(t, dp) ** 2) - &
         sum(sample * sample) / (real(n, dp) * real(t, dp))
      phi_diag = sum(y * y) / (real(n, dp) * real(t, dp) ** 2)
      do i = 1, n
         phi_diag = phi_diag - sample(i, i) ** 2 / &
            (real(n, dp) * real(t, dp))
      end do
      v1 = matmul(transpose(y), z) / real(t, dp) ** 2
      do j = 1, n
         do i = 1, n
            v1(i, j) = v1(i, j) - cov_mkt(j) * sample(i, j) / real(t, dp)
         end do
      end do
      temp = 0.0_dp
      do j = 1, n
         do i = 1, n
            temp = temp + v1(i, j) * cov_mkt(i)
         end do
      end do
      phi_off1 = temp / (real(n, dp) * var_mkt)
      do i = 1, n
         phi_off1 = phi_off1 - v1(i, i) * cov_mkt(i) / &
            (real(n, dp) * var_mkt)
      end do
      v3 = matmul(transpose(z), z) / real(t, dp) ** 2 - &
         sample * var_mkt / real(t, dp)
      temp = 0.0_dp
      do j = 1, n
         do i = 1, n
            temp = temp + v3(i, j) * cov_mkt(i) * cov_mkt(j)
         end do
      end do
      phi_off3 = temp / (real(n, dp) * var_mkt ** 2)
      do i = 1, n
         phi_off3 = phi_off3 - v3(i, i) * cov_mkt(i) ** 2 / &
            (real(n, dp) * var_mkt ** 2)
      end do
      phi_off = 2.0_dp * phi_off1 - phi_off3
      phi = phi_diag + phi_off
      if (d <= tiny(1.0_dp)) then
         shrinkage = 1.0_dp
      else
         shrinkage = max(0.0_dp, min(1.0_dp, (r2 - phi) / d))
      end if
      sigma = shrinkage * prior + (1.0_dp - shrinkage) * sample
      sigma = symmetrize(sigma)
   end function large_covariance

   function correlation_shrinkage_covariance(rets) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp), allocatable :: x(:, :), y(:, :), sample(:, :), prior(:, :)
      real(dp), allocatable :: phi_mat(:, :), rho_mat(:, :), term1(:, :)
      real(dp), allocatable :: variance(:), sqrt_variance(:)
      real(dp) :: rbar, phi, rho_value, gamma, kappa, shrinkage
      integer :: i, j, t, n

      t = size(rets, 1)
      n = size(rets, 2)
      if (t == 0 .or. n == 0) then
         sigma = 0.0_dp
         return
      end if
      allocate(x(t, n), y(t, n), sample(n, n), prior(n, n))
      allocate(phi_mat(n, n), rho_mat(n, n), term1(n, n))
      allocate(variance(n), sqrt_variance(n))
      x = rets - spread(column_means(rets), 1, t)
      y = x * x
      sample = matmul(transpose(x), x) / real(t, dp)
      do i = 1, n
         variance(i) = max(sample(i, i), tiny(1.0_dp))
      end do
      sqrt_variance = sqrt(variance)
      rbar = 0.0_dp
      if (n > 1) then
         do j = 1, n
            do i = 1, n
               rbar = rbar + sample(i, j) / &
                  (sqrt_variance(i) * sqrt_variance(j))
            end do
         end do
         rbar = (rbar - real(n, dp)) / real(n * (n - 1), dp)
      end if
      prior = 0.0_dp
      do j = 1, n
         do i = 1, n
            prior(i, j) = rbar * sqrt_variance(i) * sqrt_variance(j)
         end do
      end do
      do i = 1, n
         prior(i, i) = variance(i)
      end do
      phi_mat = matmul(transpose(y), y) / real(t, dp) - &
         2.0_dp * matmul(transpose(x), x) * sample / real(t, dp) + sample * sample
      phi = sum(phi_mat)
      term1 = matmul(transpose(x ** 3), x) / real(t, dp)
      rho_mat = term1
      do j = 1, n
         do i = 1, n
            rho_mat(i, j) = rho_mat(i, j) - variance(i) * sample(i, j)
         end do
      end do
      do i = 1, n
         rho_mat(i, i) = 0.0_dp
      end do
      rho_value = 0.0_dp
      do i = 1, n
         rho_value = rho_value + phi_mat(i, i)
      end do
      do j = 1, n
         do i = 1, n
            rho_value = rho_value + rbar * &
               sqrt_variance(j) / sqrt_variance(i) * rho_mat(i, j)
         end do
      end do
      gamma = frobenius_norm(sample - prior) ** 2
      if (gamma <= tiny(1.0_dp)) then
         shrinkage = 1.0_dp
      else
         kappa = (phi - rho_value) / gamma
         shrinkage = max(0.0_dp, min(1.0_dp, kappa / real(t, dp)))
      end if
      sigma = shrinkage * prior + (1.0_dp - shrinkage) * sample
      sigma = symmetrize(sigma)
   end function correlation_shrinkage_covariance

   function diagonal_shrinkage_covariance(rets) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp), allocatable :: x(:, :), y(:, :), sample(:, :), prior(:, :), phi_mat(:, :)
      real(dp) :: phi, rho_value, gamma, kappa, shrinkage
      integer :: i, t, n

      t = size(rets, 1)
      n = size(rets, 2)
      allocate(x(t, n), y(t, n), sample(n, n), prior(n, n), phi_mat(n, n))
      x = rets - spread(column_means(rets), 1, t)
      y = x * x
      sample = matmul(transpose(x), x) / real(max(t, 1), dp)
      prior = 0.0_dp
      do i = 1, n
         prior(i, i) = sample(i, i)
      end do
      phi_mat = matmul(transpose(y), y) / real(max(t, 1), dp) - &
         2.0_dp * matmul(transpose(x), x) * sample / real(max(t, 1), dp) + sample * sample
      phi = sum(phi_mat)
      rho_value = 0.0_dp
      do i = 1, n
         rho_value = rho_value + phi_mat(i, i)
      end do
      gamma = frobenius_norm(sample - prior) ** 2
      if (gamma <= tiny(1.0_dp)) then
         shrinkage = 1.0_dp
      else
         kappa = (phi - rho_value) / gamma
         shrinkage = max(0.0_dp, min(1.0_dp, kappa / real(max(t, 1), dp)))
      end if
      sigma = shrinkage * prior + (1.0_dp - shrinkage) * sample
      sigma = symmetrize(sigma)
   end function diagonal_shrinkage_covariance

   function one_parameter_covariance(rets) result(sigma)
      real(dp), intent(in) :: rets(:, :)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp), allocatable :: x(:, :), y(:, :), sample(:, :), prior(:, :), phi_mat(:, :)
      real(dp) :: mean_variance, phi, gamma, kappa, shrinkage
      integer :: i, t, n

      t = size(rets, 1)
      n = size(rets, 2)
      allocate(x(t, n), y(t, n), sample(n, n), prior(n, n), phi_mat(n, n))
      x = rets - spread(column_means(rets), 1, t)
      y = x * x
      sample = matmul(transpose(x), x) / real(max(t, 1), dp)
      mean_variance = 0.0_dp
      do i = 1, n
         mean_variance = mean_variance + sample(i, i)
      end do
      if (n > 0) mean_variance = mean_variance / real(n, dp)
      prior = 0.0_dp
      do i = 1, n
         prior(i, i) = mean_variance
      end do
      phi_mat = matmul(transpose(y), y) / real(max(t, 1), dp) - &
         2.0_dp * matmul(transpose(x), x) * sample / real(max(t, 1), dp) + sample * sample
      phi = sum(phi_mat)
      gamma = frobenius_norm(sample - prior) ** 2
      if (gamma <= tiny(1.0_dp)) then
         shrinkage = 1.0_dp
      else
         kappa = phi / gamma
         shrinkage = max(0.0_dp, min(1.0_dp, kappa / real(max(t, 1), dp)))
      end if
      sigma = shrinkage * prior + (1.0_dp - shrinkage) * sample
      sigma = symmetrize(sigma)
   end function one_parameter_covariance

   subroutine bayes_stein_covariance(rets, sigma_bs, info)
      real(dp), intent(in) :: rets(:, :)
      real(dp), intent(out) :: sigma_bs(:, :)
      integer, intent(out), optional :: info
      real(dp), allocatable :: inv_sigma_one(:), inv_sigma_dev(:)
      real(dp) :: sigma(size(rets, 2), size(rets, 2))
      real(dp) :: mu(size(rets, 2)), one(size(rets, 2))
      real(dp) :: w_min(size(rets, 2)), mu_min, phi, tau, denominator
      real(dp) :: a, b
      integer :: n, t, stat

      n = size(rets, 2)
      t = size(rets, 1)
      sigma = sample_covariance(rets)
      mu = column_means(rets)
      one = 1.0_dp
      call solve_linear(sigma, one, inv_sigma_one, stat)
      if (stat /= 0) then
         sigma_bs = sigma
         if (present(info)) info = stat
         return
      end if
      denominator = dot_product(one, inv_sigma_one)
      if (abs(denominator) <= tiny(1.0_dp)) then
         sigma_bs = sigma
         if (present(info)) info = -2
         return
      end if
      w_min = inv_sigma_one / denominator
      mu_min = dot_product(mu, w_min)
      call solve_linear(sigma, mu - mu_min, inv_sigma_dev, stat)
      if (stat /= 0) then
         sigma_bs = sigma
         if (present(info)) info = stat
         return
      end if
      phi = real(n + 2, dp) / (real(n + 2, dp) + real(t, dp) * &
         dot_product(mu - mu_min, inv_sigma_dev))
      phi = max(0.0_dp, min(1.0_dp - epsilon(1.0_dp), phi))
      tau = real(t, dp) * phi / (1.0_dp - phi)
      a = 1.0_dp + 1.0_dp / (real(t, dp) + tau)
      b = tau / (real(t, dp) * (real(t + 1, dp) + tau) * denominator)
      sigma_bs = a * sigma + b * outer_product(one, one)
      sigma_bs = symmetrize(sigma_bs)
      if (present(info)) info = 0
   end subroutine bayes_stein_covariance

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a), size(b))
      c = spread(a, 2, size(b)) * spread(b, 1, size(a))
   end function outer_product

end module riskportfolios_covariance
