! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_likelihood
  use multiatsm_kinds, only : dp, pi
  use multiatsm_linalg, only : inverse_matrix, logdet_spd
  use multiatsm_types, only : atsm_likelihood_result
  implicit none
  private

  public :: gaussian_log_density, atsm_log_likelihood, yield_error_variance

contains

  subroutine gaussian_log_density(residuals, sigma, log_density, info)
    real(dp), intent(in) :: residuals(:, :), sigma(:, :)
    real(dp), allocatable, intent(out) :: log_density(:)
    integer, intent(out) :: info
    real(dp), allocatable :: invsigma(:, :), solved(:, :)
    real(dp) :: logdet
    integer :: n, t, j

    n = size(residuals, 1)
    t = size(residuals, 2)
    if (size(sigma, 1) /= n .or. size(sigma, 2) /= n) then
      info = -1
      allocate(log_density(0))
      return
    end if
    call inverse_matrix(sigma, invsigma, info, 1.0e-12_dp)
    if (info /= 0) then
      allocate(log_density(0))
      return
    end if
    logdet = logdet_spd(sigma, info)
    if (info /= 0) then
      allocate(log_density(0))
      return
    end if
    allocate(solved(n, t), log_density(t))
    solved = matmul(invsigma, residuals)
    do j = 1, t
      log_density(j) = -0.5_dp * real(n, dp) * log(2.0_dp * pi) - 0.5_dp * logdet - &
        0.5_dp * dot_product(residuals(:, j), solved(:, j))
    end do
    info = 0
  end subroutine gaussian_log_density

  subroutine atsm_log_likelihood(yields, states, spanned, wpca, we, a_loadings, b_loadings, &
      k0, k1, sigma_states, n_countries, result, info)
    real(dp), intent(in) :: yields(:, :), states(:, :), spanned(:, :)
    real(dp), intent(in) :: wpca(:, :), we(:, :), a_loadings(:), b_loadings(:, :)
    real(dp), intent(in) :: k0(:), k1(:, :), sigma_states(:, :)
    integer, intent(in) :: n_countries
    type(atsm_likelihood_result), intent(out) :: result
    integer, intent(out) :: info
    real(dp), allocatable :: predicted_yields(:, :), observed_error_portfolios(:, :)
    real(dp), allocatable :: predicted_error_portfolios(:, :), e_q(:, :), e_p(:, :)
    real(dp), allocatable :: ld_q(:), ld_p(:), sigma_q(:, :), gram(:, :)
    integer :: jy, ts, ks, np, ne, nc, c, e0, e1, ne_country, t_use, i
    real(dp) :: jacobian, se_c

    jy = size(yields, 1)
    ts = size(yields, 2)
    ks = size(states, 1)
    np = size(spanned, 1)
    ne = size(we, 1)
    nc = n_countries
    if (ts < 2 .or. nc < 1 .or. size(states, 2) /= ts .or. size(spanned, 2) /= ts .or. &
        size(wpca, 1) /= np .or. size(wpca, 2) /= jy .or. size(we, 2) /= jy .or. &
        size(a_loadings) /= jy .or. size(b_loadings, 1) /= jy .or. size(b_loadings, 2) /= np .or. &
        size(k0) /= ks .or. size(k1, 1) /= ks .or. size(k1, 2) /= ks .or. &
        size(sigma_states, 1) /= ks .or. size(sigma_states, 2) /= ks .or. mod(ne, nc) /= 0) then
      info = -1
      return
    end if

    t_use = ts - 1
    ne_country = ne / nc
    allocate(predicted_yields(jy, ts), observed_error_portfolios(ne, ts))
    allocate(predicted_error_portfolios(ne, ts), e_q(ne, t_use), e_p(ks, t_use))
    predicted_yields = spread(a_loadings, 2, ts) + matmul(b_loadings, spanned)
    observed_error_portfolios = matmul(we, yields)
    predicted_error_portfolios = matmul(we, predicted_yields)
    e_q = observed_error_portfolios(:, 2:ts) - predicted_error_portfolios(:, 2:ts)
    e_p = states(:, 2:ts) - spread(k0, 2, t_use) - matmul(k1, states(:, 1:ts-1))

    allocate(result%yield_error_variance(nc), result%yield_residuals(ne, t_use), &
      result%factor_residuals(ks, t_use), sigma_q(ne, ne))
    result%yield_residuals = e_q
    result%factor_residuals = e_p
    sigma_q = 0.0_dp
    do c = 1, nc
      e0 = (c - 1) * ne_country + 1
      e1 = c * ne_country
      se_c = sqrt(sum(e_q(e0:e1, :) ** 2) / real(ne_country * t_use, dp))
      se_c = max(se_c, sqrt(epsilon(1.0_dp)))
      result%yield_error_variance(c) = se_c * se_c
    end do
    ! Replace the compact construction above with explicit diagonal assignment.
    sigma_q = 0.0_dp
    do c = 1, nc
      e0 = (c - 1) * ne_country + 1
      e1 = c * ne_country
      do i = e0, e1
        sigma_q(i, i) = result%yield_error_variance(c)
      end do
    end do

    call gaussian_log_density(e_q, sigma_q, ld_q, info)
    if (info /= 0) return
    call gaussian_log_density(e_p, sigma_states, ld_p, info)
    if (info /= 0) return

    allocate(gram(np, np))
    gram = matmul(wpca, transpose(wpca))
    jacobian = 0.5_dp * logdet_spd(gram, info)
    if (info /= 0) return

    result%negative_log_likelihood = -sum(ld_q + ld_p + jacobian)
    result%mean_negative_log_likelihood = result%negative_log_likelihood / real(t_use, dp)
    allocate(result%a_loadings(jy), result%b_loadings(jy, np))
    result%a_loadings = a_loadings
    result%b_loadings = b_loadings
    info = 0
  end subroutine atsm_log_likelihood

  subroutine yield_error_variance(wpca_full, n_spanned_per_country, measurement_variance, &
      n_countries, variances, covariance, info)
    real(dp), intent(in) :: wpca_full(:, :), measurement_variance(:)
    integer, intent(in) :: n_spanned_per_country, n_countries
    real(dp), allocatable, intent(out) :: variances(:), covariance(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: invw(:, :), sigma_portfolio(:, :)
    integer :: n, j, nc, c, r0, r1, i

    n = size(wpca_full, 1)
    nc = n_countries
    if (size(wpca_full, 2) /= n .or. nc < 1 .or. mod(n, nc) /= 0 .or. &
        size(measurement_variance) /= nc .or. n_spanned_per_country < 0) then
      info = -1
      allocate(variances(0), covariance(0, 0))
      return
    end if
    j = n / nc
    if (n_spanned_per_country > j) then
      info = -2
      allocate(variances(0), covariance(0, 0))
      return
    end if
    call inverse_matrix(wpca_full, invw, info, 1.0e-12_dp)
    if (info /= 0) then
      allocate(variances(0), covariance(0, 0))
      return
    end if
    allocate(sigma_portfolio(n, n), covariance(n, n), variances(n))
    sigma_portfolio = 0.0_dp
    do c = 1, nc
      r0 = (c - 1) * j + n_spanned_per_country + 1
      r1 = c * j
      do i = r0, r1
        sigma_portfolio(i, i) = measurement_variance(c)
      end do
    end do
    covariance = matmul(invw, matmul(sigma_portfolio, transpose(invw)))
    variances = [(covariance(i, i), i = 1, n)]
    info = 0
  end subroutine yield_error_variance

end module multiatsm_likelihood
