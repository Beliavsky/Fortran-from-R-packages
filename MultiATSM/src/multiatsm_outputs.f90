! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_outputs
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : inverse_matrix, cholesky_lower, eye
  use multiatsm_types, only : response_result, variance_decomposition_result, forecast_result
  implicit none
  private

  public :: fitted_yields, forecast_yields, forecast_rmse
  public :: impulse_responses, generalized_impulse_responses
  public :: forecast_error_variance_decomposition, generalized_fevd
  public :: expected_short_rate_component, term_premium, forward_rates

contains

  subroutine fitted_yields(a, b, factors, yields, info)
    real(dp), intent(in) :: a(:), b(:, :), factors(:, :)
    real(dp), allocatable, intent(out) :: yields(:, :)
    integer, intent(out) :: info

    if (size(b, 1) /= size(a) .or. size(b, 2) /= size(factors, 1)) then
      info = -1
      allocate(yields(0, 0))
      return
    end if
    allocate(yields(size(a), size(factors, 2)))
    yields = spread(a, 2, size(factors, 2)) + matmul(b, factors)
    info = 0
  end subroutine fitted_yields

  subroutine forecast_yields(k0, k1, a, b, last_state, horizon, result, info)
    real(dp), intent(in) :: k0(:), k1(:, :), a(:), b(:, :), last_state(:)
    integer, intent(in) :: horizon
    type(forecast_result), intent(out) :: result
    integer, intent(out) :: info
    real(dp), allocatable :: state(:)
    integer :: h, k

    k = size(k0)
    if (horizon < 1 .or. size(k1, 1) /= k .or. size(k1, 2) /= k .or. size(last_state) /= k .or. &
        size(b, 1) /= size(a) .or. size(b, 2) /= k) then
      info = -1
      return
    end if
    allocate(result%factors(k, horizon), result%yields(size(a), horizon), state(k))
    state = last_state
    do h = 1, horizon
      state = k0 + matmul(k1, state)
      result%factors(:, h) = state
      result%yields(:, h) = a + matmul(b, state)
    end do
    info = 0
  end subroutine forecast_yields

  subroutine forecast_rmse(actual, predicted, rmse, info)
    real(dp), intent(in) :: actual(:, :), predicted(:, :)
    real(dp), allocatable, intent(out) :: rmse(:)
    integer, intent(out) :: info

    if (any(shape(actual) /= shape(predicted)) .or. size(actual, 2) < 1) then
      info = -1
      allocate(rmse(0))
      return
    end if
    allocate(rmse(size(actual, 1)))
    rmse = sqrt(sum((actual - predicted) ** 2, dim=2) / real(size(actual, 2), dp))
    info = 0
  end subroutine forecast_rmse

  subroutine reduced_form_covariance(sigma, g0, cov, impact, info)
    real(dp), intent(in) :: sigma(:, :), g0(:, :)
    real(dp), allocatable, intent(out) :: cov(:, :), impact(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: invg(:, :), chol(:, :)
    integer :: k

    k = size(sigma, 1)
    if (size(sigma, 2) /= k .or. size(g0, 1) /= k .or. size(g0, 2) /= k) then
      info = -1
      allocate(cov(0, 0), impact(0, 0))
      return
    end if
    call inverse_matrix(g0, invg, info, 1.0e-12_dp)
    if (info /= 0) return
    allocate(cov(k, k))
    cov = matmul(invg, matmul(sigma, transpose(invg)))
    call cholesky_lower(cov, chol, info, 1.0e-12_dp)
    if (info /= 0) return
    allocate(impact(k, k))
    impact = chol
  end subroutine reduced_form_covariance

  subroutine impulse_responses(phi, sigma, g0, b_loadings, horizon, result, info, impact_matrix)
    real(dp), intent(in) :: phi(:, :), sigma(:, :), g0(:, :), b_loadings(:, :)
    integer, intent(in) :: horizon
    type(response_result), intent(out) :: result
    integer, intent(out) :: info
    real(dp), intent(in), optional :: impact_matrix(:, :)
    real(dp), allocatable :: cov(:, :), impact(:, :), dyn(:, :)
    integer :: k, h

    k = size(phi, 1)
    if (horizon < 1 .or. size(phi, 2) /= k .or. size(b_loadings, 2) /= k) then
      info = -1
      return
    end if
    if (present(impact_matrix)) then
      if (size(impact_matrix, 1) /= k .or. size(impact_matrix, 2) /= k) then
        info = -2
        return
      end if
      allocate(impact(k, k), cov(k, k))
      impact = impact_matrix
      cov = matmul(impact, transpose(impact))
    else
      call reduced_form_covariance(sigma, g0, cov, impact, info)
      if (info /= 0) return
    end if
    allocate(result%factors(k, k, horizon), result%yields(size(b_loadings, 1), k, horizon), dyn(k, k))
    dyn = eye(k)
    do h = 1, horizon
      result%factors(:, :, h) = matmul(dyn, impact)
      result%yields(:, :, h) = matmul(b_loadings, result%factors(:, :, h))
      dyn = matmul(phi, dyn)
    end do
    info = 0
  end subroutine impulse_responses

  subroutine generalized_impulse_responses(phi, sigma, g0, b_loadings, horizon, result, info)
    real(dp), intent(in) :: phi(:, :), sigma(:, :), g0(:, :), b_loadings(:, :)
    integer, intent(in) :: horizon
    type(response_result), intent(out) :: result
    integer, intent(out) :: info
    real(dp), allocatable :: invg(:, :), cov(:, :), dyn(:, :), impact(:, :)
    integer :: k, h, j

    k = size(phi, 1)
    if (horizon < 1 .or. size(phi, 2) /= k .or. size(sigma, 1) /= k .or. size(sigma, 2) /= k .or. &
        size(g0, 1) /= k .or. size(g0, 2) /= k .or. size(b_loadings, 2) /= k) then
      info = -1
      return
    end if
    call inverse_matrix(g0, invg, info, 1.0e-12_dp)
    if (info /= 0) return
    allocate(cov(k, k), impact(k, k), dyn(k, k))
    cov = matmul(invg, matmul(sigma, transpose(invg)))
    do j = 1, k
      impact(:, j) = cov(:, j) / sqrt(max(cov(j, j), epsilon(1.0_dp)))
    end do
    allocate(result%factors(k, k, horizon), result%yields(size(b_loadings, 1), k, horizon))
    dyn = eye(k)
    do h = 1, horizon
      result%factors(:, :, h) = matmul(dyn, impact)
      result%yields(:, :, h) = matmul(b_loadings, result%factors(:, :, h))
      dyn = matmul(phi, dyn)
    end do
    info = 0
  end subroutine generalized_impulse_responses

  subroutine forecast_error_variance_decomposition(phi, sigma, g0, b_loadings, horizon, result, info, impact_matrix)
    real(dp), intent(in) :: phi(:, :), sigma(:, :), g0(:, :), b_loadings(:, :)
    integer, intent(in) :: horizon
    type(variance_decomposition_result), intent(out) :: result
    integer, intent(out) :: info
    real(dp), intent(in), optional :: impact_matrix(:, :)
    type(response_result) :: irf
    real(dp), allocatable :: num_f(:, :), num_y(:, :)
    integer :: k, ny, h, i
    real(dp) :: den

    if (present(impact_matrix)) then
      call impulse_responses(phi, sigma, g0, b_loadings, horizon, irf, info, impact_matrix)
    else
      call impulse_responses(phi, sigma, g0, b_loadings, horizon, irf, info)
    end if
    if (info /= 0) return
    k = size(phi, 1)
    ny = size(b_loadings, 1)
    allocate(result%factors(k, k, horizon), result%yields(ny, k, horizon))
    allocate(num_f(k, k), num_y(ny, k))
    num_f = 0.0_dp
    num_y = 0.0_dp
    do h = 1, horizon
      num_f = num_f + irf%factors(:, :, h) ** 2
      num_y = num_y + irf%yields(:, :, h) ** 2
      do i = 1, k
        den = sum(num_f(i, :))
        if (den > 0.0_dp) then
          result%factors(i, :, h) = num_f(i, :) / den
        else
          result%factors(i, :, h) = 0.0_dp
        end if
      end do
      do i = 1, ny
        den = sum(num_y(i, :))
        if (den > 0.0_dp) then
          result%yields(i, :, h) = num_y(i, :) / den
        else
          result%yields(i, :, h) = 0.0_dp
        end if
      end do
    end do
    info = 0
  end subroutine forecast_error_variance_decomposition

  subroutine generalized_fevd(phi, sigma, g0, b_loadings, horizon, result, info)
    real(dp), intent(in) :: phi(:, :), sigma(:, :), g0(:, :), b_loadings(:, :)
    integer, intent(in) :: horizon
    type(variance_decomposition_result), intent(out) :: result
    integer, intent(out) :: info
    real(dp), allocatable :: invg(:, :), cov(:, :), dyn(:, :), num_f(:, :), num_y(:, :)
    real(dp), allocatable :: den_f(:), den_y(:), response(:, :), yresponse(:, :), bdyn(:, :)
    integer :: k, ny, h, i, j
    real(dp) :: scale

    k = size(phi, 1)
    ny = size(b_loadings, 1)
    if (horizon < 1 .or. size(phi, 2) /= k) then
      info = -1
      return
    end if
    call inverse_matrix(g0, invg, info, 1.0e-12_dp)
    if (info /= 0) return
    allocate(cov(k, k), dyn(k, k), num_f(k, k), num_y(ny, k), den_f(k), den_y(ny))
    allocate(response(k, k), yresponse(ny, k), bdyn(ny, k))
    cov = matmul(invg, matmul(sigma, transpose(invg)))
    allocate(result%factors(k, k, horizon), result%yields(ny, k, horizon))
    num_f = 0.0_dp
    num_y = 0.0_dp
    den_f = 0.0_dp
    den_y = 0.0_dp
    dyn = eye(k)
    do h = 1, horizon
      response = matmul(dyn, cov)
      bdyn = matmul(b_loadings, dyn)
      yresponse = matmul(bdyn, cov)
      do j = 1, k
        scale = max(cov(j, j), epsilon(1.0_dp))
        num_f(:, j) = num_f(:, j) + response(:, j) ** 2 / scale
        num_y(:, j) = num_y(:, j) + yresponse(:, j) ** 2 / scale
      end do
      den_f = den_f + [(dot_product(dyn(i, :), matmul(cov, dyn(i, :))), i = 1, k)]
      do i = 1, ny
        den_y(i) = den_y(i) + dot_product(bdyn(i, :), matmul(cov, bdyn(i, :)))
      end do
      do i = 1, k
        if (den_f(i) > 0.0_dp) then
          result%factors(i, :, h) = num_f(i, :) / den_f(i)
          if (sum(result%factors(i, :, h)) > 0.0_dp) then
            result%factors(i, :, h) = result%factors(i, :, h) / sum(result%factors(i, :, h))
          end if
        else
          result%factors(i, :, h) = 0.0_dp
        end if
      end do
      do i = 1, ny
        if (den_y(i) > 0.0_dp) then
          result%yields(i, :, h) = num_y(i, :) / den_y(i)
          if (sum(result%yields(i, :, h)) > 0.0_dp) then
            result%yields(i, :, h) = result%yields(i, :, h) / sum(result%yields(i, :, h))
          end if
        else
          result%yields(i, :, h) = 0.0_dp
        end if
      end do
      dyn = matmul(phi, dyn)
    end do
    info = 0
  end subroutine generalized_fevd

  subroutine expected_short_rate_component(k0, k1, states, rho0, rho1, maturities, expected, info, floor_zero)
    real(dp), intent(in) :: k0(:), k1(:, :), states(:, :), rho0(:), rho1(:, :)
    integer, intent(in) :: maturities(:)
    real(dp), allocatable, intent(out) :: expected(:, :, :)
    integer, intent(out) :: info
    logical, intent(in), optional :: floor_zero
    real(dp), allocatable :: path(:), avg(:)
    integer :: k, nc, nt, nm, t, c, h, j
    logical :: use_floor
    real(dp) :: short_rate

    k = size(k0)
    nc = size(rho0)
    nt = size(states, 2)
    nm = size(maturities)
    if (size(k1, 1) /= k .or. size(k1, 2) /= k .or. size(states, 1) /= k .or. &
        size(rho1, 1) /= nc .or. size(rho1, 2) /= k .or. nm < 1 .or. any(maturities < 1)) then
      info = -1
      allocate(expected(0, 0, 0))
      return
    end if
    use_floor = .true.
    if (present(floor_zero)) use_floor = floor_zero
    allocate(expected(nc, nm, nt), path(k), avg(nc))
    do t = 1, nt
      do j = 1, nm
        path = states(:, t)
        avg = 0.0_dp
        do h = 1, maturities(j)
          do c = 1, nc
            short_rate = rho0(c) + dot_product(rho1(c, :), path)
            if (use_floor) short_rate = max(short_rate, 0.0_dp)
            avg(c) = avg(c) + short_rate
          end do
          path = k0 + matmul(k1, path)
        end do
        expected(:, j, t) = avg / real(maturities(j), dp)
      end do
    end do
    info = 0
  end subroutine expected_short_rate_component

  subroutine term_premium(observed_yields, expected_component, n_countries, premiums, info)
    real(dp), intent(in) :: observed_yields(:, :), expected_component(:, :, :)
    integer, intent(in) :: n_countries
    real(dp), allocatable, intent(out) :: premiums(:, :)
    integer, intent(out) :: info
    integer :: nc, nm, nt, c, j, row

    nc = n_countries
    nm = size(expected_component, 2)
    nt = size(expected_component, 3)
    if (size(expected_component, 1) /= nc .or. size(observed_yields, 1) /= nc * nm .or. &
        size(observed_yields, 2) /= nt) then
      info = -1
      allocate(premiums(0, 0))
      return
    end if
    allocate(premiums(nc * nm, nt))
    do c = 1, nc
      do j = 1, nm
        row = (c - 1) * nm + j
        premiums(row, :) = observed_yields(row, :) - expected_component(c, j, :)
      end do
    end do
    info = 0
  end subroutine term_premium

  subroutine forward_rates(yields, maturities, forwards, info)
    real(dp), intent(in) :: yields(:, :)
    integer, intent(in) :: maturities(:)
    real(dp), allocatable, intent(out) :: forwards(:, :)
    integer, intent(out) :: info
    integer :: nm, nt, j

    nm = size(maturities)
    nt = size(yields, 2)
    if (size(yields, 1) /= nm .or. nm < 2 .or. any(maturities < 1)) then
      info = -1
      allocate(forwards(0, 0))
      return
    end if
    allocate(forwards(nm - 1, nt))
    do j = 1, nm - 1
      if (maturities(j + 1) <= maturities(j)) then
        info = -2
        return
      end if
      forwards(j, :) = (real(maturities(j + 1), dp) * yields(j + 1, :) - &
        real(maturities(j), dp) * yields(j, :)) / real(maturities(j + 1) - maturities(j), dp)
    end do
    info = 0
  end subroutine forward_rates

end module multiatsm_outputs
