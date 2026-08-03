! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_var
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
  use multiatsm_kinds, only : dp
  use multiatsm_linalg, only : least_squares, eye
  use multiatsm_types, only : var_model, varx_country_model, gvar_model, &
    VARX_UNCONSTRAINED, VARX_SPANNED_RESTRICTED, VARX_FACTOR_RESTRICTED
  implicit none
  private

  public :: fit_var, fit_restricted_ols, build_star_factors
  public :: fit_varx_system, build_gvar, fit_gvar
  public :: transition_matrix_year, transition_matrix_mean

contains

  subroutine fit_var(x, model, info, restrictions)
    real(dp), intent(in) :: x(:, :)
    type(var_model), intent(out) :: model
    integer, intent(out) :: info
    real(dp), intent(in), optional :: restrictions(:, :)
    real(dp), allocatable :: lhs(:, :), rhs(:, :), coeff(:, :), design(:, :), solution(:, :)
    integer :: k, t

    k = size(x, 1)
    t = size(x, 2)
    if (k < 1 .or. t < 2) then
      info = -1
      return
    end if
    allocate(lhs(k, t - 1), rhs(k + 1, t - 1))
    lhs = x(:, 2:t)
    rhs(1, :) = 1.0_dp
    rhs(2:, :) = x(:, 1:t-1)

    if (present(restrictions)) then
      if (size(restrictions, 1) /= k .or. size(restrictions, 2) /= k + 1) then
        info = -2
        return
      end if
      call fit_restricted_ols(lhs, rhs, restrictions, coeff, info)
      if (info /= 0) return
    else
      allocate(design(t - 1, k + 1))
      design = transpose(rhs)
      call least_squares(design, transpose(lhs), solution, info)
      if (info /= 0) return
      allocate(coeff(k, k + 1))
      coeff = transpose(solution)
    end if

    allocate(model%intercept(k), model%phi(k, k), model%residuals(k, t - 1), model%sigma(k, k))
    model%intercept = coeff(:, 1)
    model%phi = coeff(:, 2:)
    model%residuals = lhs - matmul(coeff, rhs)
    model%sigma = matmul(model%residuals, transpose(model%residuals)) / real(t - 1, dp)
    info = 0
  end subroutine fit_var

  subroutine fit_restricted_ols(lhs, rhs, restrictions, coefficients, info)
    real(dp), intent(in) :: lhs(:, :), rhs(:, :), restrictions(:, :)
    real(dp), allocatable, intent(out) :: coefficients(:, :)
    integer, intent(out) :: info
    integer, allocatable :: free_idx(:), fixed_idx(:)
    real(dp), allocatable :: y(:, :), design(:, :), sol(:, :)
    logical, allocatable :: free_mask(:)
    integer :: m, n, t, i, nf, nx

    m = size(lhs, 1)
    t = size(lhs, 2)
    n = size(rhs, 1)
    if (size(rhs, 2) /= t .or. size(restrictions, 1) /= m .or. size(restrictions, 2) /= n) then
      info = -1
      allocate(coefficients(0, 0))
      return
    end if

    allocate(coefficients(m, n))
    coefficients = 0.0_dp
    do i = 1, m
      allocate(free_mask(n))
      free_mask = ieee_is_nan(restrictions(i, :))
      free_idx = pack([(nx, nx = 1, n)], free_mask)
      fixed_idx = pack([(nx, nx = 1, n)], .not. free_mask)
      coefficients(i, fixed_idx) = restrictions(i, fixed_idx)
      nf = size(free_idx)
      allocate(y(t, 1))
      y(:, 1) = lhs(i, :)
      if (size(fixed_idx) > 0) then
        y(:, 1) = y(:, 1) - matmul(transpose(rhs(fixed_idx, :)), restrictions(i, fixed_idx))
      end if
      if (nf > 0) then
        allocate(design(t, nf))
        design = transpose(rhs(free_idx, :))
        call least_squares(design, y, sol, info)
        if (info /= 0) return
        coefficients(i, free_idx) = sol(:, 1)
        deallocate(design, sol)
      end if
      deallocate(y, free_mask, free_idx, fixed_idx)
    end do
    info = 0
  end subroutine fit_restricted_ols

  subroutine build_star_factors(z, weights, zstar, info)
    real(dp), intent(in) :: z(:, :, :)
    real(dp), intent(in) :: weights(:, :)
    real(dp), allocatable, intent(out) :: zstar(:, :, :)
    integer, intent(out) :: info
    integer :: c, j, k, t, nc, nv, nt

    nc = size(z, 1)
    nv = size(z, 2)
    nt = size(z, 3)
    if (size(weights, 1) /= nc .or. size(weights, 2) /= nc) then
      info = -1
      allocate(zstar(0, 0, 0))
      return
    end if
    allocate(zstar(nc, nv, nt))
    zstar = 0.0_dp
    do c = 1, nc
      do j = 1, nc
        do k = 1, nv
          do t = 1, nt
            zstar(c, k, t) = zstar(c, k, t) + weights(c, j) * z(j, k, t)
          end do
        end do
      end do
    end do
    info = 0
  end subroutine build_star_factors

  subroutine fit_varx_system(z, zstar, global_factors, n_spanned, constraint_kind, countries, info, constrained_index)
    real(dp), intent(in) :: z(:, :, :), zstar(:, :, :), global_factors(:, :)
    integer, intent(in) :: n_spanned, constraint_kind
    type(varx_country_model), allocatable, intent(out) :: countries(:)
    integer, intent(out) :: info
    integer, intent(in), optional :: constrained_index
    real(dp), allocatable :: lhs(:, :), rhs(:, :), restrictions(:, :), coeff(:, :)
    real(dp) :: nan_value
    integer :: nc, k, nt, g, c, m, idx, ci

    nc = size(z, 1)
    k = size(z, 2)
    nt = size(z, 3)
    g = size(global_factors, 1)
    if (size(zstar, 1) /= nc .or. size(zstar, 2) /= k .or. size(zstar, 3) /= nt) then
      info = -1
      return
    end if
    if (size(global_factors, 2) /= nt .or. n_spanned < 0 .or. n_spanned > k) then
      info = -2
      return
    end if
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    allocate(countries(nc))
    m = k - n_spanned
    ci = 0
    if (present(constrained_index)) ci = constrained_index

    do c = 1, nc
      allocate(lhs(k, nt - 1), rhs(1 + 2 * k + g, nt - 1))
      lhs = z(c, :, 2:nt)
      rhs(1, :) = 1.0_dp
      rhs(2:k+1, :) = z(c, :, 1:nt-1)
      rhs(k+2:2*k+1, :) = zstar(c, :, 1:nt-1)
      if (g > 0) rhs(2*k+2:, :) = global_factors(:, 1:nt-1)

      if (constraint_kind == VARX_UNCONSTRAINED) then
        allocate(restrictions(k, 1 + 2 * k + g))
        restrictions = nan_value
      else
        allocate(restrictions(k, 1 + 2 * k + g))
        restrictions = nan_value
        if (constraint_kind == VARX_SPANNED_RESTRICTED .and. n_spanned > 0) then
          restrictions(:, k + 2 + m:2 * k + 1) = 0.0_dp
        else if (constraint_kind == VARX_FACTOR_RESTRICTED) then
          if (ci < 1 .or. ci > k) then
            info = -3
            return
          end if
          restrictions(ci, :) = 0.0_dp
          restrictions(ci, 1) = nan_value
          restrictions(ci, 1 + ci) = nan_value
          restrictions(ci, 1 + k + ci) = nan_value
        end if
      end if

      call fit_restricted_ols(lhs, rhs, restrictions, coeff, info)
      if (info /= 0) return
      allocate(countries(c)%phi0(k), countries(c)%phi1(k, k), countries(c)%phi1_star(k, k))
      allocate(countries(c)%phi_global(k, g), countries(c)%phi0_star(k, k))
      allocate(countries(c)%residuals(k, nt - 1), countries(c)%sigma(k, k))
      countries(c)%phi0 = coeff(:, 1)
      countries(c)%phi1 = coeff(:, 2:k+1)
      countries(c)%phi1_star = coeff(:, k+2:2*k+1)
      if (g > 0) countries(c)%phi_global = coeff(:, 2*k+2:)
      countries(c)%phi0_star = 0.0_dp
      countries(c)%residuals = lhs - matmul(coeff, rhs)
      countries(c)%sigma = matmul(countries(c)%residuals, transpose(countries(c)%residuals)) / real(nt, dp)
      if (constraint_kind == VARX_FACTOR_RESTRICTED .and. ci >= 1 .and. ci <= k) then
        do idx = 1, k
          if (idx /= ci) then
            countries(c)%sigma(ci, idx) = 0.0_dp
            countries(c)%sigma(idx, ci) = 0.0_dp
          end if
        end do
      end if
      deallocate(lhs, rhs, restrictions, coeff)
    end do
    info = 0
  end subroutine fit_varx_system

  subroutine build_gvar(country_models, global_model, weights, model, info)
    type(varx_country_model), intent(in) :: country_models(:)
    type(var_model), intent(in) :: global_model
    real(dp), intent(in) :: weights(:, :)
    type(gvar_model), intent(out) :: model
    integer, intent(out) :: info
    integer :: nc, k, g, dim, c, j, r0, r1, q0, q1

    nc = size(country_models)
    if (nc < 1) then
      info = -1
      return
    end if
    k = size(country_models(1)%phi0)
    g = size(global_model%intercept)
    if (size(weights, 1) /= nc .or. size(weights, 2) /= nc) then
      info = -2
      return
    end if
    dim = g + nc * k
    allocate(model%country(nc))
    model%country = country_models
    model%global_model = global_model
    allocate(model%gy0(dim, dim), model%f0(dim), model%f1(dim, dim), model%sigma_y(dim, dim))
    model%gy0 = eye(dim)
    model%f0 = 0.0_dp
    model%f1 = 0.0_dp
    model%sigma_y = 0.0_dp

    if (g > 0) then
      model%f0(1:g) = global_model%intercept
      model%f1(1:g, 1:g) = global_model%phi
      model%sigma_y(1:g, 1:g) = global_model%sigma
    end if

    do c = 1, nc
      r0 = g + (c - 1) * k + 1
      r1 = g + c * k
      model%f0(r0:r1) = country_models(c)%phi0
      if (g > 0) model%f1(r0:r1, 1:g) = country_models(c)%phi_global
      do j = 1, nc
        q0 = g + (j - 1) * k + 1
        q1 = g + j * k
        model%f1(r0:r1, q0:q1) = weights(c, j) * country_models(c)%phi1_star
        if (j == c) model%f1(r0:r1, q0:q1) = model%f1(r0:r1, q0:q1) + country_models(c)%phi1
      end do
      model%sigma_y(r0:r1, r0:r1) = country_models(c)%sigma
    end do
    info = 0
  end subroutine build_gvar

  subroutine fit_gvar(z, global_factors, weights, n_spanned, constraint_kind, model, info, constrained_index)
    real(dp), intent(in) :: z(:, :, :), global_factors(:, :), weights(:, :)
    integer, intent(in) :: n_spanned, constraint_kind
    type(gvar_model), intent(out) :: model
    integer, intent(out) :: info
    integer, intent(in), optional :: constrained_index
    real(dp), allocatable :: zstar(:, :, :)
    type(varx_country_model), allocatable :: countries(:)
    type(var_model) :: global_model
    integer :: g

    call build_star_factors(z, weights, zstar, info)
    if (info /= 0) return
    if (present(constrained_index)) then
      call fit_varx_system(z, zstar, global_factors, n_spanned, constraint_kind, countries, info, constrained_index)
    else
      call fit_varx_system(z, zstar, global_factors, n_spanned, constraint_kind, countries, info)
    end if
    if (info /= 0) return

    g = size(global_factors, 1)
    if (g > 0) then
      call fit_var(global_factors, global_model, info)
      if (info /= 0) return
    else
      allocate(global_model%intercept(0), global_model%phi(0, 0), global_model%sigma(0, 0))
      allocate(global_model%residuals(0, max(0, size(global_factors, 2)-1)))
    end if
    call build_gvar(countries, global_model, weights, model, info)
  end subroutine fit_gvar

  subroutine transition_matrix_year(flows, weights, info)
    real(dp), intent(in) :: flows(:, :)
    real(dp), allocatable, intent(out) :: weights(:, :)
    integer, intent(out) :: info
    integer :: n, i
    real(dp) :: total

    n = size(flows, 1)
    if (size(flows, 2) /= n) then
      info = -1
      allocate(weights(0, 0))
      return
    end if
    allocate(weights(n, n))
    do i = 1, n
      total = sum(flows(i, :))
      if (.not. ieee_is_finite(total) .or. total <= 0.0_dp .or. any(.not. ieee_is_finite(flows(i, :)))) then
        weights = ieee_value(0.0_dp, ieee_quiet_nan)
        info = 1
        return
      end if
      weights(i, :) = flows(i, :) / total
    end do
    info = 0
  end subroutine transition_matrix_year

  subroutine transition_matrix_mean(flows_by_year, first_year, last_year, weights, info)
    real(dp), intent(in) :: flows_by_year(:, :, :)
    integer, intent(in) :: first_year, last_year
    real(dp), allocatable, intent(out) :: weights(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: wy(:, :)
    integer :: y, n, count

    n = size(flows_by_year, 1)
    if (size(flows_by_year, 2) /= n .or. first_year < 1 .or. last_year > size(flows_by_year, 3) .or. &
        first_year > last_year) then
      info = -1
      allocate(weights(0, 0))
      return
    end if
    allocate(weights(n, n))
    weights = 0.0_dp
    count = 0
    do y = first_year, last_year
      call transition_matrix_year(flows_by_year(:, :, y), wy, info)
      if (info /= 0) return
      weights = weights + wy
      count = count + 1
    end do
    weights = weights / real(count, dp)
    info = 0
  end subroutine transition_matrix_mean

end module multiatsm_var
