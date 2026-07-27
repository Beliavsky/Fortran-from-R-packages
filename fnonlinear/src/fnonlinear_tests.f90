! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 1988-1990 Blake LeBaron
! Copyright (C) 2000 Adrian Trapletti
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
module fnonlinear_tests
  use chaos_kinds, only : dp
  use fnonlinear_rng, only : rng_state, fill_uniform
  use fnonlinear_linalg, only : ols_fit, pca_scores_standardized
  use fnonlinear_probability, only : normal_survival, chi_square_survival, f_survival
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private

  type, public :: bds_test_result
    real(dp), allocatable :: epsilon(:)
    real(dp), allocatable :: statistic(:, :)
    real(dp), allocatable :: p_value(:, :)
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: k_value(:)
    integer :: max_dimension = 0
    integer :: effective_n = 0
  end type bds_test_result

  type, public :: neural_test_result
    real(dp) :: chi_square = 0.0_dp
    real(dp) :: f_statistic = 0.0_dp
    real(dp) :: chi_square_p = 1.0_dp
    real(dp) :: f_p = 1.0_dp
    real(dp) :: ssr_null = 0.0_dp
    real(dp) :: ssr_alternative = 0.0_dp
    integer :: numerator_df = 0
    integer :: denominator_df = 0
    integer :: observations = 0
  end type neural_test_result

  type, public :: runs_test_result
    integer :: runs = 0
    integer :: negative_count = 0
    integer :: positive_count = 0
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
  end type runs_test_result

  type, public :: generic_test_result
    real(dp), allocatable :: statistic(:)
    real(dp), allocatable :: p_value(:)
    real(dp), allocatable :: parameter(:)
    character(len=32) :: method = ""
  end type generic_test_result

  public :: bds_test, white_neural_test, terasvirta_neural_test, runs_test
  public :: ts_test
contains
  subroutine bds_test(x, max_dimension, epsilon, result, status)
    real(dp), intent(in) :: x(:), epsilon(:)
    integer, intent(in) :: max_dimension
    type(bds_test_result), intent(out) :: result
    integer, intent(out) :: status
    real(dp), allocatable :: correlations(:)
    real(dp) :: k_value, sigma, std_error
    integer :: ieps, m, n_eff, local_status

    if (max_dimension < 2 .or. size(epsilon) < 1 .or. any(epsilon <= 0.0_dp) .or. &
        size(x) <= max_dimension + 1) then
      status = 1
      return
    end if
    n_eff = size(x) - max_dimension + 1
    allocate(result%epsilon(size(epsilon)))
    allocate(result%statistic(max_dimension - 1, size(epsilon)))
    allocate(result%p_value(max_dimension - 1, size(epsilon)))
    allocate(result%correlation(max_dimension, size(epsilon)))
    allocate(result%k_value(size(epsilon)))
    result%epsilon = epsilon
    result%max_dimension = max_dimension
    result%effective_n = n_eff
    do ieps = 1, size(epsilon)
      call bds_correlations(x, max_dimension, epsilon(ieps), correlations, k_value, local_status)
      if (local_status /= 0) then
        status = 2
        return
      end if
      result%correlation(:, ieps) = correlations
      result%k_value(ieps) = k_value
      do m = 2, max_dimension
        sigma = bds_sigma(correlations(1), k_value, m)
        if (sigma <= 0.0_dp) then
          result%statistic(m - 1, ieps) = quiet_nan()
          result%p_value(m - 1, ieps) = quiet_nan()
        else
          std_error = sqrt(sigma / real(n_eff, dp))
          result%statistic(m - 1, ieps) = &
            (correlations(m) - correlations(1)**m) / std_error
          result%p_value(m - 1, ieps) = &
            2.0_dp * normal_survival(abs(result%statistic(m - 1, ieps)))
        end if
      end do
    end do
    status = 0
  end subroutine bds_test

  subroutine bds_correlations(x, max_dimension, epsilon, correlations, k_value, status)
    real(dp), intent(in) :: x(:), epsilon
    integer, intent(in) :: max_dimension
    real(dp), allocatable, intent(out) :: correlations(:)
    real(dp), intent(out) :: k_value
    integer, intent(out) :: status
    integer :: n_eff, i, j, m, close_count
    integer, allocatable :: neighbor_counts(:)
    real(dp) :: ordered_count, phi
    logical :: close_pair

    n_eff = size(x) - max_dimension + 1
    if (n_eff < 3 .or. epsilon <= 0.0_dp) then
      allocate(correlations(0))
      k_value = quiet_nan()
      status = 1
      return
    end if
    allocate(correlations(max_dimension), neighbor_counts(n_eff))
    neighbor_counts = 1
    do i = 1, n_eff - 1
      do j = i + 1, n_eff
        if (abs(x(i) - x(j)) <= epsilon) then
          neighbor_counts(i) = neighbor_counts(i) + 1
          neighbor_counts(j) = neighbor_counts(j) + 1
        end if
      end do
    end do
    ordered_count = real(sum(neighbor_counts) - n_eff, dp)
    phi = real(sum(neighbor_counts * neighbor_counts), dp) - real(n_eff, dp) - &
      3.0_dp * ordered_count
    correlations(1) = ordered_count / real(n_eff * (n_eff - 1), dp)
    k_value = phi / real(n_eff * (n_eff - 1) * (n_eff - 2), dp)
    do m = 2, max_dimension
      close_count = 0
      do i = 1, n_eff - 1
        do j = i + 1, n_eff
          close_pair = all(abs(x(i:i + m - 1) - x(j:j + m - 1)) <= epsilon)
          if (close_pair) close_count = close_count + 1
        end do
      end do
      correlations(m) = 2.0_dp * real(close_count, dp) / real(n_eff * (n_eff - 1), dp)
    end do
    status = 0
  end subroutine bds_correlations

  pure real(dp) function bds_sigma(correlation_one, k_value, dimension) result(value)
    real(dp), intent(in) :: correlation_one, k_value
    integer, intent(in) :: dimension
    integer :: j
    value = 0.0_dp
    do j = 1, dimension - 1
      value = value + 2.0_dp * k_value**(dimension - j) * correlation_one**(2 * j)
    end do
    value = value + k_value**dimension + real((dimension - 1)**2, dp) * &
      correlation_one**(2 * dimension) - real(dimension**2, dp) * k_value * &
      correlation_one**(2 * dimension - 2)
    value = 4.0_dp * value
  end function bds_sigma

  subroutine white_neural_test(x, lag, qstar, q, weight_range, result, status, rng, gamma_out)
    real(dp), intent(in) :: x(:), weight_range
    integer, intent(in) :: lag, qstar, q
    type(neural_test_result), intent(out) :: result
    integer, intent(out) :: status
    type(rng_state), intent(inout), optional :: rng
    real(dp), allocatable, intent(out), optional :: gamma_out(:, :)
    type(rng_state) :: local_rng
    real(dp), allocatable :: response(:), predictors(:, :), base_design(:, :), base_coef(:)
    real(dp), allocatable :: residuals(:), gamma(:, :), phantom(:, :), scores(:, :), eigenvalues(:)
    real(dp), allocatable :: alternative_design(:, :), alternative_coef(:), alternative_residuals(:)
    real(dp) :: ssr0, ssr, half_range
    integer :: nobs, i, j, local_status, n

    n = size(x)
    nobs = n - lag
    if (lag < 1 .or. qstar < 1 .or. q < qstar + 1 .or. weight_range <= 0.0_dp .or. &
        nobs <= lag + qstar + 1) then
      status = 1
      return
    end if
    call build_autoregression_data(x, lag, response, predictors)
    allocate(base_design(nobs, lag + 1))
    base_design(:, 1) = 1.0_dp
    base_design(:, 2:) = predictors
    call ols_fit(base_design, response, base_coef, residuals, ssr0, local_status)
    if (local_status /= 0) then
      status = 2
      return
    end if
    allocate(gamma(lag + 1, q))
    half_range = 0.5_dp * weight_range
    if (present(rng)) then
      call fill_uniform(rng, gamma, -half_range, half_range)
    else
      call fill_uniform(local_rng, gamma, -half_range, half_range)
    end if
    allocate(phantom(nobs, q))
    phantom = matmul(base_design, gamma)
    do j = 1, q
      do i = 1, nobs
        phantom(i, j) = logistic_value(phantom(i, j))
      end do
    end do
    call pca_scores_standardized(phantom, scores, eigenvalues, local_status)
    if (local_status /= 0) then
      status = 3
      return
    end if
    allocate(alternative_design(nobs, 1 + lag + qstar))
    alternative_design(:, 1) = 1.0_dp
    alternative_design(:, 2:lag + 1) = predictors
    alternative_design(:, lag + 2:) = scores(:, 2:qstar + 1)
    call ols_fit(alternative_design, residuals, alternative_coef, alternative_residuals, &
      ssr, local_status)
    if (local_status /= 0 .or. ssr <= 0.0_dp .or. ssr >= ssr0) then
      status = 4
      return
    end if
    result%chi_square = real(n, dp) * log(ssr0 / ssr)
    result%denominator_df = n - lag - qstar
    result%f_statistic = ((ssr0 - ssr) / real(qstar, dp)) / &
      (ssr / real(result%denominator_df, dp))
    result%chi_square_p = chi_square_survival(result%chi_square, real(qstar, dp))
    result%f_p = f_survival(result%f_statistic, real(qstar, dp), &
      real(result%denominator_df, dp))
    result%ssr_null = ssr0
    result%ssr_alternative = ssr
    result%numerator_df = qstar
    result%observations = n
    if (present(gamma_out)) then
      allocate(gamma_out(size(gamma, 1), size(gamma, 2)))
      gamma_out = gamma
    end if
    status = 0
  end subroutine white_neural_test

  subroutine terasvirta_neural_test(x, lag, result, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    type(neural_test_result), intent(out) :: result
    integer, intent(out) :: status
    real(dp), allocatable :: response(:), predictors(:, :), base_design(:, :), base_coef(:)
    real(dp), allocatable :: residuals(:), alternative_design(:, :), alternative_coef(:)
    real(dp), allocatable :: alternative_residuals(:)
    real(dp) :: ssr0, ssr
    integer :: n, nobs, nquad, ncubic, nonlinear_terms, i, j, k, col, local_status

    n = size(x)
    nobs = n - lag
    nquad = lag * (lag + 1) / 2
    ncubic = lag * (lag + 1) * (lag + 2) / 6
    nonlinear_terms = nquad + ncubic
    if (lag < 1 .or. nobs <= lag + nonlinear_terms + 1) then
      status = 1
      return
    end if
    call build_autoregression_data(x, lag, response, predictors)
    allocate(base_design(nobs, lag + 1))
    base_design(:, 1) = 1.0_dp
    base_design(:, 2:) = predictors
    call ols_fit(base_design, response, base_coef, residuals, ssr0, local_status)
    if (local_status /= 0) then
      status = 2
      return
    end if
    allocate(alternative_design(nobs, 1 + lag + nonlinear_terms))
    alternative_design(:, 1) = 1.0_dp
    alternative_design(:, 2:lag + 1) = predictors
    col = lag + 1
    do i = 1, lag
      do j = i, lag
        col = col + 1
        alternative_design(:, col) = predictors(:, i) * predictors(:, j)
      end do
    end do
    do i = 1, lag
      do j = i, lag
        do k = j, lag
          col = col + 1
          alternative_design(:, col) = predictors(:, i) * predictors(:, j) * predictors(:, k)
        end do
      end do
    end do
    call ols_fit(alternative_design, residuals, alternative_coef, alternative_residuals, &
      ssr, local_status)
    if (local_status /= 0 .or. ssr <= 0.0_dp .or. ssr >= ssr0) then
      status = 3
      return
    end if
    result%chi_square = real(n, dp) * log(ssr0 / ssr)
    result%denominator_df = n - lag - nonlinear_terms
    result%f_statistic = ((ssr0 - ssr) / real(nonlinear_terms, dp)) / &
      (ssr / real(result%denominator_df, dp))
    result%chi_square_p = chi_square_survival(result%chi_square, real(nonlinear_terms, dp))
    result%f_p = f_survival(result%f_statistic, real(nonlinear_terms, dp), &
      real(result%denominator_df, dp))
    result%ssr_null = ssr0
    result%ssr_alternative = ssr
    result%numerator_df = nonlinear_terms
    result%observations = n
    status = 0
  end subroutine terasvirta_neural_test

  subroutine runs_test(x, result, status)
    real(dp), intent(in) :: x(:)
    type(runs_test_result), intent(out) :: result
    integer, intent(out) :: status
    integer, allocatable :: signs(:)
    integer :: i, n, n1, n2, runs
    real(dp) :: mean_runs, sd_runs

    n = count(abs(x) > 0.0_dp)
    if (n < 2) then
      status = 1
      return
    end if
    allocate(signs(n))
    n = 0
    do i = 1, size(x)
      if (abs(x(i)) > 0.0_dp) then
        n = n + 1
        if (x(i) < 0.0_dp) then
          signs(n) = -1
        else
          signs(n) = 1
        end if
      end if
    end do
    n1 = count(signs == -1)
    n2 = count(signs == 1)
    if (n1 == 0 .or. n2 == 0) then
      status = 2
      return
    end if
    runs = 1
    do i = 2, n
      if (signs(i) /= signs(i - 1)) runs = runs + 1
    end do
    mean_runs = 1.0_dp + 2.0_dp * real(n1 * n2, dp) / real(n1 + n2, dp)
    sd_runs = sqrt(2.0_dp * real(n1 * n2, dp) * &
      real(2 * n1 * n2 - n1 - n2, dp) / &
      (real((n1 + n2)**2, dp) * real(n1 + n2 - 1, dp)))
    result%runs = runs
    result%negative_count = n1
    result%positive_count = n2
    result%statistic = (real(runs, dp) - mean_runs) / sd_runs
    result%p_value = 2.0_dp * normal_survival(abs(result%statistic))
    status = 0
  end subroutine runs_test

  subroutine ts_test(x, method, result, status, max_dimension, epsilon, lag, qstar, q, &
      weight_range, rng)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: method
    type(generic_test_result), intent(out) :: result
    integer, intent(out) :: status
    integer, intent(in), optional :: max_dimension, lag, qstar, q
    real(dp), intent(in), optional :: epsilon(:), weight_range
    type(rng_state), intent(inout), optional :: rng
    type(bds_test_result) :: bds
    type(neural_test_result) :: neural
    integer :: md, lag_value, qstar_value, q_value, i, j, idx
    real(dp) :: range_value
    real(dp), allocatable :: eps_value(:)

    select case (lowercase(trim(adjustl(method))))
    case ("bds")
      md = 3
      if (present(max_dimension)) md = max_dimension
      if (present(epsilon)) then
        allocate(eps_value(size(epsilon)))
        eps_value = epsilon
      else
        allocate(eps_value(4))
        eps_value = standard_deviation(x) * [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
      end if
      call bds_test(x, md, eps_value, bds, status)
      if (status /= 0) return
      allocate(result%statistic(size(bds%statistic)), result%p_value(size(bds%p_value)))
      idx = 0
      do j = 1, size(bds%statistic, 2)
        do i = 1, size(bds%statistic, 1)
          idx = idx + 1
          result%statistic(idx) = bds%statistic(i, j)
          result%p_value(idx) = bds%p_value(i, j)
        end do
      end do
      allocate(result%parameter(1 + size(eps_value)))
      result%parameter(1) = real(md, dp)
      result%parameter(2:) = eps_value
      result%method = "BDS"
    case ("wnn")
      lag_value = 1
      qstar_value = 2
      q_value = 10
      range_value = 4.0_dp
      if (present(lag)) lag_value = lag
      if (present(qstar)) qstar_value = qstar
      if (present(q)) q_value = q
      if (present(weight_range)) range_value = weight_range
      if (present(rng)) then
        call white_neural_test(x, lag_value, qstar_value, q_value, range_value, &
          neural, status, rng=rng)
      else
        call white_neural_test(x, lag_value, qstar_value, q_value, range_value, &
          neural, status)
      end if
      if (status /= 0) return
      call pack_neural(neural, result, "White neural network")
    case ("tnn")
      lag_value = 1
      if (present(lag)) lag_value = lag
      call terasvirta_neural_test(x, lag_value, neural, status)
      if (status /= 0) return
      call pack_neural(neural, result, "Terasvirta neural network")
    case default
      status = 1
    end select
  contains
    subroutine pack_neural(source, target, name)
      type(neural_test_result), intent(in) :: source
      type(generic_test_result), intent(inout) :: target
      character(len=*), intent(in) :: name
      allocate(target%statistic(2), target%p_value(2), target%parameter(3))
      target%statistic = [source%chi_square, source%f_statistic]
      target%p_value = [source%chi_square_p, source%f_p]
      target%parameter = [real(source%numerator_df, dp), &
        real(source%denominator_df, dp), real(source%observations, dp)]
      target%method = name
    end subroutine pack_neural
  end subroutine ts_test

  subroutine build_autoregression_data(x, lag, response, predictors)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    real(dp), allocatable, intent(out) :: response(:), predictors(:, :)
    integer :: nobs, i, j
    nobs = size(x) - lag
    allocate(response(nobs), predictors(nobs, lag))
    do i = 1, nobs
      response(i) = x(i + lag)
      do j = 1, lag
        predictors(i, j) = x(i + lag - j)
      end do
    end do
  end subroutine build_autoregression_data

  pure real(dp) function logistic_value(x) result(value)
    real(dp), intent(in) :: x
    if (x >= 0.0_dp) then
      value = 1.0_dp / (1.0_dp + exp(-x))
    else
      value = exp(x) / (1.0_dp + exp(x))
    end if
  end function logistic_value

  pure real(dp) function standard_deviation(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: mean_value
    if (size(x) < 2) then
      value = 0.0_dp
      return
    end if
    mean_value = sum(x) / real(size(x), dp)
    value = sqrt(sum((x - mean_value)**2) / real(size(x) - 1, dp))
  end function standard_deviation

  pure function lowercase(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      else
        lower(i:i) = text(i:i)
      end if
    end do
  end function lowercase

  pure real(dp) function quiet_nan() result(value)
    value = ieee_value(0.0_dp, ieee_quiet_nan)
  end function quiet_nan
end module fnonlinear_tests
