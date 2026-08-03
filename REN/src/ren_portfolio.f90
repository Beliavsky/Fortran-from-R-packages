! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren_portfolio
  use ren_kinds, only : dp
  use ren_types, only : asset_group, cluster_result, ren_success, ren_invalid_argument, &
    ren_dimension_error, ren_numerical_error, ren_dependency_error
  use ren_linalg, only : column_variances, covariance_matrix, correlation_matrix, &
    equality_minimum_variance, long_only_minimum_variance, least_squares
  use ren_regularization, only : fit_regularized_cv, fit_elastic_net_cv
  use ren_random, only : initialize_random_seed, sample_without_replacement, random_choice
  use corpcor, only : covariance_shrinkage, matrix_shrinkage_result, fast_svd, svd_result, &
    make_positive_definite, pseudoinverse, corpcor_success
  implicit none
  private
  public :: insert_at, po_avg, po_gross_exp, po_cov_shrink, po_cols, po_jm
  public :: buh_clust, po_bhu, po_tzt, po_sw, po_sw_lasso
contains
  function insert_at(a, position, value) result(out)
    real(dp), intent(in) :: a(:), value
    integer, intent(in) :: position
    real(dp), allocatable :: out(:)
    integer :: pos
    pos = max(0, min(size(a), position))
    allocate(out(size(a) + 1))
    if (pos > 0) out(1:pos) = a(1:pos)
    out(pos + 1) = value
    if (pos < size(a)) out(pos + 2:) = a(pos + 1:)
  end function insert_at

  subroutine po_cols(y0, x0, weights, status, variance_tolerance)
    real(dp), intent(in) :: y0(:), x0(:, :)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    real(dp), allocatable :: x(:, :), y(:), design(:, :), beta(:), active_weights(:)
    integer, allocatable :: active(:)
    integer :: p, j, istat
    call filter_columns(x0, x, active, variance_tolerance)
    p = size(x, 2)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    if (size(y0) /= size(x0, 1)) then
      status = ren_dimension_error
      return
    end if
    if (p == 0) then
      weights = 1.0_dp / real(size(weights), dp)
      status = ren_success
      return
    else if (p == 1) then
      weights(active(1)) = 1.0_dp
      status = ren_success
      return
    end if
    allocate(y(size(y0)), design(size(x, 1), p - 1))
    y = x(:, 1) - y0
    do j = 2, p
      design(:, j - 1) = x(:, 1) - x(:, j)
    end do
    call least_squares(design, y, beta, istat)
    active_weights = insert_at(beta, 0, 1.0_dp - sum(beta))
    weights(active) = active_weights
    status = istat
  end subroutine po_cols

  subroutine po_jm(x0, weights, status, variance_tolerance)
    real(dp), intent(in) :: x0(:, :)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    real(dp), allocatable :: x(:, :), covariance(:, :), active_weights(:)
    integer, allocatable :: active(:)
    integer :: istat
    call filter_columns(x0, x, active, variance_tolerance)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    if (size(active) == 0) then
      weights = 1.0_dp / real(size(weights), dp)
      status = ren_success
      return
    end if
    covariance = covariance_matrix(x, istat)
    call long_only_minimum_variance(covariance, active_weights, istat)
    weights(active) = active_weights
    status = istat
  end subroutine po_jm

  subroutine po_avg(y0, x0, method, weights, status, variance_tolerance, seed)
    real(dp), intent(in) :: y0(:), x0(:, :)
    character(len=*), intent(in) :: method
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    integer, intent(in), optional :: seed
    real(dp), allocatable :: x(:, :), design(:, :), response(:), beta(:), candidates(:, :), row_weights(:)
    real(dp) :: alpha_best
    integer, allocatable :: active(:)
    integer :: p, i, j, k, istat, aggregate_status, use_seed
    call filter_columns(x0, x, active, variance_tolerance)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    if (size(y0) /= size(x0, 1)) then
      status = ren_dimension_error
      return
    end if
    p = size(x, 2)
    if (p == 0) then
      weights = 1.0_dp / real(size(weights), dp)
      status = ren_success
      return
    else if (p == 1) then
      weights(active(1)) = 1.0_dp
      status = ren_success
      return
    end if
    allocate(candidates(p, p), response(size(y0)), design(size(y0), p - 1))
    candidates = 0.0_dp
    aggregate_status = ren_success
    use_seed = 1009
    if (present(seed)) use_seed = seed
    do i = 1, p
      response = x(:, i) - y0
      k = 0
      do j = 1, p
        if (j == i) cycle
        k = k + 1
        design(:, k) = x(:, i) - x(:, j)
      end do
      select case (trim(adjustl(method)))
      case ('LASSO')
        call fit_regularized_cv(design, response, 1.0_dp, beta, istat, seed=use_seed + i)
      case ('RIDGE')
        call fit_regularized_cv(design, response, 0.0_dp, beta, istat, seed=use_seed + i)
      case ('EN')
        call fit_elastic_net_cv(design, response, beta, alpha_best, istat, seed=use_seed + i)
      case default
        allocate(beta(p - 1))
        beta = 0.0_dp
        istat = ren_invalid_argument
      end select
      if (istat /= ren_success) aggregate_status = istat
      row_weights = insert_at(beta, i - 1, 1.0_dp - sum(beta))
      candidates(i, :) = row_weights
    end do
    weights(active) = sum(candidates, dim=1) / real(p, dp)
    status = aggregate_status
  end subroutine po_avg

  subroutine po_gross_exp(y0, x0, method, weights, status, variance_tolerance, seed)
    real(dp), intent(in) :: y0(:), x0(:, :)
    character(len=*), intent(in) :: method
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    integer, intent(in), optional :: seed
    real(dp), allocatable :: x(:, :), covariance(:, :), solution(:), portfolio(:)
    real(dp), allocatable :: design(:, :), response(:), beta(:), active_weights(:)
    integer, allocatable :: active(:)
    integer :: p, j, k, istat, use_seed
    call filter_columns(x0, x, active, variance_tolerance)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    if (size(y0) /= size(x0, 1)) then
      status = ren_dimension_error
      return
    end if
    p = size(x, 2)
    if (p == 0) then
      weights = 1.0_dp / real(size(weights), dp)
      status = ren_success
      return
    else if (p == 1) then
      weights(active(1)) = 1.0_dp
      status = ren_success
      return
    end if
    use_seed = 2003
    if (present(seed)) use_seed = seed
    select case (trim(adjustl(method)))
    case ('NOSHORT')
      covariance = covariance_matrix(x, istat)
      call long_only_minimum_variance(covariance, solution, istat)
      k = maxloc(solution, dim=1)
      if (abs(solution(k) - 1.0_dp) <= 5.0e-6_dp) then
        allocate(design(size(x, 1), p - 1), response(size(x, 1)))
        response = x(:, k) - y0
        j = 0
        do istat = 1, p
          if (istat == k) cycle
          j = j + 1
          design(:, j) = x(:, k) - x(:, istat)
        end do
        call fit_regularized_cv(design, response, 1.0_dp, beta, status, seed=use_seed)
        active_weights = insert_at(beta, k - 1, 1.0_dp - sum(beta))
      else
        portfolio = matmul(x, solution)
        allocate(design(size(x, 1), p), response(size(x, 1)))
        response = portfolio - y0
        do j = 1, p
          design(:, j) = portfolio - x(:, j)
        end do
        call fit_regularized_cv(design, response, 1.0_dp, beta, status, seed=use_seed)
        active_weights = (1.0_dp - sum(beta)) * solution + beta
      end if
    case ('EQUAL')
      allocate(solution(p))
      solution = 1.0_dp / real(p, dp)
      portfolio = matmul(x, solution)
      allocate(design(size(x, 1), p), response(size(x, 1)))
      response = portfolio - y0
      do j = 1, p
        design(:, j) = portfolio - x(:, j)
      end do
      call fit_regularized_cv(design, response, 1.0_dp, beta, status, seed=use_seed)
      active_weights = (1.0_dp - sum(beta)) * solution + beta
    case default
      allocate(active_weights(p))
      active_weights = 1.0_dp / real(p, dp)
      status = ren_invalid_argument
    end select
    weights(active) = active_weights
  end subroutine po_gross_exp

  subroutine po_cov_shrink(y0, x0, weights, status, variance_tolerance)
    real(dp), intent(in) :: y0(:), x0(:, :)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    real(dp), allocatable :: x(:, :), covariance(:, :), active_weights(:)
    integer, allocatable :: active(:)
    type(matrix_shrinkage_result) :: shrunk
    integer :: istat
    call filter_columns(x0, x, active, variance_tolerance)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    if (size(y0) /= size(x0, 1)) then
      status = ren_dimension_error
      return
    end if
    if (size(active) == 0) then
      weights = 1.0_dp / real(size(weights), dp)
      status = ren_success
      return
    end if
    covariance = covariance_matrix(x, istat)
    ! Source parity: REN calls corpcor::cov.shrink(cov(x0)), so the empirical
    ! covariance matrix is treated as the observations passed to corpcor.
    shrunk = covariance_shrinkage(covariance)
    if (shrunk%status /= corpcor_success) then
      status = ren_dependency_error
      return
    end if
    call equality_minimum_variance(shrunk%value, active_weights, istat)
    weights(active) = active_weights
    status = istat
  end subroutine po_cov_shrink

  function buh_clust(x0, variance_tolerance) result(result)
    real(dp), intent(in) :: x0(:, :)
    real(dp), intent(in), optional :: variance_tolerance
    type(cluster_result) :: result
    real(dp), allocatable :: x(:, :), correlation(:, :), next_correlation(:, :)
    real(dp) :: rho, best_rho
    integer, allocatable :: active(:)
    type(asset_group), allocatable :: groups(:), next_groups(:), best_groups(:)
    integer :: p, i, j, first, second, step, istat, out_index
    call filter_columns(x0, x, active, variance_tolerance)
    p = size(x, 2)
    if (p == 0) then
      allocate(result%groups(0), result%merge_correlation(0))
      result%status = ren_invalid_argument
      return
    end if
    allocate(groups(p))
    do i = 1, p
      allocate(groups(i)%index(1))
      groups(i)%index = i
    end do
    if (p == 1) then
      call copy_groups(groups, result%groups)
      allocate(result%merge_correlation(0))
      result%status = ren_success
      return
    end if
    correlation = correlation_matrix(x, istat)
    do i = 1, p
      correlation(i, i) = 0.0_dp
    end do
    allocate(result%merge_correlation(p - 1))
    best_rho = huge(1.0_dp)
    do step = 1, p - 1
      call maximum_off_diagonal(correlation, rho, first, second)
      result%merge_correlation(step) = rho
      if (rho < best_rho) then
        best_rho = rho
        call copy_groups(groups, best_groups)
      end if
      if (size(groups) == 1) exit
      allocate(next_groups(size(groups) - 1))
      allocate(next_groups(1)%index(size(groups(first)%index) + size(groups(second)%index)))
      next_groups(1)%index = [groups(first)%index, groups(second)%index]
      out_index = 1
      do i = 1, size(groups)
        if (i == first .or. i == second) cycle
        out_index = out_index + 1
        allocate(next_groups(out_index)%index(size(groups(i)%index)))
        next_groups(out_index)%index = groups(i)%index
      end do
      if (size(next_groups) > 1) then
        allocate(next_correlation(size(next_groups), size(next_groups)))
        next_correlation = 0.0_dp
        do i = 2, size(next_groups)
          next_correlation(1, i) = maximum_canonical_correlation(x(:, next_groups(1)%index), &
            x(:, next_groups(i)%index))
          next_correlation(i, 1) = next_correlation(1, i)
        end do
        do j = 2, size(next_groups)
          do i = 2, size(next_groups)
            next_correlation(i, j) = correlation(old_group_position(groups, next_groups(i)), &
              old_group_position(groups, next_groups(j)))
          end do
        end do
        do i = 1, size(next_groups)
          next_correlation(i, i) = 0.0_dp
        end do
      else
        allocate(next_correlation(1, 1))
        next_correlation = 0.0_dp
      end if
      call move_alloc(next_groups, groups)
      call move_alloc(next_correlation, correlation)
    end do
    call copy_groups(best_groups, result%groups)
    result%status = ren_success
  end function buh_clust

  subroutine po_bhu(y0, x0, groups, repetitions, weights, status, variance_tolerance, seed)
    real(dp), intent(in) :: y0(:), x0(:, :)
    type(asset_group), intent(in) :: groups(:)
    integer, intent(in) :: repetitions
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    integer, intent(in), optional :: seed
    real(dp), allocatable :: x(:, :), selected_x(:, :), covariance(:, :), solution(:), portfolio(:), &
      design(:, :), response(:), beta(:), trial(:), mean_weight(:)
    integer, allocatable :: active(:), selected(:)
    integer :: p, g, r, j, k, istat, aggregate_status, use_seed
    call filter_columns(x0, x, active, variance_tolerance)
    p = size(x, 2)
    g = size(groups)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    if (size(y0) /= size(x0, 1) .or. repetitions < 1 .or. g < 1) then
      status = ren_invalid_argument
      return
    end if
    if (p == 0) then
      weights = 1.0_dp / real(size(weights), dp)
      status = ren_success
      return
    end if
    allocate(mean_weight(p), selected(g))
    mean_weight = 0.0_dp
    use_seed = 3001
    if (present(seed)) use_seed = seed
    call initialize_random_seed(use_seed)
    aggregate_status = ren_success
    do r = 1, repetitions
      do j = 1, g
        selected(j) = random_choice(groups(j)%index)
      end do
      selected_x = x(:, selected)
      covariance = covariance_matrix(selected_x, istat)
      call long_only_minimum_variance(covariance, solution, istat)
      allocate(trial(p))
      trial = 0.0_dp
      k = maxloc(solution, dim=1)
      if (abs(solution(k) - 1.0_dp) <= 5.0e-6_dp) then
        k = selected(k)
        allocate(design(size(x, 1), p - 1), response(size(x, 1)))
        response = x(:, k) - y0
        j = 0
        do istat = 1, p
          if (istat == k) cycle
          j = j + 1
          design(:, j) = x(:, k) - x(:, istat)
        end do
        call fit_regularized_cv(design, response, 1.0_dp, beta, istat, seed=use_seed + r)
        trial = insert_at(beta, k - 1, 1.0_dp - sum(beta))
      else
        portfolio = matmul(selected_x, solution)
        allocate(design(size(x, 1), p), response(size(x, 1)))
        response = portfolio - y0
        do j = 1, p
          design(:, j) = portfolio - x(:, j)
        end do
        call fit_regularized_cv(design, response, 1.0_dp, beta, istat, seed=use_seed + r)
        trial = beta
        trial(selected) = trial(selected) + (1.0_dp - sum(beta)) * solution
      end if
      if (istat /= ren_success) aggregate_status = istat
      mean_weight = mean_weight + trial
      deallocate(trial)
      if (allocated(design)) deallocate(design)
      if (allocated(response)) deallocate(response)
      if (allocated(beta)) deallocate(beta)
    end do
    mean_weight = mean_weight / real(repetitions, dp)
    weights(active) = mean_weight
    status = aggregate_status
  end subroutine po_bhu

  subroutine po_tzt(x0, gamma, weights, status, variance_tolerance)
    real(dp), intent(in) :: x0(:, :), gamma
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    real(dp), allocatable :: x(:, :), covariance(:, :), mw(:), equal(:), mean_x(:), inverse(:, :), active_weights(:)
    real(dp) :: theta_sq, c1, pi1, pi2, delta
    integer, allocatable :: active(:)
    integer :: p, n, istat
    call filter_columns(x0, x, active, variance_tolerance)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    p = size(x, 2)
    n = size(x, 1)
    if (p == 0 .or. abs(gamma) <= epsilon(1.0_dp)) then
      if (size(weights) > 0) weights = 1.0_dp / real(size(weights), dp)
      status = ren_invalid_argument
      return
    end if
    covariance = covariance_matrix(x, istat)
    covariance = make_positive_definite(covariance)
    call equality_minimum_variance(covariance, mw, istat)
    allocate(equal(p))
    equal = 1.0_dp / real(p, dp)
    mean_x = sum(x, dim=1) / real(n, dp)
    inverse = pseudoinverse(covariance)
    theta_sq = dot_product(mean_x, matmul(inverse, mean_x))
    if ((n - p - 1) == 0 .or. (n - p - 4) == 0) then
      active_weights = equal
      status = ren_numerical_error
    else
      c1 = real((n - 2) * (n - p - 2), dp) / real((n - p - 1) * (n - p - 4), dp)
      pi1 = dot_product(equal, matmul(covariance, equal)) - 2.0_dp / gamma * dot_product(equal, mean_x) + &
        theta_sq / (gamma * gamma)
      pi2 = (c1 - 1.0_dp) * theta_sq / (gamma * gamma) + c1 * real(p, dp) / &
        (gamma * gamma * real(n, dp))
      if (abs(pi1 + pi2) <= epsilon(1.0_dp)) then
        delta = 0.0_dp
      else
        delta = pi1 / (pi1 + pi2)
      end if
      active_weights = (1.0_dp - delta) * equal + delta * mw
      status = ren_success
    end if
    weights(active) = active_weights
  end subroutine po_tzt

  subroutine po_sw(x0, b, samples, weights, status, variance_tolerance, seed)
    real(dp), intent(in) :: x0(:, :)
    integer, intent(in) :: b, samples
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    integer, intent(in), optional :: seed
    real(dp), allocatable :: x(:, :), covariance(:, :), selected_weights(:), mean_weight(:)
    integer, allocatable :: active(:), selected(:)
    integer :: p, use_b, r, istat, use_seed
    call filter_columns(x0, x, active, variance_tolerance)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    p = size(x, 2)
    if (p == 0 .or. samples < 1) then
      if (size(weights) > 0) weights = 1.0_dp / real(size(weights), dp)
      status = ren_invalid_argument
      return
    end if
    use_b = max(1, min(p, b))
    allocate(selected(use_b), mean_weight(p))
    mean_weight = 0.0_dp
    use_seed = 4001
    if (present(seed)) use_seed = seed
    call initialize_random_seed(use_seed)
    status = ren_success
    do r = 1, samples
      call sample_without_replacement(p, use_b, selected)
      covariance = covariance_matrix(x(:, selected), istat)
      call equality_minimum_variance(covariance, selected_weights, istat)
      mean_weight(selected) = mean_weight(selected) + selected_weights
      if (istat /= ren_success) status = istat
    end do
    mean_weight = mean_weight / real(samples, dp)
    weights(active) = mean_weight
  end subroutine po_sw

  subroutine po_sw_lasso(y0, x0, b, samples, weights, status, variance_tolerance, seed)
    real(dp), intent(in) :: y0(:), x0(:, :)
    integer, intent(in) :: b, samples
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: variance_tolerance
    integer, intent(in), optional :: seed
    real(dp), allocatable :: x(:, :), selected_x(:, :), covariance(:, :), solution(:), portfolio(:), &
      design(:, :), response(:), beta(:), trial(:), mean_weight(:)
    integer, allocatable :: active(:), selected(:)
    integer :: p, use_b, r, j, k, istat, aggregate_status, use_seed
    call filter_columns(x0, x, active, variance_tolerance)
    allocate(weights(size(x0, 2)))
    weights = 0.0_dp
    p = size(x, 2)
    if (size(y0) /= size(x0, 1) .or. p == 0 .or. samples < 1) then
      if (p == 0 .and. size(weights) > 0) weights = 1.0_dp / real(size(weights), dp)
      status = ren_invalid_argument
      return
    end if
    use_b = max(1, min(p, b))
    allocate(selected(use_b), mean_weight(p))
    mean_weight = 0.0_dp
    use_seed = 5003
    if (present(seed)) use_seed = seed
    call initialize_random_seed(use_seed)
    aggregate_status = ren_success
    do r = 1, samples
      call sample_without_replacement(p, use_b, selected)
      selected_x = x(:, selected)
      covariance = covariance_matrix(selected_x, istat)
      call long_only_minimum_variance(covariance, solution, istat)
      allocate(trial(p))
      trial = 0.0_dp
      k = maxloc(solution, dim=1)
      if (abs(solution(k) - 1.0_dp) <= 5.0e-6_dp) then
        k = selected(k)
        allocate(design(size(x, 1), p - 1), response(size(x, 1)))
        response = x(:, k) - y0
        j = 0
        do istat = 1, p
          if (istat == k) cycle
          j = j + 1
          design(:, j) = x(:, k) - x(:, istat)
        end do
        call fit_regularized_cv(design, response, 1.0_dp, beta, istat, seed=use_seed + r)
        trial = insert_at(beta, k - 1, 1.0_dp - sum(beta))
      else
        portfolio = matmul(selected_x, solution)
        allocate(design(size(x, 1), p), response(size(x, 1)))
        response = portfolio - y0
        do j = 1, p
          design(:, j) = portfolio - x(:, j)
        end do
        call fit_regularized_cv(design, response, 1.0_dp, beta, istat, seed=use_seed + r)
        trial = beta
        trial(selected) = trial(selected) + (1.0_dp - sum(beta)) * solution
      end if
      if (istat /= ren_success) aggregate_status = istat
      mean_weight = mean_weight + trial
      deallocate(trial)
      if (allocated(design)) deallocate(design)
      if (allocated(response)) deallocate(response)
      if (allocated(beta)) deallocate(beta)
    end do
    mean_weight = mean_weight / real(samples, dp)
    weights(active) = mean_weight
    status = aggregate_status
  end subroutine po_sw_lasso

  subroutine filter_columns(x0, x, active, variance_tolerance)
    real(dp), intent(in) :: x0(:, :)
    real(dp), allocatable, intent(out) :: x(:, :)
    integer, allocatable, intent(out) :: active(:)
    real(dp), intent(in), optional :: variance_tolerance
    real(dp), allocatable :: variance(:)
    real(dp) :: tolerance
    integer :: j, k
    tolerance = 0.0_dp
    if (present(variance_tolerance)) tolerance = max(0.0_dp, variance_tolerance)
    allocate(variance(size(x0, 2)))
    variance = column_variances(x0)
    active = pack([(j, j=1,size(x0, 2))], variance > tolerance)
    allocate(x(size(x0, 1), size(active)))
    do k = 1, size(active)
      x(:, k) = x0(:, active(k))
    end do
  end subroutine filter_columns

  subroutine maximum_off_diagonal(matrix, value, first, second)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(out) :: value
    integer, intent(out) :: first, second
    integer :: i, j
    value = -huge(1.0_dp)
    first = 1
    second = min(2, size(matrix, 1))
    do j = 2, size(matrix, 2)
      do i = 1, j - 1
        if (matrix(i, j) > value) then
          value = matrix(i, j)
          first = i
          second = j
        end if
      end do
    end do
  end subroutine maximum_off_diagonal

  function maximum_canonical_correlation(a, b) result(value)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp) :: value
    real(dp), allocatable :: ac(:, :), bc(:, :), cross(:, :)
    type(svd_result) :: sa, sb, sc
    if (size(a, 1) /= size(b, 1) .or. size(a, 1) < 2) then
      value = 0.0_dp
      return
    end if
    ac = a - spread(sum(a, dim=1) / real(size(a, 1), dp), 1, size(a, 1))
    bc = b - spread(sum(b, dim=1) / real(size(b, 1), dp), 1, size(b, 1))
    sa = fast_svd(ac)
    sb = fast_svd(bc)
    if (sa%rank == 0 .or. sb%rank == 0) then
      value = 0.0_dp
      return
    end if
    cross = matmul(transpose(sa%u), sb%u)
    sc = fast_svd(cross)
    if (sc%rank == 0) then
      value = 0.0_dp
    else
      value = max(0.0_dp, min(1.0_dp, maxval(sc%d)))
    end if
  end function maximum_canonical_correlation

  subroutine copy_groups(source, destination)
    type(asset_group), intent(in) :: source(:)
    type(asset_group), allocatable, intent(out) :: destination(:)
    integer :: i
    if (allocated(destination)) deallocate(destination)
    allocate(destination(size(source)))
    do i = 1, size(source)
      allocate(destination(i)%index(size(source(i)%index)))
      destination(i)%index = source(i)%index
    end do
  end subroutine copy_groups

  integer function old_group_position(old_groups, new_group) result(position)
    type(asset_group), intent(in) :: old_groups(:), new_group
    integer :: i
    position = 1
    do i = 1, size(old_groups)
      if (size(old_groups(i)%index) /= size(new_group%index)) cycle
      if (all(old_groups(i)%index == new_group%index)) then
        position = i
        return
      end if
    end do
  end function old_group_position
end module ren_portfolio
