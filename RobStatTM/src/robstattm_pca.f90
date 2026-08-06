! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_pca
  use robstattm_kinds, only : dp
  use robstattm_types, only : pca_result, robstattm_success, robstattm_invalid_argument, &
    robstattm_no_convergence
  use robstattm_psi, only : scale_m
  use rrcov_types, only : rrcov_pca_result => pca_result
  use rrcov_pca, only : pca_locantore
  use rrcov_linalg, only : symmetric_eigen
  implicit none
  private
  public :: pca_rob_s, sm_pca, prcomp_rob
contains
  subroutine pca_rob_s(x, result, ncomp, desired_proportion, delta_scale, max_iter, tolerance)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: ncomp, max_iter
    real(dp), intent(in), optional :: desired_proportion, delta_scale, tolerance
    type(rrcov_pca_result) :: initial
    real(dp), allocatable :: x_centered(:, :), fit_centered(:, :), residuals(:, :)
    real(dp), allocatable :: residual_norm2(:), weights(:), covariance(:, :)
    real(dp), allocatable :: values(:), vectors(:, :), old_loadings(:, :), new_center(:)
    real(dp) :: desired, delta, tol, initial_total_scale, old_scale, new_scale
    real(dp) :: scale_change, loading_change, denominator
    integer :: n, p, qmax, q, q_by_variance, maxit, status, i, iter

    n = size(x, 1)
    p = size(x, 2)
    if (n < 2 .or. p < 1) then
      result%status = robstattm_invalid_argument
      result%method = 'Residual M-scale PCA'
      return
    end if

    qmax = min(ncomp_default(ncomp, p), min(n, p))
    desired = 0.9_dp
    if (present(desired_proportion)) desired = min(1.0_dp, max(0.0_dp, desired_proportion))
    delta = 0.5_dp
    if (present(delta_scale)) delta = min(0.999999_dp, max(1.0e-6_dp, delta_scale))
    maxit = 100
    if (present(max_iter)) maxit = max(1, max_iter)
    tol = 1.0e-4_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))

    call pca_locantore(x, initial, k=min(n, p), scale_data=.false., signflip=.true.)
    if (initial%status /= 0 .or. .not. allocated(initial%loadings) .or. &
        .not. allocated(initial%eigenvalues)) then
      result%status = robstattm_invalid_argument
      result%method = 'Residual M-scale PCA'
      return
    end if

    allocate(result%initial_cumulative_variance(size(initial%eigenvalues)))
    denominator = sum(max(initial%eigenvalues, 0.0_dp))
    if (denominator > tiny(1.0_dp)) then
      result%initial_cumulative_variance = cumulative_sum(max(initial%eigenvalues, 0.0_dp) / denominator)
    else
      result%initial_cumulative_variance = 0.0_dp
    end if
    q_by_variance = 1
    do i = 1, size(result%initial_cumulative_variance)
      q_by_variance = i
      if (result%initial_cumulative_variance(i) >= desired) exit
    end do
    q = min(qmax, q_by_variance)
    q = max(1, q)

    allocate(result%center(p), result%loadings(p, q), result%scores(n, q), &
      result%fitted(n, p), result%sdev(q))
    allocate(x_centered(n, p), fit_centered(n, p), residuals(n, p), residual_norm2(n), &
      weights(n), covariance(p, p), old_loadings(p, q), new_center(p))

    result%center = initial%center
    result%loadings = initial%loadings(:, 1:q)
    old_loadings = result%loadings
    x_centered = x - spread(result%center, 1, n)
    fit_centered = matmul(matmul(x_centered, result%loadings), transpose(result%loadings))
    residuals = x_centered - fit_centered
    residual_norm2 = sum(residuals * residuals, dim=2)
    initial_total_scale = scale_m(sqrt(sum(x_centered * x_centered, dim=2)), delta, &
      'bisquare', 1.0_dp) ** 2
    old_scale = scale_m(sqrt(max(residual_norm2, 0.0_dp)), delta, 'bisquare', 1.0_dp) ** 2
    old_scale = max(old_scale, sqrt(tiny(1.0_dp)))
    weights = residual_weights(residual_norm2 / old_scale)

    result%converged = .false.
    do iter = 1, maxit
      result%iterations = iter
      denominator = sum(weights)
      if (denominator <= tiny(1.0_dp)) exit
      new_center = matmul(weights, x) / denominator
      x_centered = x - spread(new_center, 1, n)
      covariance = 0.0_dp
      do i = 1, n
        covariance = covariance + weights(i) * outer_row(x_centered(i, :))
      end do
      covariance = covariance / denominator
      call symmetric_eigen(covariance, values, vectors, status)
      if (status /= 0) exit
      result%loadings = vectors(:, 1:q)
      call align_loading_signs(result%loadings, old_loadings)
      fit_centered = matmul(matmul(x_centered, result%loadings), transpose(result%loadings))
      residuals = x_centered - fit_centered
      residual_norm2 = sum(residuals * residuals, dim=2)
      new_scale = scale_m(sqrt(max(residual_norm2, 0.0_dp)), delta, 'bisquare', 1.0_dp) ** 2
      new_scale = max(new_scale, sqrt(tiny(1.0_dp)))
      scale_change = abs(1.0_dp - new_scale / old_scale)
      loading_change = mean_abs_identity_minus_overlap(result%loadings, old_loadings)
      result%center = new_center
      if (max(scale_change, loading_change) <= tol) then
        result%converged = .true.
        old_scale = new_scale
        exit
      end if
      old_scale = new_scale
      old_loadings = result%loadings
      weights = residual_weights(residual_norm2 / old_scale)
    end do

    x_centered = x - spread(result%center, 1, n)
    result%scores = matmul(x_centered, result%loadings)
    result%fitted = matmul(result%scores, transpose(result%loadings)) + &
      spread(result%center, 1, n)
    do i = 1, q
      result%sdev(i) = sample_standard_deviation(result%scores(:, i))
    end do
    if (initial_total_scale > tiny(1.0_dp)) then
      result%explained_proportion = max(0.0_dp, min(1.0_dp, 1.0_dp - old_scale / initial_total_scale))
    else
      result%explained_proportion = 0.0_dp
    end if
    result%n_components = q
    result%method = 'Residual M-scale PCA'
    if (result%converged) then
      result%status = robstattm_success
    else
      result%status = robstattm_no_convergence
    end if
  end subroutine pca_rob_s

  subroutine sm_pca(x, result, ncomp, desired_proportion, delta_scale, max_iter, tolerance)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: ncomp, max_iter
    real(dp), intent(in), optional :: desired_proportion, delta_scale, tolerance
    call pca_rob_s(x, result, ncomp, desired_proportion, delta_scale, max_iter, tolerance)
  end subroutine sm_pca

  subroutine prcomp_rob(x, result, rank, delta_scale, max_iter, tolerance)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: rank, max_iter
    real(dp), intent(in), optional :: delta_scale, tolerance
    integer :: components
    components = size(x, 2)
    if (present(rank)) components = max(1, min(rank, size(x, 2)))
    call pca_rob_s(x, result, ncomp=components, desired_proportion=1.0_dp, &
      delta_scale=delta_scale, max_iter=max_iter, tolerance=tolerance)
    result%method = 'prcompRob-compatible PCA'
  end subroutine prcomp_rob

  pure integer function ncomp_default(ncomp, p) result(value)
    integer, intent(in), optional :: ncomp
    integer, intent(in) :: p
    value = p
    if (present(ncomp)) value = max(1, min(ncomp, p))
  end function ncomp_default

  pure function cumulative_sum(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value(size(x))
    integer :: i
    if (size(x) == 0) return
    value(1) = x(1)
    do i = 2, size(x)
      value(i) = value(i - 1) + x(i)
    end do
  end function cumulative_sum

  pure function residual_weights(ratio) result(weights)
    real(dp), intent(in) :: ratio(:)
    real(dp) :: weights(size(ratio))
    where (ratio >= 0.0_dp .and. ratio <= 1.0_dp)
      weights = (1.0_dp - ratio) ** 2
    elsewhere
      weights = 0.0_dp
    end where
  end function residual_weights

  pure function outer_row(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value(size(x), size(x))
    integer :: i
    do i = 1, size(x)
      value(i, :) = x(i) * x
    end do
  end function outer_row

  pure subroutine align_loading_signs(loadings, reference)
    real(dp), intent(inout) :: loadings(:, :)
    real(dp), intent(in) :: reference(:, :)
    integer :: j
    do j = 1, min(size(loadings, 2), size(reference, 2))
      if (dot_product(loadings(:, j), reference(:, j)) < 0.0_dp) &
        loadings(:, j) = -loadings(:, j)
    end do
  end subroutine align_loading_signs

  pure function mean_abs_identity_minus_overlap(loadings, reference) result(value)
    real(dp), intent(in) :: loadings(:, :), reference(:, :)
    real(dp) :: value
    real(dp), allocatable :: overlap(:, :)
    integer :: i, q
    q = min(size(loadings, 2), size(reference, 2))
    allocate(overlap(q, q))
    overlap = abs(matmul(transpose(loadings(:, 1:q)), reference(:, 1:q)))
    do i = 1, q
      overlap(i, i) = 1.0_dp - overlap(i, i)
    end do
    do i = 1, q
      if (i > 1) overlap(i, 1:i-1) = -overlap(i, 1:i-1)
      if (i < q) overlap(i, i+1:q) = -overlap(i, i+1:q)
    end do
    value = sum(abs(overlap)) / real(q * q, dp)
  end function mean_abs_identity_minus_overlap

  pure function sample_standard_deviation(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, mu
    if (size(x) < 2) then
      value = 0.0_dp
      return
    end if
    mu = sum(x) / real(size(x), dp)
    value = sqrt(sum((x - mu) ** 2) / real(size(x) - 1, dp))
  end function sample_standard_deviation
end module robstattm_pca
