! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_pca
  use rrcov_kinds, only : dp
  use rrcov_types, only : pca_result, covariance_result, rrcov_success, rrcov_invalid_argument
  use rrcov_random, only : seed_random, random_unit_vector
  use rrcov_linalg, only : symmetric_eigen, make_positive_definite, orthonormalize
  use rrcov_stats, only : mean_vector, covariance_matrix, spatial_median, qn_scale, &
    mad_scale, standardize_columns, robust_standardize_columns
  use rrcov_robust, only : cov_mcd
  implicit none
  private
  public :: pca_classic, pca_cov, pca_locantore, pca_grid, pca_proj, pca_hubert
  public :: pca_distances
contains
  subroutine pca_classic(x, result, k, center_data, scale_data, signflip)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: k
    logical, intent(in), optional :: center_data, scale_data, signflip
    real(dp), allocatable :: z(:, :), covariance(:, :), values(:), vectors(:, :)
    logical :: center_flag, scale_flag, flip
    integer :: components, status
    center_flag = .true.
    if (present(center_data)) center_flag = center_data
    scale_flag = .false.
    if (present(scale_data)) scale_flag = scale_data
    flip = .true.
    if (present(signflip)) flip = signflip
    call prepare_pca_data(x, z, result%center, result%scale, status, center_flag, scale_flag, .false.)
    if (status /= rrcov_success) then
      result%status = status
      result%method = "Classical PCA"
      return
    end if
    covariance = covariance_matrix(z, unbiased=.true., status=status)
    call symmetric_eigen(covariance, values, vectors, status)
    components = min(size(x, 2), size(x, 1))
    if (present(k)) components = min(components, max(1, k))
    result%loadings = vectors(:, 1:components)
    if (flip) call flip_loading_signs(result%loadings)
    result%eigenvalues = max(values(1:components), 0.0_dp)
    result%scores = matmul(z, result%loadings)
    call pca_distances(x, result)
    result%n_obs = size(x, 1)
    result%n_components = components
    result%status = status
    result%method = "Classical PCA"
  end subroutine pca_classic

  subroutine pca_cov(x, covariance_estimate, result, k, scale_data, signflip)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(in) :: covariance_estimate
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: k
    logical, intent(in), optional :: scale_data, signflip
    real(dp), allocatable :: covariance(:, :), values(:), vectors(:, :), z(:, :)
    logical :: scale_flag, flip
    integer :: components, j, status, n, p
    n = size(x, 1)
    p = size(x, 2)
    if (size(covariance_estimate%center) /= p .or. &
        size(covariance_estimate%covariance, 1) /= p) then
      result%status = rrcov_invalid_argument
      result%method = "Covariance-based PCA"
      return
    end if
    scale_flag = .false.
    if (present(scale_data)) scale_flag = scale_data
    flip = .true.
    if (present(signflip)) flip = signflip
    allocate(result%center(p), result%scale(p), z(n, p))
    result%center = covariance_estimate%center
    result%scale = 1.0_dp
    covariance = covariance_estimate%covariance
    if (scale_flag) then
      do j = 1, p
        result%scale(j) = sqrt(max(covariance(j, j), tiny(1.0_dp)))
        z(:, j) = (x(:, j) - result%center(j)) / result%scale(j)
      end do
      do j = 1, p
        covariance(j, :) = covariance(j, :) / result%scale(j)
        covariance(:, j) = covariance(:, j) / result%scale(j)
      end do
    else
      z = x - spread(result%center, 1, n)
    end if
    covariance = make_positive_definite(covariance, 1.0e-10_dp)
    call symmetric_eigen(covariance, values, vectors, status)
    components = min(p, n)
    if (present(k)) components = min(components, max(1, k))
    result%loadings = vectors(:, 1:components)
    if (flip) call flip_loading_signs(result%loadings)
    result%eigenvalues = max(values(1:components), 0.0_dp)
    result%scores = matmul(z, result%loadings)
    call pca_distances(x, result)
    result%n_obs = n
    result%n_components = components
    result%status = status
    result%method = "Covariance-based robust PCA"
  end subroutine pca_cov

  subroutine pca_locantore(x, result, k, scale_data, signflip)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: k
    logical, intent(in), optional :: scale_data, signflip
    real(dp), allocatable :: z(:, :), directions(:, :), scatter(:, :), values(:), vectors(:, :)
    real(dp) :: norm_value
    logical :: scale_flag, flip
    integer :: n, p, i, j, components, status
    n = size(x, 1)
    p = size(x, 2)
    scale_flag = .false.
    if (present(scale_data)) scale_flag = scale_data
    flip = .true.
    if (present(signflip)) flip = signflip
    if (scale_flag) then
      call robust_standardize_columns(x, z, result%center, result%scale, status)
      call spatial_median(z, values, status)
      result%center = result%center + result%scale * values
      z = z - spread(values, 1, n)
    else
      call spatial_median(x, result%center, status)
      allocate(result%scale(p), z(n, p))
      result%scale = 1.0_dp
      z = x - spread(result%center, 1, n)
    end if
    allocate(directions(n, p))
    directions = 0.0_dp
    do i = 1, n
      norm_value = sqrt(sum(z(i, :) ** 2))
      if (norm_value > tiny(1.0_dp)) directions(i, :) = z(i, :) / norm_value
    end do
    scatter = matmul(transpose(directions), directions) / real(max(1, n), dp)
    call symmetric_eigen(scatter, values, vectors, status)
    components = min(p, n)
    if (present(k)) components = min(components, max(1, k))
    result%loadings = vectors(:, 1:components)
    if (flip) call flip_loading_signs(result%loadings)
    result%scores = matmul(z, result%loadings)
    allocate(result%eigenvalues(components))
    do j = 1, components
      result%eigenvalues(j) = qn_scale(result%scores(:, j)) ** 2
      if (result%eigenvalues(j) <= tiny(1.0_dp)) &
        result%eigenvalues(j) = mad_scale(result%scores(:, j), consistency=.true.) ** 2
    end do
    call pca_distances(x, result)
    result%n_obs = n
    result%n_components = components
    result%status = status
    result%method = "Locantore spherical PCA"
  end subroutine pca_locantore

  subroutine pca_grid(x, result, k, max_directions, seed, scale_data, signflip)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: k, max_directions, seed
    logical, intent(in), optional :: scale_data, signflip
    real(dp), allocatable :: z(:, :), candidate(:), best(:), projection(:), loadings(:, :)
    real(dp) :: objective, best_objective, norm_value
    logical :: scale_flag, flip
    integer :: n, p, components, ndir, component, d, status, rank
    n = size(x, 1)
    p = size(x, 2)
    components = min(p, n)
    if (present(k)) components = min(components, max(1, k))
    ndir = max(250, 50 * p)
    if (present(max_directions)) ndir = max(p, max_directions)
    scale_flag = .false.
    if (present(scale_data)) scale_flag = scale_data
    flip = .true.
    if (present(signflip)) flip = signflip
    call prepare_pca_data(x, z, result%center, result%scale, status, .true., scale_flag, .true.)
    call seed_random(seed)
    allocate(loadings(p, components), candidate(p), best(p), projection(n))
    loadings = 0.0_dp
    allocate(result%eigenvalues(components))
    do component = 1, components
      best_objective = -1.0_dp
      best = 0.0_dp
      do d = 1, ndir
        if (d <= p) then
          candidate = 0.0_dp
          candidate(d) = 1.0_dp
        else
          call random_unit_vector(candidate)
        end if
        if (component > 1) then
          candidate = candidate - matmul(loadings(:, 1:component - 1), &
            matmul(transpose(loadings(:, 1:component - 1)), candidate))
        end if
        norm_value = sqrt(sum(candidate * candidate))
        if (norm_value <= 1.0e-10_dp) cycle
        candidate = candidate / norm_value
        projection = matmul(z, candidate)
        objective = qn_scale(projection)
        if (objective > best_objective) then
          best_objective = objective
          best = candidate
        end if
      end do
      loadings(:, component) = best
      call orthonormalize(loadings(:, 1:component), rank)
      result%eigenvalues(component) = max(best_objective, 0.0_dp) ** 2
    end do
    result%loadings = loadings
    if (flip) call flip_loading_signs(result%loadings)
    result%scores = matmul(z, result%loadings)
    call pca_distances(x, result)
    result%n_obs = n
    result%n_components = components
    result%status = rrcov_success
    result%method = "Projection-pursuit grid PCA"
  end subroutine pca_grid

  subroutine pca_proj(x, result, k, max_directions, seed, scale_data, signflip)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: k, max_directions, seed
    logical, intent(in), optional :: scale_data, signflip
    call pca_grid(x, result, k=k, max_directions=max_directions, seed=seed, &
      scale_data=scale_data, signflip=signflip)
    result%method = "Projection-pursuit PCA"
  end subroutine pca_proj

  subroutine pca_hubert(x, result, k, kmax, alpha, nsamp, seed, scale_data, signflip)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: k, kmax, nsamp, seed
    real(dp), intent(in), optional :: alpha
    logical, intent(in), optional :: scale_data, signflip
    type(pca_result) :: initial
    type(covariance_result) :: robust_scores
    real(dp), allocatable :: values(:), vectors(:,:), center_standardized(:), z(:, :)
    integer :: components, initial_components, status, n, p
    logical :: scale_flag, flip
    n = size(x, 1)
    p = size(x, 2)
    components = min(min(n, p), 10)
    if (present(k)) components = min(min(n, p), max(1, k))
    initial_components = min(min(n, p), max(components, 10))
    if (present(kmax)) initial_components = min(min(n, p), max(components, kmax))
    scale_flag = .false.
    if (present(scale_data)) scale_flag = scale_data
    flip = .true.
    if (present(signflip)) flip = signflip
    call pca_classic(x, initial, k=initial_components, scale_data=scale_flag, signflip=.false.)
    call cov_mcd(initial%scores(:, 1:initial_components), robust_scores, alpha=alpha, &
      nsamp=nsamp, seed=seed)
    call symmetric_eigen(robust_scores%covariance, values, vectors, status)
    allocate(result%center(p), result%scale(p))
    result%scale = initial%scale
    center_standardized = matmul(initial%loadings, robust_scores%center)
    result%center = initial%center + initial%scale * center_standardized
    result%loadings = matmul(initial%loadings, vectors(:, 1:components))
    if (flip) call flip_loading_signs(result%loadings)
    result%eigenvalues = max(values(1:components), 0.0_dp)
    allocate(z(n, p))
    z = (x - spread(result%center, 1, n)) / spread(result%scale, 1, n)
    result%scores = matmul(z, result%loadings)
    call pca_distances(x, result)
    result%n_obs = n
    result%n_components = components
    result%status = robust_scores%status
    result%method = "Hubert robust PCA"
  end subroutine pca_hubert

  subroutine pca_distances(x, result)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(inout) :: result
    real(dp), allocatable :: z(:, :), reconstructed(:, :), residual(:, :)
    integer :: i, j, n, p, k
    n = size(x, 1)
    p = size(x, 2)
    k = size(result%loadings, 2)
    allocate(z(n, p), reconstructed(n, p), residual(n, p))
    z = (x - spread(result%center, 1, n)) / spread(result%scale, 1, n)
    if (.not. allocated(result%scores)) result%scores = matmul(z, result%loadings)
    reconstructed = matmul(result%scores, transpose(result%loadings))
    residual = z - reconstructed
    allocate(result%score_distances(n), result%orthogonal_distances(n))
    result%score_distances = 0.0_dp
    do j = 1, k
      if (result%eigenvalues(j) > tiny(1.0_dp)) then
        result%score_distances = result%score_distances + &
          result%scores(:, j) ** 2 / result%eigenvalues(j)
      end if
    end do
    do i = 1, n
      result%orthogonal_distances(i) = sqrt(sum(residual(i, :) ** 2))
    end do
  end subroutine pca_distances

  subroutine prepare_pca_data(x, z, center, scale, status, center_data, scale_data, robust)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: z(:, :), center(:), scale(:)
    integer, intent(out) :: status
    logical, intent(in) :: center_data, scale_data, robust
    real(dp), allocatable :: covariance(:, :)
    integer :: j, n, p
    n = size(x, 1)
    p = size(x, 2)
    if (n < 2 .or. p < 1) then
      allocate(z(n, p), center(p), scale(p))
      z = 0.0_dp
      center = 0.0_dp
      scale = 1.0_dp
      status = rrcov_invalid_argument
      return
    end if
    if (robust .and. scale_data) then
      call robust_standardize_columns(x, z, center, scale, status)
      return
    end if
    allocate(z(n, p), center(p), scale(p))
    if (center_data) then
      if (robust) then
        center = 0.0_dp
        do j = 1, p
          center(j) = median_local(x(:, j))
        end do
      else
        center = mean_vector(x)
      end if
    else
      center = 0.0_dp
    end if
    scale = 1.0_dp
    if (scale_data) then
      covariance = covariance_matrix(x, unbiased=.true., status=status)
      do j = 1, p
        scale(j) = sqrt(max(covariance(j, j), tiny(1.0_dp)))
      end do
    end if
    z = (x - spread(center, 1, n)) / spread(scale, 1, n)
    status = rrcov_success
  end subroutine prepare_pca_data

  subroutine flip_loading_signs(loadings)
    real(dp), intent(inout) :: loadings(:, :)
    integer :: j, index_max
    do j = 1, size(loadings, 2)
      index_max = maxloc(abs(loadings(:, j)), dim=1)
      if (loadings(index_max, j) < 0.0_dp) loadings(:, j) = -loadings(:, j)
    end do
  end subroutine flip_loading_signs

  function median_local(x) result(value)
    use rrcov_stats, only : median
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = median(x)
  end function median_local
end module rrcov_pca
