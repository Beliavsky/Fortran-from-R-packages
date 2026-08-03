! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_robust
  use rrcov_kinds, only : dp
  use rrcov_types, only : covariance_result, rrcov_success, rrcov_invalid_argument, &
    rrcov_dimension_error, rrcov_no_convergence
  use rrcov_random, only : seed_random, random_subset, random_unit_vector
  use rrcov_sort, only : order_smallest
  use rrcov_linalg, only : symmetric_eigen, make_positive_definite, &
    mahalanobis_squared, log_determinant, identity_matrix, outer_product
  use rrcov_stats, only : mean_vector, covariance_matrix, weighted_mean_covariance, &
    median, mad_scale, qn_scale, tau_scale, chi_square_quantile, column_medians
  implicit none
  private
  public :: cov_classic, cov_mcd, cov_mve, cov_ogk, cov_mest
  public :: cov_sest, cov_mmest, cov_sde, cov_mrcd
  public :: robust_covariance, adjusted_outlyingness, medcouple
  public :: h_alpha_n, cov_mwcd
contains
  pure function h_alpha_n(alpha, n, p) result(h)
    real(dp), intent(in) :: alpha
    integer, intent(in) :: n, p
    integer :: h, n2
    n2 = floor(0.5_dp * real(n + p + 1, dp))
    h = floor(2.0_dp * real(n2, dp) - real(n, dp) + &
      2.0_dp * real(n - n2, dp) * min(1.0_dp, max(0.5_dp, alpha)))
    h = min(n, max(p + 1, h))
  end function h_alpha_n

  subroutine cov_classic(x, result, unbiased)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    logical, intent(in), optional :: unbiased
    logical :: ub
    integer :: status, i
    ub = .true.
    if (present(unbiased)) ub = unbiased
    result%n_obs = size(x, 1)
    result%rank = size(x, 2)
    result%center = mean_vector(x)
    result%covariance = covariance_matrix(x, unbiased=ub, status=status)
    call mahalanobis_squared(x, result%center, result%covariance, result%distances, status)
    allocate(result%weights(size(x, 1)), result%subset(size(x, 1)))
    result%weights = 1.0_dp
    result%subset = [(i, i=1, size(x, 1))]
    result%status = status
    result%method = "Classical covariance"
    result%objective = log_determinant(make_positive_definite(result%covariance))
  end subroutine cov_classic

  subroutine cov_mcd(x, result, alpha, nsamp, seed, max_csteps, reweight)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: nsamp, seed, max_csteps
    logical, intent(in), optional :: reweight
    real(dp), allocatable :: center(:), covariance(:, :), distances(:), &
      best_center(:), best_covariance(:, :), work_cov(:, :), weights(:)
    integer, allocatable :: initial(:), subset(:), best_subset(:)
    real(dp) :: a, objective, best_objective, factor, cutoff, q50
    integer :: n, p, h, trials, maxsteps, trial, status, iterations, i
    logical :: do_reweight

    n = size(x, 1)
    p = size(x, 2)
    if (n <= p .or. p < 1) then
      call empty_covariance_result(result, n, p, rrcov_invalid_argument, "Minimum covariance determinant")
      return
    end if
    a = 0.5_dp
    if (present(alpha)) a = min(1.0_dp, max(0.5_dp, alpha))
    h = h_alpha_n(a, n, p)
    trials = min(500, max(50, 10 * p))
    if (present(nsamp)) trials = max(1, nsamp)
    maxsteps = 50
    if (present(max_csteps)) maxsteps = max(1, max_csteps)
    do_reweight = .true.
    if (present(reweight)) do_reweight = reweight
    call seed_random(seed)
    allocate(initial(min(p + 1, h)), subset(h), best_subset(h))
    allocate(best_center(p), best_covariance(p, p))
    best_objective = huge(1.0_dp)
    iterations = 0

    do trial = 1, trials
      if (trial == 1) then
        initial = [(i, i=1, size(initial))]
      else
        call random_subset(n, size(initial), initial)
      end if
      call subset_covariance(x, initial, center, covariance, status)
      if (status /= rrcov_success) cycle
      covariance = make_positive_definite(covariance, 1.0e-7_dp)
      call c_step(x, h, center, covariance, subset, maxsteps, status, iterations)
      if (status /= rrcov_success .and. status /= rrcov_no_convergence) cycle
      objective = log_determinant(make_positive_definite(covariance, 1.0e-10_dp), status)
      if (objective < best_objective) then
        best_objective = objective
        best_center = center
        best_covariance = covariance
        best_subset = subset
      end if
    end do
    if (best_objective >= 0.5_dp * huge(1.0_dp)) then
      call cov_classic(x, result)
      result%status = rrcov_no_convergence
      result%method = "Minimum covariance determinant (fallback)"
      return
    end if

    call mahalanobis_squared(x, best_center, best_covariance, distances, status)
    q50 = chi_square_quantile(0.5_dp, real(p, dp))
    factor = median(distances) / max(q50, tiny(1.0_dp))
    if (factor > tiny(1.0_dp)) best_covariance = factor * best_covariance
    call mahalanobis_squared(x, best_center, best_covariance, distances, status)
    allocate(weights(n))
    weights = 0.0_dp
    if (do_reweight) then
      cutoff = chi_square_quantile(0.975_dp, real(p, dp))
      where (distances <= cutoff) weights = 1.0_dp
      if (sum(weights) > real(p + 1, dp)) then
        call weighted_mean_covariance(x, weights, center, work_cov, status, normalize=.false.)
        if (status == rrcov_success) then
          best_center = center
          best_covariance = make_positive_definite(work_cov, 1.0e-10_dp)
          call mahalanobis_squared(x, best_center, best_covariance, distances, status)
        end if
      end if
    else
      weights(best_subset) = 1.0_dp
    end if
    result%center = best_center
    result%covariance = best_covariance
    result%distances = distances
    result%weights = weights
    result%subset = best_subset
    result%n_obs = n
    result%rank = p
    result%iterations = iterations
    result%status = rrcov_success
    result%objective = best_objective
    result%method = "Minimum covariance determinant"
  end subroutine cov_mcd

  subroutine cov_mve(x, result, alpha, nsamp, seed, reweight)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: nsamp, seed
    logical, intent(in), optional :: reweight
    real(dp), allocatable :: center(:), covariance(:, :), distances(:), &
      best_center(:), best_covariance(:, :), weights(:), reweighted(:, :)
    integer, allocatable :: initial(:), subset(:), best_subset(:)
    real(dp) :: a, objective, best_objective, radius, cutoff, factor
    integer :: n, p, h, trials, trial, status
    logical :: do_reweight
    n = size(x, 1)
    p = size(x, 2)
    if (n <= p .or. p < 1) then
      call empty_covariance_result(result, n, p, rrcov_invalid_argument, "Minimum volume ellipsoid")
      return
    end if
    a = 0.5_dp
    if (present(alpha)) a = min(1.0_dp, max(0.5_dp, alpha))
    h = h_alpha_n(a, n, p)
    trials = min(1000, max(100, 20 * p))
    if (present(nsamp)) trials = max(1, nsamp)
    do_reweight = .true.
    if (present(reweight)) do_reweight = reweight
    call seed_random(seed)
    allocate(initial(p + 1), subset(h), best_subset(h), best_center(p), best_covariance(p, p))
    best_objective = huge(1.0_dp)
    do trial = 1, trials
      call random_subset(n, p + 1, initial)
      call subset_covariance(x, initial, center, covariance, status)
      if (status /= rrcov_success) cycle
      covariance = make_positive_definite(covariance, 1.0e-6_dp)
      call mahalanobis_squared(x, center, covariance, distances, status)
      call order_smallest(distances, h, subset)
      radius = maxval(distances(subset))
      call subset_covariance(x, subset, center, covariance, status)
      covariance = make_positive_definite(covariance, 1.0e-8_dp)
      objective = log_determinant(covariance, status) + real(p, dp) * log(max(radius, tiny(1.0_dp)))
      if (objective < best_objective) then
        best_objective = objective
        best_center = center
        best_covariance = covariance
        best_subset = subset
      end if
    end do
    if (best_objective >= 0.5_dp * huge(1.0_dp)) then
      call cov_classic(x, result)
      result%status = rrcov_no_convergence
      result%method = "Minimum volume ellipsoid (fallback)"
      return
    end if
    call mahalanobis_squared(x, best_center, best_covariance, distances, status)
    factor = median(distances) / max(chi_square_quantile(0.5_dp, real(p, dp)), tiny(1.0_dp))
    if (factor > tiny(1.0_dp)) best_covariance = factor * best_covariance
    call mahalanobis_squared(x, best_center, best_covariance, distances, status)
    allocate(weights(n))
    weights = 0.0_dp
    if (do_reweight) then
      cutoff = chi_square_quantile(0.975_dp, real(p, dp))
      where (distances <= cutoff) weights = 1.0_dp
      if (sum(weights) > real(p + 1, dp)) then
        call weighted_mean_covariance(x, weights, center, reweighted, status, normalize=.false.)
        best_center = center
        best_covariance = make_positive_definite(reweighted, 1.0e-10_dp)
        call mahalanobis_squared(x, best_center, best_covariance, distances, status)
      end if
    else
      weights(best_subset) = 1.0_dp
    end if
    result%center = best_center
    result%covariance = best_covariance
    result%distances = distances
    result%weights = weights
    result%subset = best_subset
    result%n_obs = n
    result%rank = p
    result%iterations = trials
    result%status = rrcov_success
    result%objective = best_objective
    result%method = "Minimum volume ellipsoid"
  end subroutine cov_mve

  subroutine cov_ogk(x, result, niter, use_tau, use_quadrant)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    integer, intent(in), optional :: niter
    logical, intent(in), optional :: use_tau, use_quadrant
    real(dp), allocatable :: z(:, :), y(:, :), location(:), scales(:), &
      correlation(:, :), values(:), vectors(:, :), transform(:, :), &
      old_transform(:, :), final_location(:), final_scales(:), distances(:)
    real(dp) :: pair_value
    integer :: n, p, iterations, iteration, i, j, status
    logical :: tau, quadrant
    n = size(x, 1)
    p = size(x, 2)
    if (n < 2 .or. p < 1) then
      call empty_covariance_result(result, n, p, rrcov_invalid_argument, "Orthogonalized Gnanadesikan-Kettenring")
      return
    end if
    iterations = 2
    if (present(niter)) iterations = max(1, niter)
    tau = .false.
    if (present(use_tau)) tau = use_tau
    quadrant = .false.
    if (present(use_quadrant)) quadrant = use_quadrant
    allocate(z(n, p), transform(p, p))
    z = x
    transform = identity_matrix(p)
    allocate(result%center(p))
    result%center = 0.0_dp
    do iteration = 1, iterations
      call marginal_location_scale(z, location, scales, tau)
      allocate(y(n, p), correlation(p, p))
      do j = 1, p
        y(:, j) = (z(:, j) - location(j)) / scales(j)
      end do
      correlation = identity_matrix(p)
      do j = 2, p
        do i = 1, j - 1
          if (quadrant) then
            pair_value = quadrant_correlation(y(:, i), y(:, j))
          else
            pair_value = gk_covariance(y(:, i), y(:, j), tau)
          end if
          correlation(i, j) = pair_value
          correlation(j, i) = pair_value
        end do
      end do
      correlation = make_positive_definite(correlation, 1.0e-8_dp)
      call symmetric_eigen(correlation, values, vectors, status)
      result%center = result%center + matmul(transform, location)
      old_transform = transform
      transform = matmul(old_transform, matmul(diagonal_matrix(scales), vectors))
      z = matmul(y, vectors)
      deallocate(y, correlation, location, scales, values, vectors)
    end do
    call marginal_location_scale(z, final_location, final_scales, tau)
    result%center = result%center + matmul(transform, final_location)
    result%covariance = matmul(transform, matmul(diagonal_matrix(final_scales ** 2), transpose(transform)))
    result%covariance = make_positive_definite(result%covariance, 1.0e-10_dp)
    call mahalanobis_squared(x, result%center, result%covariance, distances, status)
    result%distances = distances
    allocate(result%weights(n), result%subset(n))
    result%weights = 1.0_dp
    result%subset = [(i, i=1, n)]
    result%n_obs = n
    result%rank = p
    result%iterations = iterations
    result%status = rrcov_success
    result%objective = log_determinant(result%covariance)
    result%method = "Orthogonalized Gnanadesikan-Kettenring"
  end subroutine cov_ogk

  subroutine cov_mest(x, result, arp, tolerance, max_iterations)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: arp, tolerance
    integer, intent(in), optional :: max_iterations
    type(covariance_result) :: initial
    real(dp), allocatable :: center(:), covariance(:, :), distances(:), weights(:), &
      next_center(:), next_covariance(:, :)
    real(dp) :: alpha, cutoff, tol, error
    integer :: n, p, iteration, maxit, status
    n = size(x, 1)
    p = size(x, 2)
    alpha = 0.05_dp
    if (present(arp)) alpha = min(0.49_dp, max(1.0e-6_dp, arp))
    tol = 1.0e-6_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 120
    if (present(max_iterations)) maxit = max(1, max_iterations)
    call cov_ogk(x, initial)
    center = initial%center
    covariance = initial%covariance
    allocate(weights(n))
    cutoff = sqrt(chi_square_quantile(1.0_dp - alpha, real(p, dp)))
    do iteration = 1, maxit
      call mahalanobis_squared(x, center, covariance, distances, status)
      weights = 1.0_dp
      where (sqrt(distances) > cutoff)
        weights = cutoff / max(sqrt(distances), tiny(1.0_dp))
      end where
      call weighted_mean_covariance(x, weights, next_center, next_covariance, status, normalize=.true.)
      next_covariance = make_positive_definite(next_covariance, 1.0e-9_dp)
      error = max(maxval(abs(next_center - center)), maxval(abs(next_covariance - covariance)))
      center = next_center
      covariance = next_covariance
      if (error <= tol * max(1.0_dp, maxval(abs(covariance)))) exit
    end do
    call finish_iterative_result(x, center, covariance, weights, iteration, maxit, result, &
      "Multivariate M-estimator")
  end subroutine cov_mest

  subroutine cov_sest(x, result, bdp, tolerance, max_iterations, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: bdp, tolerance
    integer, intent(in), optional :: max_iterations, nsamp, seed
    type(covariance_result) :: initial
    real(dp), allocatable :: center(:), covariance(:, :), distances(:), weights(:), &
      next_center(:), next_covariance(:, :)
    real(dp) :: breakdown, cutoff, tol, error, u, correction
    integer :: n, p, i, iteration, maxit, samples
    n = size(x, 1)
    p = size(x, 2)
    breakdown = 0.5_dp
    if (present(bdp)) breakdown = min(0.5_dp, max(0.05_dp, bdp))
    tol = 1.0e-6_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 150
    if (present(max_iterations)) maxit = max(1, max_iterations)
    samples = min(500, max(50, 10 * p))
    if (present(nsamp)) samples = max(1, nsamp)
    call cov_mcd(x, initial, alpha=max(0.5_dp, 1.0_dp - breakdown), nsamp=samples, seed=seed)
    center = initial%center
    covariance = initial%covariance
    allocate(weights(n))
    cutoff = sqrt(chi_square_quantile(min(0.999_dp, 1.0_dp - 0.5_dp * breakdown), real(p, dp)))
    cutoff = max(cutoff, sqrt(real(p, dp)) + 1.5_dp)
    do iteration = 1, maxit
      call mahalanobis_squared(x, center, covariance, distances, i)
      do i = 1, n
        u = sqrt(distances(i)) / cutoff
        if (u < 1.0_dp) then
          weights(i) = (1.0_dp - u * u) ** 2
        else
          weights(i) = 0.0_dp
        end if
      end do
      if (sum(weights) <= real(p + 1, dp)) exit
      call weighted_mean_covariance(x, weights, next_center, next_covariance, i, normalize=.true.)
      next_covariance = make_positive_definite(next_covariance, 1.0e-9_dp)
      call mahalanobis_squared(x, next_center, next_covariance, distances, i)
      correction = median(distances) / max(chi_square_quantile(0.5_dp, real(p, dp)), tiny(1.0_dp))
      if (correction > tiny(1.0_dp)) next_covariance = correction * next_covariance
      error = max(maxval(abs(next_center - center)), maxval(abs(next_covariance - covariance)))
      center = next_center
      covariance = next_covariance
      if (error <= tol * max(1.0_dp, maxval(abs(covariance)))) exit
    end do
    call finish_iterative_result(x, center, covariance, weights, iteration, maxit, result, &
      "Multivariate S-estimator")
  end subroutine cov_sest

  subroutine cov_mmest(x, result, bdp, efficiency, tolerance, max_iterations, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: bdp, efficiency, tolerance
    integer, intent(in), optional :: max_iterations, nsamp, seed
    type(covariance_result) :: initial
    real(dp), allocatable :: center(:), covariance(:, :), distances(:), weights(:), &
      next_center(:), next_covariance(:, :)
    real(dp) :: eff, tol, cutoff, u, error
    integer :: n, p, i, iteration, maxit
    n = size(x, 1)
    p = size(x, 2)
    eff = 0.95_dp
    if (present(efficiency)) eff = min(0.999_dp, max(0.5_dp, efficiency))
    tol = 1.0e-6_dp
    if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
    maxit = 100
    if (present(max_iterations)) maxit = max(1, max_iterations)
    call cov_sest(x, initial, bdp=bdp, nsamp=nsamp, seed=seed)
    center = initial%center
    covariance = initial%covariance
    allocate(weights(n))
    cutoff = max(sqrt(chi_square_quantile(eff, real(p, dp))), sqrt(real(p, dp)) + 2.5_dp)
    do iteration = 1, maxit
      call mahalanobis_squared(x, center, covariance, distances, i)
      do i = 1, n
        u = sqrt(distances(i)) / cutoff
        if (u < 1.0_dp) then
          weights(i) = (1.0_dp - u * u) ** 2
        else
          weights(i) = 0.0_dp
        end if
      end do
      if (sum(weights) <= real(p + 1, dp)) exit
      call weighted_mean_covariance(x, weights, next_center, next_covariance, i, normalize=.true.)
      next_covariance = make_positive_definite(next_covariance, 1.0e-9_dp)
      error = max(maxval(abs(next_center - center)), maxval(abs(next_covariance - covariance)))
      center = next_center
      covariance = next_covariance
      if (error <= tol * max(1.0_dp, maxval(abs(covariance)))) exit
    end do
    call finish_iterative_result(x, center, covariance, weights, iteration, maxit, result, &
      "Multivariate MM-estimator")
  end subroutine cov_mmest

  subroutine cov_sde(x, result, alpha, ndirections, seed, reweight)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: ndirections, seed
    logical, intent(in), optional :: reweight
    real(dp), allocatable :: center0(:), direction(:), projection(:), outlyingness(:), &
      center(:), covariance(:, :), distances(:), weights(:)
    integer, allocatable :: subset(:)
    real(dp) :: a, loc, scale, cutoff
    integer :: n, p, h, ndir, d, i, j, status
    logical :: do_reweight
    n = size(x, 1)
    p = size(x, 2)
    if (n <= p .or. p < 1) then
      call empty_covariance_result(result, n, p, rrcov_invalid_argument, "Stahel-Donoho estimator")
      return
    end if
    a = 0.5_dp
    if (present(alpha)) a = min(1.0_dp, max(0.5_dp, alpha))
    h = h_alpha_n(a, n, p)
    ndir = max(250, 25 * p)
    if (present(ndirections)) ndir = max(p, ndirections)
    do_reweight = .true.
    if (present(reweight)) do_reweight = reweight
    call seed_random(seed)
    center0 = column_medians(x)
    allocate(direction(p), projection(n), outlyingness(n))
    outlyingness = 0.0_dp
    do d = 1, ndir
      if (d <= p) then
        direction = 0.0_dp
        direction(d) = 1.0_dp
      else if (d <= p + min(n - 1, ndir - p)) then
        i = 1 + modulo(37 * d, n)
        j = 1 + modulo(97 * d + 11, n)
        if (i == j) j = modulo(j, n) + 1
        direction = x(i, :) - x(j, :)
        scale = sqrt(sum(direction * direction))
        if (scale <= tiny(1.0_dp)) call random_unit_vector(direction)
        if (scale > tiny(1.0_dp)) direction = direction / scale
      else
        call random_unit_vector(direction)
      end if
      projection = matmul(x, direction)
      loc = median(projection)
      scale = mad_scale(projection, consistency=.true.)
      if (scale <= tiny(1.0_dp)) cycle
      outlyingness = max(outlyingness, abs(projection - loc) / scale)
    end do
    allocate(subset(h))
    call order_smallest(outlyingness, h, subset)
    call subset_covariance(x, subset, center, covariance, status)
    covariance = make_positive_definite(covariance, 1.0e-9_dp)
    call mahalanobis_squared(x, center, covariance, distances, status)
    allocate(weights(n))
    weights = 0.0_dp
    if (do_reweight) then
      cutoff = chi_square_quantile(0.975_dp, real(p, dp))
      where (distances <= cutoff) weights = 1.0_dp
      if (sum(weights) > real(p + 1, dp)) then
        call weighted_mean_covariance(x, weights, center, covariance, status, normalize=.false.)
        covariance = make_positive_definite(covariance, 1.0e-9_dp)
        call mahalanobis_squared(x, center, covariance, distances, status)
      end if
    else
      weights(subset) = 1.0_dp
    end if
    result%center = center
    result%covariance = covariance
    result%distances = distances
    result%weights = weights
    result%subset = subset
    result%n_obs = n
    result%rank = p
    result%iterations = ndir
    result%status = rrcov_success
    result%objective = log_determinant(covariance)
    result%method = "Stahel-Donoho estimator"
  end subroutine cov_sde

  subroutine cov_mrcd(x, result, alpha, rho, nsamp, seed, max_csteps)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha, rho
    integer, intent(in), optional :: nsamp, seed, max_csteps
    real(dp), allocatable :: center(:), covariance(:, :), regularized(:, :), target(:, :), &
      distances(:), best_center(:), best_covariance(:, :), best_regularized(:, :), scales(:)
    integer, allocatable :: initial(:), subset(:), best_subset(:)
    real(dp) :: a, r, objective, best_objective, condition_limit
    integer :: n, p, h, trials, trial, status, maxsteps, iterations, j
    n = size(x, 1)
    p = size(x, 2)
    if (n <= p .or. p < 1) then
      call empty_covariance_result(result, n, p, rrcov_invalid_argument, &
        "Minimum regularized covariance determinant")
      return
    end if
    a = 0.75_dp
    if (present(alpha)) a = min(1.0_dp, max(0.5_dp, alpha))
    r = -1.0_dp
    if (present(rho)) r = min(1.0_dp, max(0.0_dp, rho))
    h = h_alpha_n(a, n, p)
    trials = min(500, max(50, 10 * p))
    if (present(nsamp)) trials = max(1, nsamp)
    maxsteps = 50
    if (present(max_csteps)) maxsteps = max(1, max_csteps)
    call seed_random(seed)
    allocate(scales(p), target(p, p))
    target = 0.0_dp
    do j = 1, p
      scales(j) = qn_scale(x(:, j))
      if (scales(j) <= tiny(1.0_dp)) scales(j) = 1.0_dp
      target(j, j) = scales(j) ** 2
    end do
    allocate(initial(min(p + 1, h)), subset(h), best_subset(h))
    allocate(best_center(p), best_covariance(p, p), best_regularized(p, p))
    best_objective = huge(1.0_dp)
    condition_limit = 50.0_dp
    iterations = 0
    do trial = 1, trials
      call random_subset(n, size(initial), initial)
      call subset_covariance(x, initial, center, covariance, status)
      if (status /= rrcov_success) cycle
      call regularize_covariance(covariance, target, r, condition_limit, regularized)
      call c_step_regularized(x, h, center, covariance, target, r, condition_limit, &
        subset, maxsteps, status, iterations)
      call regularize_covariance(covariance, target, r, condition_limit, regularized)
      objective = log_determinant(regularized, status)
      if (objective < best_objective) then
        best_objective = objective
        best_center = center
        best_covariance = covariance
        best_regularized = regularized
        best_subset = subset
      end if
    end do
    if (best_objective >= 0.5_dp * huge(1.0_dp)) then
      call cov_ogk(x, result)
      result%status = rrcov_no_convergence
      result%method = "MRCD (fallback OGK)"
      return
    end if
    call mahalanobis_squared(x, best_center, best_regularized, distances, status)
    allocate(result%weights(n))
    result%weights = 0.0_dp
    result%weights(best_subset) = 1.0_dp
    result%center = best_center
    result%covariance = best_regularized
    result%distances = distances
    result%subset = best_subset
    result%n_obs = n
    result%rank = p
    result%iterations = iterations
    result%status = rrcov_success
    result%objective = best_objective
    result%method = "Minimum regularized covariance determinant"
  end subroutine cov_mrcd

  subroutine robust_covariance(x, method, result, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    character(len=*), intent(in) :: method
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: nsamp, seed
    character(len=:), allocatable :: key
    key = lowercase(trim(adjustl(method)))
    select case (key)
    case ("classic", "classical", "cov")
      call cov_classic(x, result)
    case ("mcd")
      call cov_mcd(x, result, alpha=alpha, nsamp=nsamp, seed=seed)
    case ("mve")
      call cov_mve(x, result, alpha=alpha, nsamp=nsamp, seed=seed)
    case ("ogk")
      call cov_ogk(x, result)
    case ("mest", "m")
      call cov_mest(x, result)
    case ("sest", "s")
      call cov_sest(x, result, nsamp=nsamp, seed=seed)
    case ("mmest", "mm")
      call cov_mmest(x, result, nsamp=nsamp, seed=seed)
    case ("sde", "donostah")
      call cov_sde(x, result, alpha=alpha, seed=seed)
    case ("mrcd")
      call cov_mrcd(x, result, alpha=alpha, nsamp=nsamp, seed=seed)
    case default
      call cov_classic(x, result)
      result%status = rrcov_invalid_argument
      result%method = "Unknown method; classical fallback"
    end select
  end subroutine robust_covariance

  subroutine cov_mwcd(x, grouping, result, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(covariance_result), intent(out) :: result
    real(dp), intent(in), optional :: alpha
    integer, intent(in), optional :: nsamp, seed
    real(dp), allocatable :: residuals(:, :), group_data(:, :)
    integer, allocatable :: labels(:), rows(:)
    type(covariance_result) :: group_result
    integer :: n, p, g, ng, i, count_rows
    n = size(x, 1)
    p = size(x, 2)
    if (size(grouping) /= n) then
      call empty_covariance_result(result, n, p, rrcov_dimension_error, "Multi-group MCD")
      return
    end if
    call unique_labels(grouping, labels)
    ng = size(labels)
    allocate(residuals(n, p))
    do g = 1, ng
      count_rows = count(grouping == labels(g))
      allocate(rows(count_rows), group_data(count_rows, p))
      rows = pack([(i, i=1, n)], grouping == labels(g))
      group_data = x(rows, :)
      call cov_mcd(group_data, group_result, alpha=alpha, nsamp=nsamp, seed=seed)
      if (group_result%status == rrcov_invalid_argument) then
        group_result%center = mean_vector(group_data)
      end if
      do i = 1, count_rows
        residuals(rows(i), :) = group_data(i, :) - group_result%center
      end do
      deallocate(rows, group_data)
    end do
    call cov_mcd(residuals, result, alpha=alpha, nsamp=nsamp, seed=seed)
    result%method = "Multi-group minimum covariance determinant"
  end subroutine cov_mwcd

  subroutine adjusted_outlyingness(x, outlyingness, ndirections, seed)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: outlyingness(:)
    integer, intent(in), optional :: ndirections, seed
    real(dp), allocatable :: direction(:), projection(:)
    real(dp) :: med, mc, lower_scale, upper_scale, value
    integer :: n, p, ndir, d, i
    n = size(x, 1)
    p = size(x, 2)
    ndir = max(250, 25 * p)
    if (present(ndirections)) ndir = max(p, ndirections)
    call seed_random(seed)
    allocate(outlyingness(n), direction(p), projection(n))
    outlyingness = 0.0_dp
    do d = 1, ndir
      if (d <= p) then
        direction = 0.0_dp
        direction(d) = 1.0_dp
      else
        call random_unit_vector(direction)
      end if
      projection = matmul(x, direction)
      med = median(projection)
      mc = medcouple(projection)
      lower_scale = max(med - minval(projection), mad_scale(projection, consistency=.true.))
      upper_scale = max(maxval(projection) - med, mad_scale(projection, consistency=.true.))
      if (mc >= 0.0_dp) then
        lower_scale = lower_scale * exp(-4.0_dp * mc)
        upper_scale = upper_scale * exp(3.0_dp * mc)
      else
        lower_scale = lower_scale * exp(-3.0_dp * mc)
        upper_scale = upper_scale * exp(4.0_dp * mc)
      end if
      do i = 1, n
        if (projection(i) >= med) then
          value = (projection(i) - med) / max(upper_scale, tiny(1.0_dp))
        else
          value = (med - projection(i)) / max(lower_scale, tiny(1.0_dp))
        end if
        outlyingness(i) = max(outlyingness(i), value)
      end do
    end do
  end subroutine adjusted_outlyingness

  function medcouple(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, med, denominator
    real(dp), allocatable :: lower(:), upper(:), kernels(:)
    integer :: i, j, k, nl, nu
    med = median(x)
    lower = pack(x, x <= med)
    upper = pack(x, x >= med)
    nl = size(lower)
    nu = size(upper)
    allocate(kernels(nl * nu))
    k = 0
    do i = 1, nl
      do j = 1, nu
        k = k + 1
        denominator = upper(j) - lower(i)
        if (abs(denominator) <= tiny(1.0_dp)) then
          kernels(k) = 0.0_dp
        else
          kernels(k) = (upper(j) + lower(i) - 2.0_dp * med) / denominator
        end if
      end do
    end do
    value = median(kernels)
  end function medcouple

  subroutine c_step(x, h, center, covariance, subset, maxsteps, status, iterations)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: h, maxsteps
    real(dp), allocatable, intent(inout) :: center(:), covariance(:, :)
    integer, intent(out) :: subset(h)
    integer, intent(out) :: status, iterations
    real(dp), allocatable :: distances(:), next_center(:), next_covariance(:, :)
    integer, allocatable :: previous(:)
    integer :: step
    allocate(previous(h))
    previous = 0
    status = rrcov_success
    do step = 1, maxsteps
      call mahalanobis_squared(x, center, covariance, distances, status)
      call order_smallest(distances, h, subset)
      if (all(subset == previous)) exit
      previous = subset
      call subset_covariance(x, subset, next_center, next_covariance, status)
      if (status /= rrcov_success) return
      center = next_center
      covariance = make_positive_definite(next_covariance, 1.0e-9_dp)
    end do
    iterations = step
    if (step > maxsteps) status = rrcov_no_convergence
  end subroutine c_step

  subroutine c_step_regularized(x, h, center, covariance, target, rho, condition_limit, &
      subset, maxsteps, status, iterations)
    real(dp), intent(in) :: x(:, :), target(:, :), rho, condition_limit
    integer, intent(in) :: h, maxsteps
    real(dp), allocatable, intent(inout) :: center(:), covariance(:, :)
    integer, intent(out) :: subset(h)
    integer, intent(out) :: status, iterations
    real(dp), allocatable :: distances(:), regularized(:, :), next_center(:), next_covariance(:, :)
    integer, allocatable :: previous(:)
    integer :: step
    allocate(previous(h))
    previous = 0
    status = rrcov_success
    do step = 1, maxsteps
      call regularize_covariance(covariance, target, rho, condition_limit, regularized)
      call mahalanobis_squared(x, center, regularized, distances, status)
      call order_smallest(distances, h, subset)
      if (all(subset == previous)) exit
      previous = subset
      call subset_covariance(x, subset, next_center, next_covariance, status)
      center = next_center
      covariance = next_covariance
    end do
    iterations = step
    if (step > maxsteps) status = rrcov_no_convergence
  end subroutine c_step_regularized

  subroutine regularize_covariance(covariance, target, requested_rho, condition_limit, regularized)
    real(dp), intent(in) :: covariance(:, :), target(:, :), requested_rho, condition_limit
    real(dp), allocatable, intent(out) :: regularized(:, :)
    real(dp), allocatable :: values(:), vectors(:, :)
    real(dp) :: rho_value, condition
    integer :: status, trial
    if (requested_rho >= 0.0_dp) then
      rho_value = requested_rho
    else
      rho_value = 0.0_dp
    end if
    do trial = 1, 101
      regularized = (1.0_dp - rho_value) * covariance + rho_value * target
      regularized = make_positive_definite(regularized, 1.0e-10_dp)
      call symmetric_eigen(regularized, values, vectors, status)
      condition = maxval(values) / max(minval(values), tiny(1.0_dp))
      if (condition <= condition_limit .or. requested_rho >= 0.0_dp) exit
      rho_value = min(1.0_dp, rho_value + 0.01_dp)
    end do
  end subroutine regularize_covariance

  subroutine subset_covariance(x, subset, center, covariance, status)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: subset(:)
    real(dp), allocatable, intent(out) :: center(:), covariance(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: data(:, :)
    if (size(subset) < 2 .or. any(subset < 1) .or. any(subset > size(x, 1))) then
      allocate(center(size(x, 2)), covariance(size(x, 2), size(x, 2)))
      center = 0.0_dp
      covariance = 0.0_dp
      status = rrcov_invalid_argument
      return
    end if
    data = x(subset, :)
    center = mean_vector(data)
    covariance = covariance_matrix(data, unbiased=.false., status=status)
  end subroutine subset_covariance

  subroutine finish_iterative_result(x, center, covariance, weights, iteration, maxit, result, method)
    real(dp), intent(in) :: x(:, :), center(:), covariance(:, :), weights(:)
    integer, intent(in) :: iteration, maxit
    type(covariance_result), intent(out) :: result
    character(len=*), intent(in) :: method
    integer :: status, i, count_positive
    call mahalanobis_squared(x, center, covariance, result%distances, status)
    result%center = center
    result%covariance = covariance
    result%weights = weights
    count_positive = count(weights > 0.0_dp)
    allocate(result%subset(count_positive))
    result%subset = pack([(i, i=1, size(weights))], weights > 0.0_dp)
    result%n_obs = size(x, 1)
    result%rank = size(x, 2)
    result%iterations = min(iteration, maxit)
    if (iteration > maxit) then
      result%status = rrcov_no_convergence
    else
      result%status = status
    end if
    result%objective = log_determinant(covariance)
    result%method = method
  end subroutine finish_iterative_result

  subroutine marginal_location_scale(x, location, scale, tau)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: location(:), scale(:)
    logical, intent(in) :: tau
    integer :: j, p
    p = size(x, 2)
    allocate(location(p), scale(p))
    do j = 1, p
      if (tau) then
        scale(j) = tau_scale(x(:, j), location(j))
      else
        location(j) = median(x(:, j))
        scale(j) = mad_scale(x(:, j), consistency=.true.)
      end if
      if (scale(j) <= tiny(1.0_dp)) scale(j) = 1.0_dp
    end do
  end subroutine marginal_location_scale

  function gk_covariance(x, y, tau) result(value)
    real(dp), intent(in) :: x(:), y(:)
    logical, intent(in) :: tau
    real(dp) :: value, plus, minus
    if (tau) then
      plus = tau_scale(x + y)
      minus = tau_scale(x - y)
    else
      plus = mad_scale(x + y, consistency=.true.)
      minus = mad_scale(x - y, consistency=.true.)
    end if
    value = 0.25_dp * (plus * plus - minus * minus)
    value = min(1.0_dp, max(-1.0_dp, value))
  end function gk_covariance

  function quadrant_correlation(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: value, mx, my
    integer :: i, concordant, discordant
    mx = median(x)
    my = median(y)
    concordant = 0
    discordant = 0
    do i = 1, size(x)
      if ((x(i) > mx .and. y(i) > my) .or. (x(i) < mx .and. y(i) < my)) then
        concordant = concordant + 1
      else if ((x(i) > mx .and. y(i) < my) .or. (x(i) < mx .and. y(i) > my)) then
        discordant = discordant + 1
      end if
    end do
    if (concordant + discordant == 0) then
      value = 0.0_dp
    else
      value = sin(0.5_dp * acos(-1.0_dp) * real(concordant - discordant, dp) / &
        real(concordant + discordant, dp))
    end if
  end function quadrant_correlation

  function diagonal_matrix(diagonal) result(value)
    real(dp), intent(in) :: diagonal(:)
    real(dp) :: value(size(diagonal), size(diagonal))
    integer :: i
    value = 0.0_dp
    do i = 1, size(diagonal)
      value(i, i) = diagonal(i)
    end do
  end function diagonal_matrix

  subroutine empty_covariance_result(result, n, p, status, method)
    type(covariance_result), intent(out) :: result
    integer, intent(in) :: n, p, status
    character(len=*), intent(in) :: method
    allocate(result%center(max(0, p)), result%covariance(max(0, p), max(0, p)))
    allocate(result%distances(max(0, n)), result%weights(max(0, n)), result%subset(0))
    result%center = 0.0_dp
    result%covariance = 0.0_dp
    result%distances = 0.0_dp
    result%weights = 0.0_dp
    result%n_obs = n
    result%rank = 0
    result%status = status
    result%method = method
  end subroutine empty_covariance_result

  subroutine unique_labels(grouping, labels)
    integer, intent(in) :: grouping(:)
    integer, allocatable, intent(out) :: labels(:)
    integer, allocatable :: work(:)
    integer :: i, count_labels
    allocate(work(size(grouping)))
    count_labels = 0
    do i = 1, size(grouping)
      if (count_labels == 0 .or. .not. any(work(1:count_labels) == grouping(i))) then
        count_labels = count_labels + 1
        work(count_labels) = grouping(i)
      end if
    end do
    allocate(labels(count_labels))
    labels = work(1:count_labels)
  end subroutine unique_labels

  pure function lowercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: i, code
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        value(i:i) = achar(code + 32)
      else
        value(i:i) = text(i:i)
      end if
    end do
  end function lowercase
end module rrcov_robust
