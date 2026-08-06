! SPDX-License-Identifier: GPL-3.0-only
module mass_mds
  use rrcov_kinds, only : dp
  use rrcov_types, only : rrcov_success
  use rrcov_linalg, only : symmetric_eigen
  use mass_types, only : mds_result, mass_success, mass_invalid_argument, &
    mass_dimension_error, mass_no_convergence
  use mass_math, only : isotonic_increasing, sort_real_mass
  implicit none
  private
  public :: classical_mds, sammon, iso_mds, shepard
contains

  subroutine classical_mds(distance, dimensions, points, status)
    real(dp), intent(in) :: distance(:, :)
    integer, intent(in) :: dimensions
    real(dp), allocatable, intent(out) :: points(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: d2(:, :), centering(:, :), gram(:, :), &
      values(:), vectors(:, :)
    integer :: n, i, st, k

    n = size(distance, 1)
    if (size(distance, 2) /= n .or. dimensions < 1 .or. dimensions >= n) then
      allocate(points(0, 0))
      status = mass_dimension_error
      return
    end if
    allocate(d2(n, n), centering(n, n))
    d2 = distance * distance
    centering = -1.0_dp / real(n, dp)
    do i = 1, n
      centering(i, i) = centering(i, i) + 1.0_dp
    end do
    gram = -0.5_dp * matmul(centering, matmul(d2, centering))
    call symmetric_eigen(gram, values, vectors, st)
    k = min(dimensions, count(values > 0.0_dp))
    allocate(points(n, dimensions))
    points = 0.0_dp
    do i = 1, k
      points(:, i) = vectors(:, i) * sqrt(values(i))
    end do
    status = merge(mass_success, mass_no_convergence, st == rrcov_success)
  end subroutine classical_mds

  subroutine sammon(distance, result, initial, dimensions, maxit, tolerance, magic)
    real(dp), intent(in) :: distance(:, :)
    type(mds_result), intent(out) :: result
    real(dp), intent(in), optional :: initial(:, :), tolerance, magic
    integer, intent(in), optional :: dimensions, maxit
    real(dp), allocatable :: y(:, :), candidate(:, :), gradient(:), curvature(:)
    real(dp) :: total, stress, previous, checkpoint, step, tol, target, fitted
    integer :: n, k, mit, iter, i, j, m, status

    n = size(distance, 1)
    k = 2
    if (present(dimensions)) k = dimensions
    mit = 100
    if (present(maxit)) mit = maxit
    tol = 1.0e-4_dp
    if (present(tolerance)) tol = tolerance
    step = 0.2_dp
    if (present(magic)) step = magic
    if (size(distance, 2) /= n .or. k < 1 .or. n <= k .or. &
        any(distance < 0.0_dp)) then
      result%status = mass_invalid_argument
      return
    end if
    if (present(initial)) then
      if (size(initial, 1) /= n .or. size(initial, 2) /= k) then
        result%status = mass_dimension_error
        return
      end if
      y = initial
    else
      call classical_mds(distance, k, y, status)
      if (status /= mass_success) then
        result%status = status
        return
      end if
    end if
    allocate(candidate(n, k), gradient(k), curvature(k))
    call sammon_stress(distance, y, stress, total)
    previous = stress
    checkpoint = stress
    do iter = 1, mit
      do
        do i = 1, n
          gradient = 0.0_dp
          curvature = 0.0_dp
          do j = 1, n
            if (i == j .or. distance(i, j) <= tiny(1.0_dp)) cycle
            fitted = sqrt(sum((y(i, :) - y(j, :))**2))
            if (fitted <= tiny(1.0_dp)) then
              result%status = mass_invalid_argument
              return
            end if
            target = distance(i, j) - fitted
            do m = 1, k
              gradient(m) = gradient(m) + &
                (y(i, m) - y(j, m)) * target / (distance(i, j) * fitted)
              curvature(m) = curvature(m) + &
                (target - (y(i, m) - y(j, m))**2 * &
                (1.0_dp + target / fitted) / fitted) / &
                (distance(i, j) * fitted)
            end do
          end do
          candidate(i, :) = y(i, :) + step * gradient / &
            max(abs(curvature), 1.0e-12_dp)
        end do
        call sammon_stress(distance, candidate, stress, total)
        if (stress <= previous) exit
        stress = previous
        step = step * 0.2_dp
        if (step <= 1.0e-3_dp) exit
      end do
      if (step <= 1.0e-3_dp .and. stress >= previous) exit
      step = min(0.5_dp, 1.5_dp * step)
      previous = stress
      do m = 1, k
        y(:, m) = candidate(:, m) - sum(candidate(:, m)) / real(n, dp)
      end do
      if (mod(iter, 10) == 0) then
        if (stress > checkpoint - tol) exit
        checkpoint = stress
      end if
    end do
    result%points = y
    result%stress = stress
    result%iterations = iter
    result%status = mass_success
    result%method = "Sammon mapping"
    call pairwise_fitted(y, result%fitted_distances)
  end subroutine sammon

  subroutine iso_mds(distance, result, initial, dimensions, maxit, &
      tolerance, minkowski_power)
    real(dp), intent(in) :: distance(:, :)
    type(mds_result), intent(out) :: result
    real(dp), intent(in), optional :: initial(:, :), tolerance, minkowski_power
    integer, intent(in), optional :: dimensions, maxit
    real(dp), allocatable :: y(:, :), gradient(:, :), candidate(:, :)
    real(dp), allocatable :: target(:), fitted(:), disparities(:)
    real(dp) :: stress, new_stress, tol, alpha, power
    integer :: n, k, mit, iter, status, i, j, pair, n_pairs

    n = size(distance, 1)
    k = 2
    if (present(dimensions)) k = dimensions
    mit = 100
    if (present(maxit)) mit = maxit
    tol = 1.0e-5_dp
    if (present(tolerance)) tol = tolerance
    power = 2.0_dp
    if (present(minkowski_power)) power = minkowski_power
    if (size(distance, 2) /= n .or. n <= k .or. k < 1 .or. power <= 0.0_dp) then
      result%status = mass_invalid_argument
      return
    end if
    if (present(initial)) then
      if (size(initial, 1) /= n .or. size(initial, 2) /= k) then
        result%status = mass_dimension_error
        return
      end if
      y = initial
    else
      call classical_mds(distance, k, y, status)
      if (status /= mass_success) then
        result%status = status
        return
      end if
    end if
    n_pairs = n * (n - 1) / 2
    allocate(target(n_pairs), fitted(n_pairs), disparities(n_pairs), &
      gradient(n, k), candidate(n, k))
    pair = 0
    do i = 1, n - 1
      do j = i + 1, n
        pair = pair + 1
        target(pair) = distance(i, j)
      end do
    end do
    call monotone_disparities(target, y, power, fitted, disparities)
    stress = sqrt(sum((fitted - disparities)**2) / &
      max(sum(fitted**2), tiny(1.0_dp)))
    alpha = 0.05_dp
    candidate = y
    new_stress = stress
    do iter = 1, mit
      call mds_gradient(y, disparities, power, gradient)
      do
        candidate = y - alpha * gradient
        call center_columns(candidate)
        call monotone_disparities(target, candidate, power, fitted, disparities)
        new_stress = sqrt(sum((fitted - disparities)**2) / &
          max(sum(fitted**2), tiny(1.0_dp)))
        if (new_stress <= stress .or. alpha < 1.0e-10_dp) exit
        alpha = 0.5_dp * alpha
      end do
      if (abs(stress - new_stress) < tol * max(1.0_dp, stress)) exit
      y = candidate
      stress = new_stress
      alpha = min(0.2_dp, 1.2_dp * alpha)
    end do
    result%points = candidate
    result%stress = new_stress
    result%fitted_distances = disparities
    result%iterations = iter
    result%status = merge(mass_success, mass_no_convergence, iter <= mit)
    result%method = "Kruskal nonmetric MDS"
  end subroutine iso_mds

  subroutine shepard(distance, points, target, fitted, disparities, &
      status, minkowski_power)
    real(dp), intent(in) :: distance(:, :), points(:, :)
    real(dp), allocatable, intent(out) :: target(:), fitted(:), disparities(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: minkowski_power
    real(dp) :: power
    integer :: n, i, j, pair
    n = size(distance, 1)
    power = 2.0_dp
    if (present(minkowski_power)) power = minkowski_power
    if (size(distance, 2) /= n .or. size(points, 1) /= n) then
      allocate(target(0), fitted(0), disparities(0))
      status = mass_dimension_error
      return
    end if
    allocate(target(n * (n - 1) / 2))
    pair = 0
    do i = 1, n - 1
      do j = i + 1, n
        pair = pair + 1
        target(pair) = distance(i, j)
      end do
    end do
    call monotone_disparities(target, points, power, fitted, disparities)
    status = mass_success
  end subroutine shepard

  subroutine sammon_stress(distance, y, stress, total)
    real(dp), intent(in) :: distance(:, :), y(:, :)
    real(dp), intent(out) :: stress, total
    real(dp) :: fitted, error
    integer :: i, j, n
    n = size(distance, 1)
    stress = 0.0_dp
    total = 0.0_dp
    do i = 1, n - 1
      do j = i + 1, n
        if (distance(i, j) <= tiny(1.0_dp)) cycle
        fitted = sqrt(sum((y(i, :) - y(j, :))**2))
        error = distance(i, j) - fitted
        stress = stress + error * error / distance(i, j)
        total = total + distance(i, j)
      end do
    end do
    stress = stress / max(total, tiny(1.0_dp))
  end subroutine sammon_stress

  subroutine pairwise_fitted(points, distances)
    real(dp), intent(in) :: points(:, :)
    real(dp), allocatable, intent(out) :: distances(:)
    integer :: n, i, j, pair
    n = size(points, 1)
    allocate(distances(n * (n - 1) / 2))
    pair = 0
    do i = 1, n - 1
      do j = i + 1, n
        pair = pair + 1
        distances(pair) = sqrt(sum((points(i, :) - points(j, :))**2))
      end do
    end do
  end subroutine pairwise_fitted

  subroutine monotone_disparities(target, points, power, fitted, disparities)
    real(dp), intent(in) :: target(:), points(:, :), power
    real(dp), allocatable, intent(out) :: fitted(:), disparities(:)
    real(dp), allocatable :: sorted_target(:), sorted_fit(:), iso(:)
    integer, allocatable :: order(:)
    integer :: n, i, j, pair, k, idx
    n = size(points, 1)
    allocate(fitted(size(target)), disparities(size(target)), &
      sorted_target(size(target)), sorted_fit(size(target)), &
      iso(size(target)), order(size(target)))
    pair = 0
    do i = 1, n - 1
      do j = i + 1, n
        pair = pair + 1
        fitted(pair) = sum(abs(points(i, :) - points(j, :))**power)**(1.0_dp / power)
      end do
    end do
    sorted_target = target
    order = [(i, i = 1, size(target))]
    call sort_with_index(sorted_target, order)
    sorted_fit = fitted(order)
    call isotonic_increasing(sorted_fit, fit=iso)
    disparities = 0.0_dp
    do k = 1, size(target)
      idx = order(k)
      disparities(idx) = iso(k)
    end do
  end subroutine monotone_disparities

  subroutine sort_with_index(values, order)
    real(dp), intent(inout) :: values(:)
    integer, intent(inout) :: order(:)
    real(dp) :: key
    integer :: key_order, i, j
    do i = 2, size(values)
      key = values(i)
      key_order = order(i)
      j = i - 1
      do while (j >= 1)
        if (values(j) <= key) exit
        values(j + 1) = values(j)
        order(j + 1) = order(j)
        j = j - 1
      end do
      values(j + 1) = key
      order(j + 1) = key_order
    end do
  end subroutine sort_with_index

  subroutine mds_gradient(points, disparities, power, gradient)
    real(dp), intent(in) :: points(:, :), disparities(:), power
    real(dp), intent(out) :: gradient(:, :)
    real(dp) :: fitted, factor
    integer :: n, i, j, pair, m
    n = size(points, 1)
    gradient = 0.0_dp
    pair = 0
    do i = 1, n - 1
      do j = i + 1, n
        pair = pair + 1
        fitted = sum(abs(points(i, :) - points(j, :))**power)**(1.0_dp / power)
        if (fitted <= tiny(1.0_dp)) cycle
        factor = 2.0_dp * (fitted - disparities(pair)) / fitted
        do m = 1, size(points, 2)
          gradient(i, m) = gradient(i, m) + factor * &
            sign(abs(points(i, m) - points(j, m))**(power - 1.0_dp), &
            points(i, m) - points(j, m))
          gradient(j, m) = gradient(j, m) - factor * &
            sign(abs(points(i, m) - points(j, m))**(power - 1.0_dp), &
            points(i, m) - points(j, m))
        end do
      end do
    end do
  end subroutine mds_gradient

  subroutine center_columns(x)
    real(dp), intent(inout) :: x(:, :)
    integer :: j
    do j = 1, size(x, 2)
      x(:, j) = x(:, j) - sum(x(:, j)) / real(size(x, 1), dp)
    end do
  end subroutine center_columns

end module mass_mds
