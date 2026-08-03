! SPDX-License-Identifier: GPL-2.0-or-later
module cluster_partition
  use, intrinsic :: iso_fortran_env, only: int64
  use fastcluster_kinds, only: dp
  use cluster_types, only: partition_result, cluster_success, cluster_invalid_argument, &
    cluster_allocation_failure, cluster_not_converged
  use cluster_daisy, only: daisy
  implicit none
  private

  public :: pam
  public :: pam_distance
  public :: clara
  public :: fanny
  public :: assign_to_medoids

contains

  subroutine pam(x, k, result, metric, max_iter)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: k
    type(partition_result), intent(out) :: result
    character(len=*), intent(in), optional :: metric
    integer, intent(in), optional :: max_iter

    real(dp), allocatable :: distances(:, :)
    integer :: status
    character(len=:), allocatable :: message

    if (present(metric)) then
      call daisy(x, distances, metric=metric, status=status, message=message)
    else
      call daisy(x, distances, status=status, message=message)
    end if
    if (status /= cluster_success) then
      call fail_partition(result, status, message)
      return
    end if
    if (present(max_iter)) then
      call pam_distance(distances, k, result, max_iter=max_iter)
    else
      call pam_distance(distances, k, result)
    end if
  end subroutine pam

  subroutine pam_distance(distances, k, result, max_iter, initial_medoids)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(in) :: k
    type(partition_result), intent(out) :: result
    integer, intent(in), optional :: max_iter
    integer, intent(in), optional :: initial_medoids(:)

    integer, allocatable :: medoids(:), trial(:), labels(:)
    real(dp), allocatable :: nearest(:), trial_nearest(:)
    real(dp) :: best_objective, objective, tol
    integer :: best_h, best_m, h, iter, m, n, niter
    logical :: improved

    call initialize_partition(result)
    n = size(distances, 1)
    if (n < 2 .or. size(distances, 2) /= n .or. k < 1 .or. k > n) then
      call fail_partition(result, cluster_invalid_argument, &
        'distances must be square and 1 <= k <= n')
      return
    end if
    if (any(distances < 0.0_dp)) then
      call fail_partition(result, cluster_invalid_argument, 'distances must be nonnegative')
      return
    end if
    niter = 100
    if (present(max_iter)) niter = max(0, max_iter)
    allocate(medoids(k), trial(k), labels(n), nearest(n), trial_nearest(n))
    if (present(initial_medoids)) then
      if (size(initial_medoids) /= k .or. any(initial_medoids < 1) .or. &
          any(initial_medoids > n) .or. has_duplicates(initial_medoids)) then
        call fail_partition(result, cluster_invalid_argument, 'invalid initial medoids')
        return
      end if
      medoids = initial_medoids
    else
      call build_medoids(distances, k, medoids)
    end if
    call assign_to_medoids(distances, medoids, labels, nearest, best_objective)
    tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, best_objective)

    do iter = 1, niter
      improved = .false.
      best_h = 0
      best_m = 0
      objective = best_objective
      do m = 1, k
        do h = 1, n
          if (any(medoids == h)) cycle
          trial = medoids
          trial(m) = h
          call assign_to_medoids(distances, trial, labels, trial_nearest, objective)
          if (objective < best_objective - tol) then
            best_objective = objective
            best_h = h
            best_m = m
            improved = .true.
          end if
        end do
      end do
      if (.not. improved) exit
      medoids(best_m) = best_h
      call assign_to_medoids(distances, medoids, labels, nearest, best_objective)
    end do

    call assign_to_medoids(distances, medoids, labels, nearest, best_objective)
    result%n = n
    result%k = k
    result%medoids = medoids
    result%clustering = labels
    result%objective = best_objective / real(n, dp)
    result%iterations = min(iter, niter)
    if (iter > niter .and. niter > 0) then
      result%status = cluster_not_converged
      result%message = 'maximum PAM swap iterations reached'
    else
      result%status = cluster_success
      result%message = 'ok'
    end if
  end subroutine pam_distance

  subroutine clara(x, k, result, samples, sample_size, metric, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: k
    type(partition_result), intent(out) :: result
    integer, intent(in), optional :: samples, sample_size, seed
    character(len=*), intent(in), optional :: metric

    real(dp), allocatable :: all_dist(:, :), sub_dist(:, :), nearest(:)
    integer, allocatable :: subset(:), medoids(:), labels(:), sub_medoids(:)
    type(partition_result) :: sub_result
    integer(int64) :: state
    integer :: b, best_b, i, j, n, nsamp, ssize, st
    real(dp) :: objective, best_objective
    character(len=:), allocatable :: message

    call initialize_partition(result)
    n = size(x, 1)
    nsamp = 5
    if (present(samples)) nsamp = samples
    ssize = min(n, max(40 + 2*k, 2*k))
    if (present(sample_size)) ssize = sample_size
    if (k < 1 .or. k > n .or. nsamp < 1 .or. ssize < k .or. ssize > n) then
      call fail_partition(result, cluster_invalid_argument, 'invalid CLARA dimensions')
      return
    end if
    if (present(metric)) then
      call daisy(x, all_dist, metric=metric, status=st, message=message)
    else
      call daisy(x, all_dist, status=st, message=message)
    end if
    if (st /= cluster_success) then
      call fail_partition(result, st, message)
      return
    end if
    allocate(subset(ssize), sub_dist(ssize, ssize), medoids(k), labels(n), nearest(n))
    state = 104729_int64
    if (present(seed)) state = max(1_int64, int(seed, int64))
    best_objective = huge(1.0_dp)
    best_b = 0
    do b = 1, nsamp
      call sample_without_replacement(n, ssize, state, subset)
      do i = 1, ssize
        do j = 1, ssize
          sub_dist(i, j) = all_dist(subset(i), subset(j))
        end do
      end do
      call pam_distance(sub_dist, k, sub_result)
      if (.not. sub_result%ok()) cycle
      allocate(sub_medoids(k))
      do i = 1, k
        sub_medoids(i) = subset(sub_result%medoids(i))
      end do
      call assign_to_medoids(all_dist, sub_medoids, labels, nearest, objective)
      if (objective < best_objective) then
        best_objective = objective
        medoids = sub_medoids
        best_b = b
      end if
      deallocate(sub_medoids)
    end do
    if (best_b == 0) then
      call fail_partition(result, cluster_not_converged, 'no CLARA sample produced a solution')
      return
    end if
    call assign_to_medoids(all_dist, medoids, labels, nearest, objective)
    result%n = n
    result%k = k
    result%medoids = medoids
    result%clustering = labels
    result%objective = objective / real(n, dp)
    result%iterations = nsamp
    result%status = cluster_success
    result%message = 'ok'
  end subroutine clara

  subroutine fanny(x, k, result, membership_exponent, metric, tolerance, max_iter)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: k
    type(partition_result), intent(out) :: result
    real(dp), intent(in), optional :: membership_exponent, tolerance
    character(len=*), intent(in), optional :: metric
    integer, intent(in), optional :: max_iter

    real(dp), allocatable :: distances(:, :), u(:, :), unew(:, :), nearest(:)
    integer, allocatable :: medoids(:), labels(:)
    type(partition_result) :: initial
    real(dp) :: exponent, tol, power, objective, best, value
    integer :: c, i, iter, j, n, niter, st
    character(len=:), allocatable :: message

    call initialize_partition(result)
    n = size(x, 1)
    exponent = 2.0_dp
    if (present(membership_exponent)) exponent = membership_exponent
    tol = 1.0e-6_dp
    if (present(tolerance)) tol = tolerance
    niter = 200
    if (present(max_iter)) niter = max_iter
    if (k < 2 .or. k > n .or. exponent <= 1.0_dp .or. tol <= 0.0_dp .or. niter < 1) then
      call fail_partition(result, cluster_invalid_argument, 'invalid FANNY arguments')
      return
    end if
    if (present(metric)) then
      call daisy(x, distances, metric=metric, status=st, message=message)
    else
      call daisy(x, distances, status=st, message=message)
    end if
    if (st /= cluster_success) then
      call fail_partition(result, st, message)
      return
    end if
    call pam_distance(distances, k, initial)
    if (.not. initial%ok()) then
      result = initial
      return
    end if
    allocate(medoids(k), labels(n), nearest(n), u(n, k), unew(n, k))
    medoids = initial%medoids
    power = 2.0_dp / (exponent - 1.0_dp)
    call update_membership(distances, medoids, power, u)

    do iter = 1, niter
      do c = 1, k
        best = huge(1.0_dp)
        medoids(c) = 1
        do j = 1, n
          value = sum((u(:, c)**exponent) * distances(:, j))
          if (value < best) then
            best = value
            medoids(c) = j
          end if
        end do
      end do
      call make_unique_medoids(distances, medoids)
      call update_membership(distances, medoids, power, unew)
      if (maxval(abs(unew - u)) <= tol) exit
      u = unew
    end do
    u = unew
    do i = 1, n
      labels(i) = maxloc(u(i, :), dim=1)
    end do
    objective = 0.0_dp
    do i = 1, n
      do c = 1, k
        objective = objective + u(i, c)**exponent * distances(i, medoids(c))
      end do
    end do
    result%n = n
    result%k = k
    result%medoids = medoids
    result%clustering = labels
    result%membership = u
    result%objective = objective / real(n, dp)
    result%iterations = min(iter, niter)
    if (iter > niter) then
      result%status = cluster_not_converged
      result%message = 'maximum FANNY iterations reached'
    else
      result%status = cluster_success
      result%message = 'ok'
    end if
  end subroutine fanny

  subroutine build_medoids(distances, k, medoids)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(in) :: k
    integer, intent(out) :: medoids(k)

    real(dp), allocatable :: nearest(:)
    real(dp) :: best_gain, gain
    integer :: c, candidate, i, n

    n = size(distances, 1)
    medoids(1) = minloc(sum(distances, dim=1), dim=1)
    allocate(nearest(n))
    nearest = distances(:, medoids(1))
    do c = 2, k
      best_gain = -huge(1.0_dp)
      medoids(c) = 1
      do candidate = 1, n
        if (any(medoids(1:c-1) == candidate)) cycle
        gain = 0.0_dp
        do i = 1, n
          gain = gain + max(0.0_dp, nearest(i) - distances(i, candidate))
        end do
        if (gain > best_gain) then
          best_gain = gain
          medoids(c) = candidate
        end if
      end do
      nearest = min(nearest, distances(:, medoids(c)))
    end do
  end subroutine build_medoids

  subroutine assign_to_medoids(distances, medoids, labels, nearest, objective)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(in) :: medoids(:)
    integer, intent(out) :: labels(:)
    real(dp), intent(out) :: nearest(:)
    real(dp), intent(out) :: objective

    integer :: c, i

    do i = 1, size(distances, 1)
      labels(i) = 1
      nearest(i) = distances(i, medoids(1))
      do c = 2, size(medoids)
        if (distances(i, medoids(c)) < nearest(i)) then
          nearest(i) = distances(i, medoids(c))
          labels(i) = c
        end if
      end do
    end do
    objective = sum(nearest)
  end subroutine assign_to_medoids

  subroutine update_membership(distances, medoids, power, u)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(in) :: medoids(:)
    real(dp), intent(in) :: power
    real(dp), intent(out) :: u(:, :)

    real(dp) :: denom
    integer :: c, i, k, zcount

    k = size(medoids)
    do i = 1, size(distances, 1)
      zcount = 0
      do c = 1, k
        if (distances(i, medoids(c)) <= sqrt(epsilon(1.0_dp))) zcount = zcount + 1
      end do
      if (zcount > 0) then
        do c = 1, k
          if (distances(i, medoids(c)) <= sqrt(epsilon(1.0_dp))) then
            u(i, c) = 1.0_dp / real(zcount, dp)
          else
            u(i, c) = 0.0_dp
          end if
        end do
      else
        denom = 0.0_dp
        do c = 1, k
          u(i, c) = distances(i, medoids(c))**(-power)
          denom = denom + u(i, c)
        end do
        u(i, :) = u(i, :) / denom
      end if
    end do
  end subroutine update_membership

  subroutine make_unique_medoids(distances, medoids)
    real(dp), intent(in) :: distances(:, :)
    integer, intent(inout) :: medoids(:)

    integer :: c, candidate
    real(dp) :: best

    do c = 2, size(medoids)
      if (any(medoids(1:c-1) == medoids(c))) then
        best = -1.0_dp
        do candidate = 1, size(distances, 1)
          if (any(medoids(1:c-1) == candidate)) cycle
          if (minval(distances(candidate, medoids(1:c-1))) > best) then
            best = minval(distances(candidate, medoids(1:c-1)))
            medoids(c) = candidate
          end if
        end do
      end if
    end do
  end subroutine make_unique_medoids

  subroutine sample_without_replacement(n, m, state, sample)
    integer, intent(in) :: n, m
    integer(int64), intent(inout) :: state
    integer, intent(out) :: sample(m)

    integer, allocatable :: pool(:)
    integer :: i, j, tmp
    real(dp) :: u

    allocate(pool(n))
    pool = [(i, i=1,n)]
    do i = 1, m
      call lcg_uniform(state, u)
      j = i + int(u * real(n - i + 1, dp))
      j = min(n, max(i, j))
      tmp = pool(i)
      pool(i) = pool(j)
      pool(j) = tmp
      sample(i) = pool(i)
    end do
  end subroutine sample_without_replacement

  subroutine lcg_uniform(state, u)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: u

    integer(int64), parameter :: modulus = 2147483647_int64
    integer(int64), parameter :: multiplier = 48271_int64

    state = modulo(multiplier * state, modulus)
    if (state == 0_int64) state = 1_int64
    u = real(state, dp) / real(modulus, dp)
  end subroutine lcg_uniform

  logical function has_duplicates(values) result(yes)
    integer, intent(in) :: values(:)
    integer :: i, j
    yes = .false.
    do i = 1, size(values) - 1
      do j = i + 1, size(values)
        if (values(i) == values(j)) then
          yes = .true.
          return
        end if
      end do
    end do
  end function has_duplicates

  subroutine initialize_partition(result)
    type(partition_result), intent(out) :: result
    result%status = cluster_success
    result%message = 'ok'
  end subroutine initialize_partition

  subroutine fail_partition(result, status, message)
    type(partition_result), intent(out) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    result%status = status
    result%message = message
  end subroutine fail_partition

end module cluster_partition
