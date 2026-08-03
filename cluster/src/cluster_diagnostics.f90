! SPDX-License-Identifier: GPL-2.0-or-later
module cluster_diagnostics
  use, intrinsic :: iso_fortran_env, only: int64
  use fastcluster_kinds, only: dp
  use cluster_types, only: silhouette_result, gap_result, clustering_callback, &
    cluster_success, cluster_invalid_argument, cluster_not_converged
  use cluster_daisy, only: daisy
  implicit none
  private

  public :: silhouette
  public :: sort_silhouette
  public :: medoids
  public :: meanabsdev
  public :: size_diss
  public :: lower_to_upper_tri_inds
  public :: upper_to_lower_tri_inds
  public :: clus_gap
  public :: max_se
  public :: within_cluster_dispersion

contains

  subroutine silhouette(labels, distances, result)
    integer, intent(in) :: labels(:)
    real(dp), intent(in) :: distances(:, :)
    type(silhouette_result), intent(out) :: result

    integer, allocatable :: counts(:)
    real(dp) :: a, b, average
    integer :: c, i, j, k, n, own, best_cluster

    n = size(labels)
    if (n < 1 .or. size(distances, 1) /= n .or. size(distances, 2) /= n .or. &
        any(labels < 1)) then
      call fail_silhouette(result, cluster_invalid_argument, 'invalid labels or distance matrix')
      return
    end if
    k = maxval(labels)
    if (k < 1) then
      call fail_silhouette(result, cluster_invalid_argument, 'at least one cluster is required')
      return
    end if
    allocate(result%width(n), result%neighbor_distance(n), result%neighbor_cluster(n), counts(k))
    counts = 0
    do i = 1, n
      counts(labels(i)) = counts(labels(i)) + 1
    end do
    do i = 1, n
      own = labels(i)
      if (counts(own) <= 1) then
        result%width(i) = 0.0_dp
        result%neighbor_distance(i) = 0.0_dp
        result%neighbor_cluster(i) = 0
        cycle
      end if
      a = 0.0_dp
      do j = 1, n
        if (j /= i .and. labels(j) == own) a = a + distances(i, j)
      end do
      a = a / real(counts(own) - 1, dp)
      b = huge(1.0_dp)
      best_cluster = 0
      do c = 1, k
        if (c == own .or. counts(c) == 0) cycle
        average = 0.0_dp
        do j = 1, n
          if (labels(j) == c) average = average + distances(i, j)
        end do
        average = average / real(counts(c), dp)
        if (average < b) then
          b = average
          best_cluster = c
        end if
      end do
      if (best_cluster == 0) then
        result%width(i) = 0.0_dp
        result%neighbor_distance(i) = 0.0_dp
      else
        result%neighbor_distance(i) = b
        if (max(a, b) > 0.0_dp) then
          result%width(i) = (b - a) / max(a, b)
        else
          result%width(i) = 0.0_dp
        end if
      end if
      result%neighbor_cluster(i) = best_cluster
    end do
    result%average_width = sum(result%width) / real(n, dp)
    result%status = cluster_success
    result%message = 'ok'
  end subroutine silhouette

  subroutine sort_silhouette(labels, widths, order)
    integer, intent(in) :: labels(:)
    real(dp), intent(in) :: widths(:)
    integer, allocatable, intent(out) :: order(:)

    integer :: i, j, n, tmp

    n = size(labels)
    if (size(widths) /= n) then
      allocate(order(0))
      return
    end if
    allocate(order(n))
    order = [(i, i=1,n)]
    do i = 1, n - 1
      do j = i + 1, n
        if (labels(order(j)) < labels(order(i)) .or. &
            (labels(order(j)) == labels(order(i)) .and. widths(order(j)) > widths(order(i)))) then
          tmp = order(i)
          order(i) = order(j)
          order(j) = tmp
        end if
      end do
    end do
  end subroutine sort_silhouette

  subroutine medoids(labels, distances, indices, objective, status)
    integer, intent(in) :: labels(:)
    real(dp), intent(in) :: distances(:, :)
    integer, allocatable, intent(out) :: indices(:)
    real(dp), intent(out), optional :: objective
    integer, intent(out), optional :: status

    real(dp) :: best, value, total
    integer :: c, i, j, k, n

    n = size(labels)
    if (n < 1 .or. size(distances, 1) /= n .or. size(distances, 2) /= n .or. any(labels < 1)) then
      allocate(indices(0))
      if (present(objective)) objective = huge(1.0_dp)
      if (present(status)) status = cluster_invalid_argument
      return
    end if
    k = maxval(labels)
    allocate(indices(k))
    total = 0.0_dp
    do c = 1, k
      best = huge(1.0_dp)
      indices(c) = 0
      do i = 1, n
        if (labels(i) /= c) cycle
        value = 0.0_dp
        do j = 1, n
          if (labels(j) == c) value = value + distances(i, j)
        end do
        if (value < best) then
          best = value
          indices(c) = i
        end if
      end do
      if (indices(c) == 0) then
        if (present(status)) status = cluster_invalid_argument
        return
      end if
      total = total + best
    end do
    if (present(objective)) objective = total
    if (present(status)) status = cluster_success
  end subroutine medoids

  real(dp) function meanabsdev(x, center) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: center

    real(dp) :: location

    if (size(x) == 0) then
      value = 0.0_dp
      return
    end if
    if (present(center)) then
      location = center
    else
      location = sum(x) / real(size(x), dp)
    end if
    value = sum(abs(x - location)) / real(size(x), dp)
  end function meanabsdev

  integer function size_diss(length_value) result(n)
    integer, intent(in) :: length_value
    real(dp) :: root
    if (length_value < 0) then
      n = -1
      return
    end if
    root = (1.0_dp + sqrt(1.0_dp + 8.0_dp * real(length_value, dp))) / 2.0_dp
    n = nint(root)
    if (n * (n - 1) / 2 /= length_value) n = -1
  end function size_diss

  subroutine lower_to_upper_tri_inds(n, indices)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: indices(:)

    integer :: i, j, k, lower_index

    if (n < 2) then
      allocate(indices(0))
      return
    end if
    allocate(indices(n*(n-1)/2))
    k = 0
    do j = 2, n
      do i = 1, j - 1
        k = k + 1
        lower_index = j + (i - 1)*n - i*(i + 1)/2
        indices(k) = lower_index
      end do
    end do
  end subroutine lower_to_upper_tri_inds

  subroutine upper_to_lower_tri_inds(n, indices)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: indices(:)

    integer, allocatable :: inverse(:)
    integer :: i

    call lower_to_upper_tri_inds(n, inverse)
    allocate(indices(size(inverse)))
    do i = 1, size(inverse)
      indices(inverse(i)) = i
    end do
  end subroutine upper_to_lower_tri_inds

  subroutine clus_gap(x, k_max, b_references, cluster_fun, result, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: k_max, b_references
    procedure(clustering_callback) :: cluster_fun
    type(gap_result), intent(out) :: result
    integer, intent(in), optional :: seed

    real(dp), allocatable :: reference(:, :), distances(:, :), log_refs(:, :)
    integer, allocatable :: labels(:)
    real(dp), allocatable :: xmin(:), xmax(:)
    integer(int64) :: state
    real(dp) :: dispersion, u
    integer :: b, i, j, k, n, p, status
    character(len=:), allocatable :: message

    n = size(x, 1)
    p = size(x, 2)
    if (n < 2 .or. p < 1 .or. k_max < 1 .or. k_max > n .or. b_references < 1) then
      call fail_gap(result, cluster_invalid_argument, 'invalid gap-statistic arguments')
      return
    end if
    allocate(result%k(k_max), result%log_w(k_max), result%gap(k_max), result%se(k_max))
    allocate(log_refs(b_references, k_max), reference(n, p), xmin(p), xmax(p))
    result%k = [(k, k=1,k_max)]
    do j = 1, p
      xmin(j) = minval(x(:, j))
      xmax(j) = maxval(x(:, j))
    end do
    call daisy(x, distances, status=status, message=message)
    if (status /= cluster_success) then
      call fail_gap(result, status, message)
      return
    end if
    do k = 1, k_max
      call cluster_fun(x, k, labels, status)
      if (status /= cluster_success .or. size(labels) /= n) then
        call fail_gap(result, cluster_not_converged, 'clustering callback failed on observed data')
        return
      end if
      dispersion = max(tiny(1.0_dp), within_cluster_dispersion(labels, distances))
      result%log_w(k) = log(dispersion)
    end do
    state = 65537_int64
    if (present(seed)) state = max(1_int64, int(seed, int64))
    do b = 1, b_references
      do i = 1, n
        do j = 1, p
          call lcg_uniform(state, u)
          reference(i, j) = xmin(j) + u * (xmax(j) - xmin(j))
        end do
      end do
      call daisy(reference, distances, status=status, message=message)
      if (status /= cluster_success) then
        call fail_gap(result, status, message)
        return
      end if
      do k = 1, k_max
        call cluster_fun(reference, k, labels, status)
        if (status /= cluster_success .or. size(labels) /= n) then
          call fail_gap(result, cluster_not_converged, 'clustering callback failed on reference data')
          return
        end if
        dispersion = max(tiny(1.0_dp), within_cluster_dispersion(labels, distances))
        log_refs(b, k) = log(dispersion)
      end do
    end do
    do k = 1, k_max
      result%gap(k) = sum(log_refs(:, k)) / real(b_references, dp) - result%log_w(k)
      if (b_references > 1) then
        result%se(k) = sqrt(sum((log_refs(:, k) - &
          sum(log_refs(:, k))/real(b_references, dp))**2) / real(b_references-1, dp)) * &
          sqrt(1.0_dp + 1.0_dp/real(b_references, dp))
      else
        result%se(k) = 0.0_dp
      end if
    end do
    result%selected_k = max_se(result%gap, result%se, 'Tibs2001SEmax', 1.0_dp)
    result%status = cluster_success
    result%message = 'ok'
  end subroutine clus_gap

  integer function max_se(f, se, method, se_factor) result(index_value)
    real(dp), intent(in) :: f(:), se(:)
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: se_factor

    character(len=:), allocatable :: name
    real(dp) :: factor, threshold
    integer :: i, imax

    if (size(f) < 1 .or. size(se) /= size(f)) then
      index_value = 0
      return
    end if
    name = 'tibs2001semax'
    if (present(method)) name = lower_ascii(trim(adjustl(method)))
    factor = 1.0_dp
    if (present(se_factor)) factor = se_factor
    imax = maxloc(f, dim=1)
    select case (name)
    case ('firstmax')
      index_value = imax
    case ('globalmax')
      index_value = imax
    case ('globalsemax')
      threshold = f(imax) - factor * se(imax)
      index_value = 1
      do i = 1, imax
        if (f(i) >= threshold) then
          index_value = i
          exit
        end if
      end do
    case ('firstsemax')
      index_value = imax
      do i = 1, size(f) - 1
        if (f(i) >= f(i+1) - factor*se(i+1)) then
          index_value = i
          exit
        end if
      end do
    case default
      index_value = size(f)
      do i = 1, size(f) - 1
        if (f(i) >= f(i+1) - factor*se(i+1)) then
          index_value = i
          exit
        end if
      end do
    end select
  end function max_se

  real(dp) function within_cluster_dispersion(labels, distances) result(value)
    integer, intent(in) :: labels(:)
    real(dp), intent(in) :: distances(:, :)

    integer, allocatable :: counts(:)
    integer :: c, i, j, k

    k = maxval(labels)
    allocate(counts(k))
    counts = 0
    do i = 1, size(labels)
      counts(labels(i)) = counts(labels(i)) + 1
    end do
    value = 0.0_dp
    do c = 1, k
      if (counts(c) <= 1) cycle
      do i = 1, size(labels) - 1
        if (labels(i) /= c) cycle
        do j = i + 1, size(labels)
          if (labels(j) == c) value = value + distances(i, j) / real(counts(c), dp)
        end do
      end do
    end do
  end function within_cluster_dispersion

  subroutine lcg_uniform(state, u)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: u
    integer(int64), parameter :: modulus = 2147483647_int64
    integer(int64), parameter :: multiplier = 48271_int64
    state = modulo(multiplier * state, modulus)
    if (state == 0_int64) state = 1_int64
    u = real(state, dp) / real(modulus, dp)
  end subroutine lcg_uniform

  subroutine fail_silhouette(result, status, message)
    type(silhouette_result), intent(out) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    result%status = status
    result%message = message
  end subroutine fail_silhouette

  subroutine fail_gap(result, status, message)
    type(gap_result), intent(out) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    result%status = status
    result%message = message
  end subroutine fail_gap

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code
    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    end do
  end function lower_ascii

end module cluster_diagnostics
