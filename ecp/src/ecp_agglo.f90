module ecp_agglo
  use ecp_kinds, only: dp
  use ecp_types, only: agglo_result
  use ecp_energy, only: get_within, get_between
  implicit none
  private

  public :: e_agglo

contains

  pure function e_agglo(x, member, alpha, penalty) result(res)
    real(dp), intent(in) :: x(:,:) !! Time series matrix with observations in rows.
    integer, intent(in), optional :: member(:) !! Initial contiguous cluster labels; default gives each row its own cluster.
    real(dp), intent(in), optional :: alpha !! Energy-distance exponent in (0, 2]; default is one.
    real(dp), intent(in), optional :: penalty(:) !! Optional additive penalty for each clustering stage from N to one clusters.
    type(agglo_result) :: res
    integer, allocatable :: mem(:), labels(:), sizes(:), left(:), right(:), lm(:)
    integer, allocatable :: progression(:,:), merged(:,:)
    logical, allocatable :: open(:)
    real(dp), allocatable :: dmat(:,:), within(:), fit(:)
    real(dp) :: aa, best, candidate
    integer :: nobs, ncl, i, j, stage, new_index, ll, rr, n1, n2, n3, total
    integer :: best_i, best_j, best_row, count_est

    nobs = size(x, 1)
    aa = 1.0_dp
    if (present(alpha)) aa = alpha
    if (present(member)) then
      mem = member
    else
      allocate(mem(nobs))
      mem = [(i, i=1, nobs)]
    end if
    if (size(mem) /= nobs .or. nobs == 0) then
      allocate(res%estimates(0), res%cluster(0), res%merged(0, 0), res%progression(0, 0), res%fit(0))
      return
    end if

    call relabel_contiguous(mem, labels, ncl)
    if (ncl <= 0) then
      allocate(res%estimates(0), res%cluster(0), res%merged(0, 0), res%progression(0, 0), res%fit(0))
      return
    end if
    if (.not. contiguous_membership(labels)) then
      allocate(res%estimates(0), res%cluster(0), res%merged(0, 0), res%progression(0, 0), res%fit(0))
      return
    end if

    allocate(sizes(2*ncl), source=0)
    allocate(left(2*ncl - 1), source=0)
    allocate(right(2*ncl - 1), source=0)
    allocate(lm(2*ncl - 1), source=0)
    allocate(open(2*ncl - 1), source=.false.)
    allocate(dmat(2*ncl, 2*ncl), source=huge(1.0_dp))
    allocate(within(ncl), source=0.0_dp)
    allocate(fit(ncl), source=0.0_dp)
    allocate(merged(max(0, ncl - 1), 2), source=0)
    allocate(progression(ncl, ncl + 1), source=0)

    do i = 1, ncl
      sizes(i) = count(labels == i)
      open(i) = .true.
      lm(i) = i
      within(i) = cluster_within(x, labels, i, aa)
    end do
    if (ncl == 1) then
      left(1) = 1
      right(1) = 1
    else
      do i = 2, ncl - 1
        left(i) = i - 1
        right(i) = i + 1
      end do
      left(ncl) = ncl - 1
      right(ncl) = 1
      left(1) = ncl
      right(1) = 2
    end if

    dmat = huge(1.0_dp)
    do i = 1, 2*ncl
      dmat(i, i) = 0.0_dp
    end do
    do i = 1, ncl - 1
      do j = i + 1, ncl
        dmat(i, j) = cluster_between(x, labels, i, j, aa) - within(i) - within(j)
        dmat(j, i) = dmat(i, j)
      end do
    end do

    progression(1, 1) = 1
    do i = 2, ncl + 1
      progression(1, i) = progression(1, i - 1) + sizes(i - 1)
    end do
    do i = 1, ncl
      fit(1) = fit(1) + dmat(i, left(i)) + dmat(i, right(i))
    end do

    do stage = 1, ncl - 1
      new_index = ncl + stage
      best = -huge(1.0_dp)
      best_i = 0
      best_j = 0
      do i = 1, new_index - 1
        if (.not. open(i)) cycle
        j = right(i)
        ll = left(i)
        rr = right(j)
        candidate = fit(stage) - 2.0_dp*(dmat(i, j) + dmat(i, ll) + dmat(j, rr))
        n1 = sizes(i)
        n2 = sizes(j)
        n3 = sizes(ll)
        total = n1 + n2 + n3
        candidate = candidate + 2.0_dp*((n1 + n3)*dmat(i, ll) + &
          (n2 + n3)*dmat(j, ll) - n3*dmat(i, j))/real(total, dp)
        n3 = sizes(rr)
        total = n1 + n2 + n3
        candidate = candidate + 2.0_dp*((n1 + n3)*dmat(i, rr) + &
          (n2 + n3)*dmat(j, rr) - n3*dmat(i, j))/real(total, dp)
        if (candidate > best) then
          best = candidate
          best_i = i
          best_j = j
        end if
      end do
      if (best_i == 0) exit
      fit(stage + 1) = best
      i = best_i
      j = best_j
      if (i <= ncl) then
        merged(stage, 1) = -i
      else
        merged(stage, 1) = i - ncl
      end if
      if (j <= ncl) then
        merged(stage, 2) = -j
      else
        merged(stage, 2) = j - ncl
      end if
      ll = left(i)
      rr = right(j)
      left(new_index) = ll
      right(new_index) = rr
      right(ll) = new_index
      left(rr) = new_index
      open(i) = .false.
      open(j) = .false.
      open(new_index) = .true.
      n1 = sizes(i)
      n2 = sizes(j)
      sizes(new_index) = n1 + n2
      progression(stage + 1, :) = progression(stage, :)
      progression(stage + 1, lm(j)) = 0
      lm(new_index) = lm(i)
      do j = 1, new_index - 1
        if (.not. open(j)) cycle
        n3 = sizes(j)
        total = n1 + n2 + n3
        dmat(new_index, j) = ((total - n2)*dmat(i, j) + (total - n1)*dmat(best_j, j) - &
          n3*dmat(i, best_j))/real(total, dp)
        dmat(j, new_index) = dmat(new_index, j)
      end do
    end do

    if (present(penalty)) then
      do i = 1, min(size(penalty), ncl)
        fit(i) = fit(i) + penalty(i)
      end do
    end if
    best_row = maxloc(fit, dim=1)
    count_est = count(progression(best_row, :) > 0)
    allocate(res%estimates(count_est))
    res%estimates = pack(progression(best_row, :), progression(best_row, :) > 0)
    call sort_int(res%estimates)
    if (size(res%estimates) > 0) then
      if (res%estimates(1) /= 1 .and. res%estimates(size(res%estimates)) == nobs + 1) then
        res%estimates = res%estimates(:size(res%estimates) - 1)
      end if
    end if
    res%cluster = cluster_from_boundaries(res%estimates, nobs)
    res%merged = merged
    res%progression = progression
    res%fit = fit
  end function e_agglo

  pure subroutine relabel_contiguous(member, labels, ncl)
    integer, intent(in) :: member(:) !! Original integer cluster labels.
    integer, allocatable, intent(out) :: labels(:) !! Labels recoded by first occurrence to 1, 2, ... .
    integer, intent(out) :: ncl !! Number of distinct labels encountered.
    integer, allocatable :: seen(:)
    integer :: i, j, idx

    allocate(labels(size(member)), source=0)
    allocate(seen(size(member)), source=0)
    ncl = 0
    do i = 1, size(member)
      idx = 0
      do j = 1, ncl
        if (seen(j) == member(i)) then
          idx = j
          exit
        end if
      end do
      if (idx == 0) then
        ncl = ncl + 1
        seen(ncl) = member(i)
        idx = ncl
      end if
      labels(i) = idx
    end do
  end subroutine relabel_contiguous

  pure logical function contiguous_membership(labels) result(ok)
    integer, intent(in) :: labels(:) !! Recoded cluster labels to verify for contiguous segments.
    integer :: i

    ok = .true.
    do i = 2, size(labels)
      if (labels(i) < labels(i - 1)) then
        ok = .false.
        return
      end if
      if (labels(i) > labels(i - 1) + 1) then
        ok = .false.
        return
      end if
    end do
  end function contiguous_membership

  pure real(dp) function cluster_within(x, labels, label, alpha) result(value)
    real(dp), intent(in) :: x(:,:) !! Full observation matrix.
    integer, intent(in) :: labels(:) !! Cluster label for each observation row.
    integer, intent(in) :: label !! Cluster whose within distance is requested.
    real(dp), intent(in) :: alpha !! Energy-distance exponent.
    real(dp), allocatable :: tmp(:,:)
    integer :: i, r, n

    n = count(labels == label)
    allocate(tmp(n, size(x, 2)))
    r = 0
    do i = 1, size(labels)
      if (labels(i) == label) then
        r = r + 1
        tmp(r, :) = x(i, :)
      end if
    end do
    value = get_within(alpha, tmp)
  end function cluster_within

  pure real(dp) function cluster_between(x, labels, left_label, right_label, alpha) result(value)
    real(dp), intent(in) :: x(:,:) !! Full observation matrix.
    integer, intent(in) :: labels(:) !! Cluster label for each observation row.
    integer, intent(in) :: left_label !! First cluster label.
    integer, intent(in) :: right_label !! Second cluster label.
    real(dp), intent(in) :: alpha !! Energy-distance exponent.
    real(dp), allocatable :: a(:,:), b(:,:)
    integer :: i, r, n1, n2

    n1 = count(labels == left_label)
    n2 = count(labels == right_label)
    allocate(a(n1, size(x, 2)), b(n2, size(x, 2)))
    r = 0
    do i = 1, size(labels)
      if (labels(i) == left_label) then
        r = r + 1
        a(r, :) = x(i, :)
      end if
    end do
    r = 0
    do i = 1, size(labels)
      if (labels(i) == right_label) then
        r = r + 1
        b(r, :) = x(i, :)
      end if
    end do
    value = get_between(alpha, a, b)
  end function cluster_between

  pure function cluster_from_boundaries(boundaries, nobs) result(cluster)
    integer, intent(in) :: boundaries(:) !! Segment start indices, optionally including nobs+1 as an end boundary.
    integer, intent(in) :: nobs !! Total number of observations to label.
    integer, allocatable :: cluster(:)
    integer, allocatable :: b(:)
    integer :: i, j, nb, start_i, end_i, label

    allocate(cluster(nobs), source=1)
    if (size(boundaries) == 0) return
    if (boundaries(1) == 1) then
      b = boundaries
    else
      allocate(b(size(boundaries) + 1))
      b(1) = 1
      b(2:) = boundaries
    end if
    nb = size(b)
    label = 0
    do i = 1, nb
      start_i = b(i)
      if (i < nb) then
        end_i = b(i + 1) - 1
      else if (b(nb) == nobs + 1) then
        cycle
      else
        end_i = nobs
      end if
      if (start_i > nobs) cycle
      label = label + 1
      do j = start_i, min(end_i, nobs)
        cluster(j) = label
      end do
    end do
  end function cluster_from_boundaries

  pure subroutine sort_int(x)
    integer, intent(inout) :: x(:) !! Integer vector sorted in ascending order in place.
    integer :: i, j, key

    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_int

end module ecp_agglo
