module ecp_energy
  use ecp_kinds, only: dp
  implicit none
  private

  public :: point_distance_alpha, distance_matrix_alpha
  public :: get_within, get_between, energy_split_stat
  public :: split_point, ks_split_stat, windowed_energy_split_stat

contains

  pure real(dp) function point_distance_alpha(x, y, alpha) result(d)
    real(dp), intent(in) :: x(:) !! First observation vector.
    real(dp), intent(in) :: y(:) !! Second observation vector with the same dimension as x.
    real(dp), intent(in) :: alpha !! Distance exponent in the interval (0, 2].

    d = sum((x - y)**2)**(0.5_dp*alpha)
  end function point_distance_alpha

  pure function distance_matrix_alpha(x, alpha) result(d)
    real(dp), intent(in) :: x(:,:) !! Observations, with rows indexing time and columns indexing variables.
    real(dp), intent(in) :: alpha !! Distance exponent in the interval (0, 2].
    real(dp), allocatable :: d(:,:)
    integer :: i, j, n

    n = size(x, 1)
    allocate(d(n, n), source=0.0_dp)
    do j = 2, n
      do i = 1, j - 1
        d(i, j) = point_distance_alpha(x(i, :), x(j, :), alpha)
        d(j, i) = d(i, j)
      end do
    end do
  end function distance_matrix_alpha

  pure real(dp) function get_within(alpha, x) result(value)
    real(dp), intent(in) :: alpha !! Distance exponent in the interval (0, 2].
    real(dp), intent(in) :: x(:,:) !! Sample matrix whose rows are observations.
    integer :: i, j, n

    n = size(x, 1)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    value = 0.0_dp
    do j = 1, n
      do i = 1, n
        value = value + point_distance_alpha(x(i, :), x(j, :), alpha)
      end do
    end do
    value = value/real(n*n, dp)
  end function get_within

  pure real(dp) function get_between(alpha, x, y) result(value)
    real(dp), intent(in) :: alpha !! Distance exponent in the interval (0, 2].
    real(dp), intent(in) :: x(:,:) !! First sample matrix whose rows are observations.
    real(dp), intent(in) :: y(:,:) !! Second sample matrix with the same number of columns as x.
    integer :: i, j, n, m

    n = size(x, 1)
    m = size(y, 1)
    if (n == 0 .or. m == 0) then
      value = 0.0_dp
      return
    end if
    value = 0.0_dp
    do j = 1, m
      do i = 1, n
        value = value + point_distance_alpha(x(i, :), y(j, :), alpha)
      end do
    end do
    value = 2.0_dp*value/real(n*m, dp)
  end function get_between

  pure real(dp) function energy_split_stat(d, first, split, last, cp3o_scale) result(value)
    real(dp), intent(in) :: d(:,:) !! Pairwise alpha-powered distance matrix.
    integer, intent(in) :: first !! First observation index in the combined interval.
    integer, intent(in) :: split !! Last observation index of the left segment.
    integer, intent(in) :: last !! Last observation index of the right segment.
    logical, intent(in) :: cp3o_scale !! Use CP3O scaling n*m/(n+m)^2 instead of divisive n*m/(n+m).
    integer :: i, j, n, m
    real(dp) :: ll, rr, lr, scale

    n = split - first + 1
    m = last - split
    if (n < 2 .or. m < 2) then
      value = -huge(1.0_dp)
      return
    end if

    ll = 0.0_dp
    do j = first + 1, split
      do i = first, j - 1
        ll = ll + d(i, j)
      end do
    end do
    rr = 0.0_dp
    do j = split + 2, last
      do i = split + 1, j - 1
        rr = rr + d(i, j)
      end do
    end do
    lr = 0.0_dp
    do j = split + 1, last
      do i = first, split
        lr = lr + d(i, j)
      end do
    end do

    value = 2.0_dp*lr/real(n*m, dp)
    value = value - 2.0_dp*ll/real(n*(n - 1), dp)
    value = value - 2.0_dp*rr/real(m*(m - 1), dp)
    if (cp3o_scale) then
      scale = real(n*m, dp)/real((n + m)*(n + m), dp)
    else
      scale = real(n*m, dp)/real(n + m, dp)
    end if
    value = value*scale
  end function energy_split_stat

  pure subroutine split_point(d, first, last, min_size, split, statistic)
    real(dp), intent(in) :: d(:,:) !! Pairwise alpha-powered distance matrix for the full series.
    integer, intent(in) :: first !! First observation index of the interval to split.
    integer, intent(in) :: last !! Last observation index of the interval to split.
    integer, intent(in) :: min_size !! Minimum number of observations required in each child segment.
    integer, intent(out) :: split !! First observation index of the right child, or -1 if no split is feasible.
    real(dp), intent(out) :: statistic !! Largest divisive energy statistic among feasible splits.
    integer :: s
    real(dp) :: candidate

    split = -1
    statistic = -huge(1.0_dp)
    if (last - first + 1 < 2*min_size) return

    do s = first + min_size - 1, last - min_size
      candidate = energy_split_stat(d, first, s, last, .false.)
      if (candidate > statistic) then
        statistic = candidate
        split = s + 1
      end if
    end do
  end subroutine split_point

  pure real(dp) function ks_split_stat(x, first, split, last, window) result(value)
    real(dp), intent(in) :: x(:,:) !! Observations with rows indexing time; each column is compared marginally.
    integer, intent(in) :: first !! First observation index of the combined interval.
    integer, intent(in) :: split !! Last observation index of the left segment.
    integer, intent(in) :: last !! Last observation index of the right segment.
    integer, intent(in), optional :: window !! Optional maximum observations retained on each side near the split.
    integer :: left_first, right_last, j
    real(dp) :: candidate

    left_first = first
    right_last = last
    if (present(window)) then
      left_first = max(first, split - window + 1)
      right_last = min(last, split + window)
    end if

    value = 0.0_dp
    do j = 1, size(x, 2)
      candidate = ks_two_sample(x(left_first:split, j), x(split + 1:right_last, j))
      value = max(value, candidate)
    end do
  end function ks_split_stat

  pure real(dp) function ks_two_sample(x, y) result(value)
    real(dp), intent(in) :: x(:) !! First scalar sample.
    real(dp), intent(in) :: y(:) !! Second scalar sample.
    real(dp), allocatable :: xs(:), ys(:)
    integer :: i, j, n, m, ii, jj
    real(dp) :: diff, current, s

    n = size(x)
    m = size(y)
    if (n == 0 .or. m == 0) then
      value = 0.0_dp
      return
    end if
    xs = x
    ys = y
    call sort_real(xs)
    call sort_real(ys)
    i = 1
    j = 1
    current = 0.0_dp
    diff = 0.0_dp
    do while (i <= n .and. j <= m)
      if (xs(i) < ys(j)) then
        current = current + 1.0_dp/real(n, dp)
        i = i + 1
      else if (xs(i) > ys(j)) then
        current = current - 1.0_dp/real(m, dp)
        j = j + 1
      else
        s = xs(i)
        ii = i
        do while (ii <= n)
          if (xs(ii) > s) exit
          current = current + 1.0_dp/real(n, dp)
          ii = ii + 1
        end do
        jj = j
        do while (jj <= m)
          if (ys(jj) > s) exit
          current = current - 1.0_dp/real(m, dp)
          jj = jj + 1
        end do
        i = ii
        j = jj
      end if
      diff = max(diff, abs(current))
    end do
    do while (i <= n)
      current = current + 1.0_dp/real(n, dp)
      diff = max(diff, abs(current))
      i = i + 1
    end do
    do while (j <= m)
      current = current - 1.0_dp/real(m, dp)
      diff = max(diff, abs(current))
      j = j + 1
    end do
    value = diff*real(n*m, dp)/real((n + m)*(n + m), dp)
  end function ks_two_sample

  pure real(dp) function windowed_energy_split_stat(d, first, split, last, delta) result(value)
    real(dp), intent(in) :: d(:,:) !! Pairwise alpha-powered distance matrix for the full series.
    integer, intent(in) :: first !! First observation index of the combined interval.
    integer, intent(in) :: split !! Last observation index of the left segment.
    integer, intent(in) :: last !! Last observation index of the right segment.
    integer, intent(in) :: delta !! Complete window size used on each side of the split.
    integer :: i, j, n, m, lf, rr0, cll, crr, clr
    real(dp) :: dll, drr, dlr

    n = split - first + 1
    m = last - split
    if (n < delta + 1 .or. m < delta + 1) then
      value = -huge(1.0_dp)
      return
    end if
    lf = split - delta + 1
    rr0 = split + 1
    dll = 0.0_dp
    do j = lf + 1, split
      do i = lf, j - 1
        dll = dll + d(i, j)
      end do
    end do
    drr = 0.0_dp
    do j = rr0 + 1, split + delta
      do i = rr0, j - 1
        drr = drr + d(i, j)
      end do
    end do
    dlr = 0.0_dp
    do j = rr0, split + delta
      do i = lf, split
        dlr = dlr + d(i, j)
      end do
    end do

    do i = first, lf - 1
      dll = dll + d(i, i + 1)
    end do
    do i = split + delta, last - 1
      drr = drr + d(i, i + 1)
    end do
    cll = delta*(delta - 1)/2 + (n - delta)
    crr = delta*(delta - 1)/2 + (m - delta)
    clr = delta*delta
    value = 2.0_dp*dlr/real(clr, dp) - dll/real(cll, dp) - drr/real(crr, dp)
    value = value*real(n*m, dp)/real((n + m)*(n + m), dp)
  end function windowed_energy_split_stat

  pure subroutine sort_real(x)
    real(dp), intent(inout) :: x(:) !! Real vector sorted into ascending order in place.
    integer :: i, j
    real(dp) :: key

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
  end subroutine sort_real

end module ecp_energy
