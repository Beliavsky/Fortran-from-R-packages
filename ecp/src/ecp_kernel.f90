module ecp_kernel
  use ecp_kinds, only: dp
  implicit none
  private

  public :: kcpa

contains

  pure function kcpa(x, l, c, bandwidth_rows) result(cps)
    real(dp), intent(in) :: x(:,:) !! Time series matrix with observations in rows.
    integer, intent(in) :: l !! Maximum number of change points.
    real(dp), intent(in) :: c !! Model-selection penalty constant.
    integer, intent(in), optional :: bandwidth_rows(:) !! Optional 1-based rows used to estimate kernel bandwidth.
    integer, allocatable :: cps(:)
    real(dp), allocatable :: kernel(:,:), vcost(:,:), ii(:,:), score(:)
    integer, allocatable :: h(:,:)
    real(dp) :: bandwidth, vmax, d2, tmp
    integer :: n, max_segments, seg, i, j, best_seg, current

    n = size(x, 1)
    if (n == 0) then
      allocate(cps(0))
      return
    end if
    max_segments = min(max(1, l + 1), n)
    bandwidth = kernel_bandwidth(x, bandwidth_rows)
    bandwidth = max(bandwidth, sqrt(tiny(1.0_dp)))
    allocate(kernel(n, n))
    do j = 1, n
      do i = 1, n
        d2 = sum((x(i, :) - x(j, :))**2)
        kernel(i, j) = exp(-d2/bandwidth)
      end do
    end do
    vcost = segment_costs(kernel)
    vmax = get_vmax(x)

    allocate(ii(max_segments, n), source=huge(1.0_dp))
    allocate(h(max_segments, n), source=0)
    ii(1, :) = vcost(1, :)
    do seg = 2, max_segments
      do i = seg, n
        do j = seg - 1, i - 1
          tmp = ii(seg - 1, j) + vcost(j + 1, i)
          if (tmp < ii(seg, i)) then
            ii(seg, i) = tmp
            h(seg, i) = j + 1
          end if
        end do
      end do
    end do

    allocate(score(max_segments))
    do seg = 1, max_segments
      score(seg) = ii(seg, n) + c*vmax*real(seg, dp)/real(n, dp)* &
        (1.0_dp + log(real(n, dp)/real(seg, dp)))
    end do
    best_seg = minloc(score, dim=1)
    allocate(cps(best_seg + 1), source=0)
    cps(1) = 1
    cps(best_seg + 1) = n + 1
    current = n
    do seg = best_seg, 2, -1
      cps(seg) = h(seg, current)
      current = cps(seg) - 1
    end do
  end function kcpa

  pure real(dp) function kernel_bandwidth(x, rows) result(value)
    real(dp), intent(in) :: x(:,:) !! Observation matrix used to form pairwise squared distances.
    integer, intent(in), optional :: rows(:) !! Optional 1-based subset of rows used in the bandwidth median.
    real(dp), allocatable :: distances(:)
    integer, allocatable :: use_rows(:)
    integer :: i, j, n, nr, idx

    n = size(x, 1)
    if (present(rows)) then
      if (size(rows) > 0 .and. all(rows >= 1) .and. all(rows <= n)) then
        use_rows = rows
      else
        allocate(use_rows(n))
        use_rows = [(i, i=1, n)]
      end if
    else
      allocate(use_rows(n))
      use_rows = [(i, i=1, n)]
    end if
    nr = size(use_rows)
    allocate(distances(nr*nr))
    idx = 0
    do i = 1, nr
      do j = 1, nr
        idx = idx + 1
        distances(idx) = sum((x(use_rows(i), :) - x(use_rows(j), :))**2)
      end do
    end do
    call sort_real(distances)
    if (mod(size(distances), 2) == 1) then
      value = distances((size(distances) + 1)/2)
    else
      value = 0.5_dp*(distances(size(distances)/2) + distances(size(distances)/2 + 1))
    end if
  end function kernel_bandwidth

  pure function segment_costs(kernel) result(v)
    real(dp), intent(in) :: kernel(:,:) !! Symmetric kernel matrix.
    real(dp), allocatable :: v(:,:)
    real(dp), allocatable :: prefix(:,:), diag_prefix(:)
    real(dp) :: diag_sum, total
    integer :: i, j, n, len

    n = size(kernel, 1)
    allocate(prefix(0:n, 0:n), source=0.0_dp)
    do j = 1, n
      do i = 1, n
        prefix(i, j) = kernel(i, j) + prefix(i - 1, j) + prefix(i, j - 1) - prefix(i - 1, j - 1)
      end do
    end do
    allocate(diag_prefix(0:n), source=0.0_dp)
    do i = 1, n
      diag_prefix(i) = diag_prefix(i - 1) + kernel(i, i)
    end do

    allocate(v(n, n), source=0.0_dp)
    do i = 1, n
      do j = i, n
        diag_sum = diag_prefix(j) - diag_prefix(i - 1)
        total = block_sum(prefix, i, j, i, j)
        len = j - i + 1
        v(i, j) = diag_sum - total/real(len, dp)
        v(j, i) = v(i, j)
      end do
    end do
  end function segment_costs

  pure real(dp) function block_sum(prefix, row_first, row_last, col_first, col_last) result(value)
    real(dp), intent(in) :: prefix(0:,0:) !! Two-dimensional prefix table with a zero lower border.
    integer, intent(in) :: row_first !! First row index of the requested block.
    integer, intent(in) :: row_last !! Last row index of the requested block.
    integer, intent(in) :: col_first !! First column index of the requested block.
    integer, intent(in) :: col_last !! Last column index of the requested block.

    value = prefix(row_last, col_last) - prefix(row_first - 1, col_last)
    value = value - prefix(row_last, col_first - 1) + prefix(row_first - 1, col_first - 1)
  end function block_sum

  pure real(dp) function get_vmax(x) result(value)
    real(dp), intent(in) :: x(:,:) !! Observation matrix used by the KCPA penalty bound.
    integer :: n, lower, upper
    real(dp) :: a, b

    n = size(x, 1)
    lower = ceiling(0.05_dp*real(n, dp))
    upper = floor(0.95_dp*real(n, dp))
    lower = max(1, min(lower, n))
    upper = max(1, min(upper, n))
    if (lower == 1) then
      a = 1.0_dp
    else
      a = trace_sample_cov(x(1:lower, :))
    end if
    b = trace_sample_cov(x(upper:n, :))
    value = max(a, b)
  end function get_vmax

  pure real(dp) function trace_sample_cov(x) result(value)
    real(dp), intent(in) :: x(:,:) !! Sample matrix whose covariance trace is requested.
    real(dp), allocatable :: mean(:)
    integer :: i, n

    n = size(x, 1)
    if (n <= 1) then
      value = 0.0_dp
      return
    end if
    mean = sum(x, dim=1)/real(n, dp)
    value = 0.0_dp
    do i = 1, n
      value = value + sum((x(i, :) - mean)**2)
    end do
    value = value/real(n - 1, dp)
  end function trace_sample_cov

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

end module ecp_kernel
