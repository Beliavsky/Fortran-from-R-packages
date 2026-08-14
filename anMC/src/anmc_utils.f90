! SPDX-License-Identifier: GPL-3.0-only
module anmc_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_is_finite
  use anmc_kinds, only : dp, i8
  implicit none
  private
  public :: chronotime_ns, wall_time_seconds, positive_infinity, negative_infinity
  public :: mean_value, sample_variance, mean_finite, linear_fit, slope_through_origin
  public :: weighted_sample_without_replacement, weighted_sample_one, pairwise_distances
  public :: sort_indices_descending, complement_indices, gather_vector, gather_matrix
  public :: seed_fortran_rng, clamp_probability, normal_cdf_local

contains

  real(dp) function chronotime_ns() result(t)
    integer(i8) :: count, rate
    call system_clock(count=count, count_rate=rate)
    if (rate > 0_i8) then
      t = real(count, dp) * 1.0e9_dp / real(rate, dp)
    else
      t = 0.0_dp
    end if
  end function chronotime_ns

  real(dp) function wall_time_seconds() result(t)
    t = chronotime_ns() * 1.0e-9_dp
  end function wall_time_seconds

  real(dp) function positive_infinity() result(x)
    x = ieee_value(0.0_dp, ieee_positive_inf)
  end function positive_infinity

  real(dp) function negative_infinity() result(x)
    x = ieee_value(0.0_dp, ieee_negative_inf)
  end function negative_infinity

  real(dp) function mean_value(x) result(v)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      v = 0.0_dp
    else
      v = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  real(dp) function sample_variance(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    integer :: n
    n = size(x)
    if (n <= 1) then
      v = 0.0_dp
      return
    end if
    m = sum(x) / real(n, dp)
    v = sum((x - m)**2) / real(n - 1, dp)
  end function sample_variance

  real(dp) function mean_finite(x, fallback) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: fallback
    integer :: i, n
    real(dp) :: s, fb
    fb = 0.0_dp
    if (present(fallback)) fb = fallback
    s = 0.0_dp
    n = 0
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        s = s + x(i)
        n = n + 1
      end if
    end do
    if (n > 0) then
      v = s / real(n, dp)
    else
      v = fb
    end if
  end function mean_finite

  subroutine linear_fit(x, y, intercept, slope)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: intercept, slope
    real(dp) :: xm, ym, den
    if (size(x) /= size(y) .or. size(x) == 0) then
      intercept = 0.0_dp
      slope = 0.0_dp
      return
    end if
    xm = mean_value(x)
    ym = mean_value(y)
    den = sum((x - xm)**2)
    if (den <= tiny(1.0_dp)) then
      slope = 0.0_dp
      intercept = ym
    else
      slope = sum((x - xm) * (y - ym)) / den
      intercept = ym - slope * xm
    end if
  end subroutine linear_fit

  real(dp) function slope_through_origin(x, y) result(slope)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: den
    if (size(x) /= size(y) .or. size(x) == 0) then
      slope = 0.0_dp
      return
    end if
    den = sum(x*x)
    if (den <= tiny(1.0_dp)) then
      slope = 0.0_dp
    else
      slope = sum(x*y) / den
    end if
  end function slope_through_origin

  subroutine seed_fortran_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1) + 1
    end do
    call random_seed(put=put)
  end subroutine seed_fortran_rng

  integer function weighted_sample_one(weights) result(idx)
    real(dp), intent(in) :: weights(:)
    real(dp) :: total, u, accum
    integer :: i
    total = sum(max(weights, 0.0_dp))
    if (total <= 0.0_dp) then
      idx = 0
      return
    end if
    call random_number(u)
    u = u * total
    accum = 0.0_dp
    idx = size(weights)
    do i = 1, size(weights)
      accum = accum + max(weights(i), 0.0_dp)
      if (u <= accum) then
        idx = i
        return
      end if
    end do
  end function weighted_sample_one

  function weighted_sample_without_replacement(weights, q) result(indices)
    real(dp), intent(in) :: weights(:)
    integer, intent(in) :: q
    integer, allocatable :: indices(:)
    real(dp), allocatable :: w(:)
    integer :: i, k
    if (q < 0 .or. q > size(weights)) then
      allocate(indices(0))
      return
    end if
    allocate(indices(q), w(size(weights)))
    w = max(weights, 0.0_dp)
    do i = 1, q
      k = weighted_sample_one(w)
      if (k == 0) then
        ! Fall back to uniform selection among remaining positions.
        do k = 1, size(w)
          if (w(k) >= 0.0_dp .and. .not. any(indices(1:max(0,i-1)) == k)) exit
        end do
      end if
      indices(i) = k
      w(k) = 0.0_dp
      ! Distinguish already selected zero-weight items from unselected zeros.
      if (sum(w) <= 0.0_dp .and. i < q) then
        w = 1.0_dp
        do k = 1, i
          w(indices(k)) = 0.0_dp
        end do
      end if
    end do
  end function weighted_sample_without_replacement

  function pairwise_distances(e) result(d)
    real(dp), intent(in) :: e(:,:)
    real(dp), allocatable :: d(:,:)
    integer :: i, j, n
    n = size(e,1)
    allocate(d(n,n))
    do j = 1, n
      do i = 1, n
        d(i,j) = sqrt(sum((e(i,:) - e(j,:))**2))
      end do
    end do
  end function pairwise_distances

  function sort_indices_descending(x) result(idx)
    real(dp), intent(in) :: x(:)
    integer, allocatable :: idx(:)
    integer :: i, j, key
    allocate(idx(size(x)))
    idx = [(i, i=1,size(x))]
    do i = 2, size(x)
      key = idx(i)
      j = i - 1
      do while (j >= 1)
        if (x(idx(j)) >= x(key)) exit
        idx(j+1) = idx(j)
        j = j - 1
      end do
      idx(j+1) = key
    end do
  end function sort_indices_descending

  function complement_indices(n, selected) result(comp)
    integer, intent(in) :: n, selected(:)
    integer, allocatable :: comp(:)
    logical, allocatable :: mark(:)
    integer :: i, k
    allocate(mark(n))
    mark = .false.
    do i = 1, size(selected)
      if (selected(i) >= 1 .and. selected(i) <= n) mark(selected(i)) = .true.
    end do
    allocate(comp(count(.not. mark)))
    k = 0
    do i = 1, n
      if (.not. mark(i)) then
        k = k + 1
        comp(k) = i
      end if
    end do
  end function complement_indices

  function gather_vector(x, idx) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: idx(:)
    real(dp), allocatable :: y(:)
    integer :: i
    allocate(y(size(idx)))
    do i = 1, size(idx)
      y(i) = x(idx(i))
    end do
  end function gather_vector

  function gather_matrix(a, rows, cols) result(b)
    real(dp), intent(in) :: a(:,:)
    integer, intent(in) :: rows(:), cols(:)
    real(dp), allocatable :: b(:,:)
    integer :: i, j
    allocate(b(size(rows), size(cols)))
    do j = 1, size(cols)
      do i = 1, size(rows)
        b(i,j) = a(rows(i), cols(j))
      end do
    end do
  end function gather_matrix

  pure real(dp) function clamp_probability(p) result(v)
    real(dp), intent(in) :: p
    v = min(1.0_dp, max(0.0_dp, p))
  end function clamp_probability

  pure real(dp) function normal_cdf_local(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf_local

end module anmc_utils
