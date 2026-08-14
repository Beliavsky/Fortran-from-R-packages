module mlr_utils
  use mlr_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: mean_dp, variance_dp, median_dp, ranks_average, solve_linear, normal_cdf
  public :: sigmoid, log1pexp, is_finite_vector, argsort_real, unique_int_count
contains
  pure real(dp) function mean_dp(x) result(v)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      v = sum(x) / real(size(x), dp)
    end if
  end function mean_dp

  pure real(dp) function variance_dp(x, sample) result(v)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: sample
    real(dp) :: m, den
    logical :: s
    s = .true.
    if (present(sample)) s = sample
    if (size(x) == 0 .or. (s .and. size(x) < 2)) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    m = mean_dp(x)
    den = real(size(x), dp)
    if (s) den = den - 1.0_dp
    v = sum((x - m)**2) / den
  end function variance_dp

  real(dp) function median_dp(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    integer :: i, j, n
    real(dp) :: key
    n = size(x)
    if (n == 0) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    allocate(y(n)); y = x
    do i = 2, n
      key = y(i); j = i - 1
      do while (j >= 1)
        if (y(j) <= key) exit
        y(j+1) = y(j); j = j - 1
      end do
      y(j+1) = key
    end do
    if (mod(n,2) == 1) then
      v = y((n+1)/2)
    else
      v = 0.5_dp * (y(n/2) + y(n/2+1))
    end if
  end function median_dp

  subroutine argsort_real(x, idx)
    real(dp), intent(in) :: x(:)
    integer, allocatable, intent(out) :: idx(:)
    integer :: i, j, key
    allocate(idx(size(x)))
    idx = [(i, i=1,size(x))]
    do i = 2, size(idx)
      key = idx(i); j = i - 1
      do while (j >= 1)
        if (x(idx(j)) <= x(key)) exit
        idx(j+1) = idx(j); j = j - 1
      end do
      idx(j+1) = key
    end do
  end subroutine argsort_real

  subroutine ranks_average(x, r)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: r(:)
    integer, allocatable :: idx(:)
    integer :: i, j, k, n
    real(dp) :: rr, tol
    n = size(x)
    if (size(r) /= n) error stop "ranks_average: size mismatch"
    call argsort_real(x, idx)
    i = 1
    do while (i <= n)
      j = i
      tol = 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(idx(i))))
      do while (j < n)
        if (abs(x(idx(j+1))-x(idx(i))) > tol) exit
        j = j + 1
      end do
      rr = 0.5_dp * real(i+j, dp)
      do k = i, j
        r(idx(k)) = rr
      end do
      i = j + 1
    end do
  end subroutine ranks_average

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: pivot, factor, tmp
    integer :: n, i, j, k, p, stat
    n = size(b)
    stat = 0
    if (size(a,1) /= n .or. size(a,2) /= n) then
      if (present(info)) info = -1
      allocate(x(0)); return
    end if
    allocate(aa(n,n), bb(n), x(n)); aa = a; bb = b
    do k = 1, n-1
      p = k
      do i = k+1, n
        if (abs(aa(i,k)) > abs(aa(p,k))) p = i
      end do
      if (abs(aa(p,k)) <= tiny(1.0_dp)) then
        stat = k; exit
      end if
      if (p /= k) then
        do j = k, n
          tmp = aa(k,j); aa(k,j) = aa(p,j); aa(p,j) = tmp
        end do
        tmp = bb(k); bb(k) = bb(p); bb(p) = tmp
      end if
      do i = k+1, n
        factor = aa(i,k) / aa(k,k)
        aa(i,k:n) = aa(i,k:n) - factor * aa(k,k:n)
        bb(i) = bb(i) - factor * bb(k)
      end do
    end do
    if (stat == 0) then
      if (abs(aa(n,n)) <= tiny(1.0_dp)) stat = n
    end if
    if (stat /= 0) then
      x = 0.0_dp
      if (present(info)) info = stat
      return
    end if
    do i = n, 1, -1
      pivot = bb(i)
      if (i < n) pivot = pivot - dot_product(aa(i,i+1:n), x(i+1:n))
      x(i) = pivot / aa(i,i)
    end do
    if (present(info)) info = 0
  end subroutine solve_linear

  pure elemental real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure elemental real(dp) function sigmoid(x) result(p)
    real(dp), intent(in) :: x
    if (x >= 0.0_dp) then
      p = 1.0_dp / (1.0_dp + exp(-x))
    else
      p = exp(x) / (1.0_dp + exp(x))
    end if
  end function sigmoid

  pure elemental real(dp) function log1pexp(x) result(v)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      v = x + log(1.0_dp + exp(-x))
    else
      v = log(1.0_dp + exp(x))
    end if
  end function log1pexp

  pure logical function is_finite_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i
    ok = .true.
    do i = 1, size(x)
      if (.not. ieee_is_finite_local(x(i))) then
        ok = .false.; return
      end if
    end do
  contains
    pure logical function ieee_is_finite_local(v)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: v
      ieee_is_finite_local = ieee_is_finite(v)
    end function ieee_is_finite_local
  end function is_finite_vector

  pure integer function unique_int_count(x) result(nu)
    integer, intent(in) :: x(:)
    integer :: i, j
    logical :: seen
    nu = 0
    do i = 1, size(x)
      seen = .false.
      do j = 1, i-1
        if (x(j) == x(i)) then
          seen = .true.; exit
        end if
      end do
      if (.not. seen) nu = nu + 1
    end do
  end function unique_int_count
end module mlr_utils
