module trading_stats
  use trading_kinds, only : dp
  implicit none
  private

  public :: mean_dp
  public :: variance_dp
  public :: sd_dp
  public :: covariance_dp
  public :: correlation_dp
  public :: correlation_matrix
  public :: quantile_type7
  public :: normal_cdf
  public :: sort_real
  public :: set_random_seed
  public :: random_normal

contains

  pure real(dp) function mean_dp(x) result(value)
    real(dp), intent(in) :: x(:)

    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_dp

  pure real(dp) function variance_dp(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: center

    if (size(x) < 2) then
      value = 0.0_dp
      return
    end if

    center = mean_dp(x)
    value = sum((x - center)**2) / real(size(x) - 1, dp)
  end function variance_dp

  pure real(dp) function sd_dp(x) result(value)
    real(dp), intent(in) :: x(:)

    value = sqrt(max(variance_dp(x), 0.0_dp))
  end function sd_dp

  pure real(dp) function covariance_dp(x, y) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp) :: mean_x
    real(dp) :: mean_y
    integer :: n

    n = min(size(x), size(y))
    if (n < 2) then
      value = 0.0_dp
      return
    end if

    mean_x = sum(x(:n)) / real(n, dp)
    mean_y = sum(y(:n)) / real(n, dp)
    value = sum((x(:n) - mean_x) * (y(:n) - mean_y)) / real(n - 1, dp)
  end function covariance_dp

  pure real(dp) function correlation_dp(x, y) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp) :: denom

    denom = sd_dp(x) * sd_dp(y)
    if (denom <= tiny(1.0_dp)) then
      value = 0.0_dp
    else
      value = covariance_dp(x, y) / denom
      value = max(-1.0_dp, min(1.0_dp, value))
    end if
  end function correlation_dp

  subroutine correlation_matrix(x, corr)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: corr(:, :)
    integer :: i
    integer :: j
    integer :: p

    p = size(x, 2)
    if (size(corr, 1) /= p .or. size(corr, 2) /= p) then
      error stop "correlation_matrix: result has the wrong shape"
    end if

    do j = 1, p
      do i = 1, p
        corr(i, j) = correlation_dp(x(:, i), x(:, j))
      end do
    end do
  end subroutine correlation_matrix

  real(dp) function quantile_type7(x, probability) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: probability
    real(dp), allocatable :: work(:)
    real(dp) :: h
    real(dp) :: fraction
    integer :: lower
    integer :: n

    n = size(x)
    if (n == 0) error stop "quantile_type7: empty input"
    if (probability < 0.0_dp .or. probability > 1.0_dp) then
      error stop "quantile_type7: probability must be in [0, 1]"
    end if

    allocate(work(n))
    work = x
    call sort_real(work)

    if (n == 1) then
      value = work(1)
      return
    end if

    h = 1.0_dp + real(n - 1, dp) * probability
    lower = floor(h)
    fraction = h - real(lower, dp)

    if (lower >= n) then
      value = work(n)
    else
      value = (1.0_dp - fraction) * work(lower) + fraction * work(lower + 1)
    end if
  end function quantile_type7

  pure real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x

    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)

    if (size(x) > 1) call quicksort_real(x, 1, size(x))
  end subroutine sort_real

  recursive subroutine quicksort_real(x, left, right)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: left
    integer, intent(in) :: right
    real(dp) :: pivot
    real(dp) :: temp
    integer :: i
    integer :: j

    i = left
    j = right
    pivot = x((left + right) / 2)

    do
      do while (x(i) < pivot)
        i = i + 1
      end do
      do while (x(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        temp = x(i)
        x(i) = x(j)
        x(j) = temp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do

    if (left < j) call quicksort_real(x, left, j)
    if (i < right) call quicksort_real(x, i, right)
  end subroutine quicksort_real

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer, allocatable :: seed_values(:)
    integer :: i
    integer :: n

    call random_seed(size=n)
    allocate(seed_values(n))
    do i = 1, n
      seed_values(i) = modulo(seed + 104729 * i, huge(1) - 1)
      if (seed_values(i) <= 0) seed_values(i) = i
    end do
    call random_seed(put=seed_values)
  end subroutine set_random_seed

  subroutine random_normal(x)
    real(dp), intent(out) :: x(:)
    real(dp) :: u1
    real(dp) :: u2
    real(dp) :: radius
    real(dp) :: angle
    integer :: i

    i = 1
    do while (i <= size(x))
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      radius = sqrt(-2.0_dp * log(u1))
      angle = 2.0_dp * acos(-1.0_dp) * u2
      x(i) = radius * cos(angle)
      if (i + 1 <= size(x)) x(i + 1) = radius * sin(angle)
      i = i + 2
    end do
  end subroutine random_normal

end module trading_stats
