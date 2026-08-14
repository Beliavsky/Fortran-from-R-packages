module dice_design_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  use dice_design_kinds, only : dp
  implicit none
  private

  public :: euclidean_distance, pairwise_distances, rescale_unit_cube
  public :: determinant, sort_real, empirical_cdf_column, empirical_quantile
  public :: normal_quantile, mean_real, sample_sd

contains

  pure function euclidean_distance(x, y) result(d)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: d
    d = sqrt(sum((x - y)**2))
  end function euclidean_distance

  subroutine pairwise_distances(x, dmat)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: dmat(:, :)
    integer :: n, i, j

    n = size(x, 1)
    allocate(dmat(n, n))
    dmat = 0.0_dp
    do j = 2, n
      do i = 1, j - 1
        dmat(i, j) = euclidean_distance(x(i, :), x(j, :))
        dmat(j, i) = dmat(i, j)
      end do
    end do
  end subroutine pairwise_distances

  subroutine rescale_unit_cube(x, y, changed)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: y(:, :)
    logical, intent(out) :: changed
    integer :: j
    real(dp) :: xmin, xmax

    allocate(y(size(x, 1), size(x, 2)))
    y = x
    changed = minval(x) < 0.0_dp .or. maxval(x) > 1.0_dp
    if (.not. changed) return

    do j = 1, size(x, 2)
      xmin = minval(x(:, j))
      xmax = maxval(x(:, j))
      if (xmax > xmin) then
        y(:, j) = (x(:, j) - xmin) / (xmax - xmin)
      else
        y(:, j) = 0.5_dp
      end if
    end do
  end subroutine rescale_unit_cube

  function determinant(a) result(det)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: det
    real(dp), allocatable :: b(:, :)
    real(dp) :: pivot, factor
    integer :: n, i, j, k, p
    real(dp) :: maxabs

    n = size(a, 1)
    if (size(a, 2) /= n) error stop 'determinant: matrix must be square'
    allocate(b(n, n))
    b = a
    det = 1.0_dp
    do k = 1, n
      p = k
      maxabs = abs(b(k, k))
      do i = k + 1, n
        if (abs(b(i, k)) > maxabs) then
          p = i
          maxabs = abs(b(i, k))
        end if
      end do
      if (maxabs <= tiny(1.0_dp)) then
        det = 0.0_dp
        return
      end if
      if (p /= k) then
        do j = 1, n
          pivot = b(k, j)
          b(k, j) = b(p, j)
          b(p, j) = pivot
        end do
        det = -det
      end if
      pivot = b(k, k)
      det = det * pivot
      do i = k + 1, n
        factor = b(i, k) / pivot
        b(i, k) = 0.0_dp
        b(i, k + 1:n) = b(i, k + 1:n) - factor * b(k, k + 1:n)
      end do
    end do
  end function determinant

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
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

  subroutine empirical_cdf_column(x, u)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: u(:)
    integer :: n, i, j, count_le

    n = size(x)
    if (size(u) /= n) error stop 'empirical_cdf_column: size mismatch'
    do i = 1, n
      count_le = 0
      do j = 1, n
        if (x(j) <= x(i)) count_le = count_le + 1
      end do
      u(i) = real(count_le, dp) / real(n, dp)
    end do
  end subroutine empirical_cdf_column

  function empirical_quantile(x, p) result(q)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: p
    real(dp) :: q
    real(dp), allocatable :: s(:)
    real(dp) :: h, frac
    integer :: n, lo

    n = size(x)
    if (n == 0) then
      q = ieee_value(0.0_dp, ieee_positive_inf)
      return
    end if
    allocate(s(n))
    s = x
    call sort_real(s)
    if (p <= 0.0_dp) then
      q = s(1)
    else if (p >= 1.0_dp) then
      q = s(n)
    else
      h = 1.0_dp + real(n - 1, dp) * p
      lo = int(floor(h))
      frac = h - real(lo, dp)
      if (lo >= n) then
        q = s(n)
      else
        q = (1.0_dp - frac) * s(lo) + frac * s(lo + 1)
      end if
    end if
  end function empirical_quantile

  function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a1 = -3.969683028665376e1_dp
    real(dp), parameter :: a2 = 2.209460984245205e2_dp
    real(dp), parameter :: a3 = -2.759285104469687e2_dp
    real(dp), parameter :: a4 = 1.383577518672690e2_dp
    real(dp), parameter :: a5 = -3.066479806614716e1_dp
    real(dp), parameter :: a6 = 2.506628277459239_dp
    real(dp), parameter :: b1 = -5.447609879822406e1_dp
    real(dp), parameter :: b2 = 1.615858368580409e2_dp
    real(dp), parameter :: b3 = -1.556989798598866e2_dp
    real(dp), parameter :: b4 = 6.680131188771972e1_dp
    real(dp), parameter :: b5 = -1.328068155288572e1_dp
    real(dp), parameter :: c1 = -7.784894002430293e-3_dp
    real(dp), parameter :: c2 = -3.223964580411365e-1_dp
    real(dp), parameter :: c3 = -2.400758277161838_dp
    real(dp), parameter :: c4 = -2.549732539343734_dp
    real(dp), parameter :: c5 = 4.374664141464968_dp
    real(dp), parameter :: c6 = 2.938163982698783_dp
    real(dp), parameter :: d1 = 7.784695709041462e-3_dp
    real(dp), parameter :: d2 = 3.224671290700398e-1_dp
    real(dp), parameter :: d3 = 2.445134137142996_dp
    real(dp), parameter :: d4 = 3.754408661907416_dp
    real(dp), parameter :: plow = 0.02425_dp
    real(dp), parameter :: phigh = 1.0_dp - plow

    if (p <= 0.0_dp) then
      x = -ieee_value(0.0_dp, ieee_positive_inf)
    else if (p >= 1.0_dp) then
      x = ieee_value(0.0_dp, ieee_positive_inf)
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
          ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
          (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
           ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
  end function normal_quantile

  pure function mean_real(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) == 0) then
      m = 0.0_dp
    else
      m = sum(x) / real(size(x), dp)
    end if
  end function mean_real

  pure function sample_sd(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: s, m
    if (size(x) <= 1) then
      s = 0.0_dp
    else
      m = sum(x) / real(size(x), dp)
      s = sqrt(sum((x - m)**2) / real(size(x) - 1, dp))
    end if
  end function sample_sd

end module dice_design_utils
