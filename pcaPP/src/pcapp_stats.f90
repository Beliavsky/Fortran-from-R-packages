module pcapp_stats
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use pcapp_kinds, only: dp
  implicit none
  private

  real(dp), parameter :: qn_asymptotic = 2.219144465985076_dp
  real(dp), parameter :: mad_constant = 1.482602218505602_dp

  public :: mean_vec, sample_sd, median_vec, mad_scale, qn_scale, robust_scale
  public :: column_center, column_scale, cor_fk_vector, cor_fk_matrix
  public :: spatial_objective, spatial_gradient, spatial_hessian
  public :: symmetric_eigen_jacobi, sort_eigen_desc, orient_columns
  public :: random_normal_vector, random_normal_matrix

contains

  pure recursive subroutine quicksort_real(a, left, right)
    real(dp), intent(inout) :: a(:) !! Array sorted in ascending order in place.
    integer, intent(in) :: left !! First one-based element index in the active partition.
    integer, intent(in) :: right !! Last one-based element index in the active partition.
    integer :: i, j
    real(dp) :: pivot, temp

    if (left >= right) return
    i = left
    j = right
    pivot = a((left + right) / 2)
    do
      do while (a(i) < pivot)
        i = i + 1
      end do
      do while (a(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        temp = a(i)
        a(i) = a(j)
        a(j) = temp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call quicksort_real(a, left, j)
    if (i < right) call quicksort_real(a, i, right)
  end subroutine quicksort_real

  pure subroutine sort_real(a)
    real(dp), intent(inout) :: a(:) !! Array sorted in ascending order in place.

    if (size(a) > 1) call quicksort_real(a, 1, size(a))
  end subroutine sort_real

  pure real(dp) function mean_vec(x) result(value)
    real(dp), intent(in) :: x(:) !! Numeric observations whose arithmetic mean is requested.

    if (size(x) == 0) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_vec

  pure real(dp) function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:) !! Numeric observations; at least two are required for sample SD.
    real(dp) :: mu

    if (size(x) < 2) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    mu = mean_vec(x)
    value = sqrt(sum((x - mu)**2) / real(size(x) - 1, dp))
  end function sample_sd

  pure real(dp) function median_vec(x) result(value)
    real(dp), intent(in) :: x(:) !! Numeric observations whose sample median is requested.
    real(dp), allocatable :: work(:)
    integer :: n

    n = size(x)
    if (n == 0) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    work = x
    call sort_real(work)
    if (mod(n, 2) == 1) then
      value = work((n + 1) / 2)
    else
      value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
    end if
  end function median_vec

  pure real(dp) function mad_scale(x) result(value)
    real(dp), intent(in) :: x(:) !! Numeric observations for the median absolute deviation scale.
    real(dp), allocatable :: dev(:)
    real(dp) :: center

    if (size(x) == 0) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    center = median_vec(x)
    dev = abs(x - center)
    value = mad_constant * median_vec(dev)
  end function mad_scale

  pure real(dp) function qn_finite_factor(n) result(value)
    integer, intent(in) :: n !! Sample size used for the pcaPP finite-sample Qn correction.
    real(dp), parameter :: small(2:9) = [ &
      0.400_dp, 0.993_dp, 0.514_dp, 0.845_dp, &
      0.612_dp, 0.859_dp, 0.670_dp, 0.874_dp]

    if (n < 2) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (n <= 9) then
      value = qn_asymptotic * small(n)
    else if (mod(n, 2) == 1) then
      value = qn_asymptotic * real(n, dp) / (real(n, dp) + 1.4_dp)
    else
      value = qn_asymptotic * real(n, dp) / (real(n, dp) + 3.8_dp)
    end if
  end function qn_finite_factor

  pure real(dp) function qn_scale(x, corr_fact) result(value)
    real(dp), intent(in) :: x(:) !! Numeric sample for the robust Qn scale estimator.
    real(dp), intent(in), optional :: corr_fact !! Optional asymptotic normalization replacing pcaPP's default constant.
    real(dp), allocatable :: diffs(:)
    real(dp) :: factor
    integer :: h, i, j, k, n, pos

    n = size(x)
    if (n < 2) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    allocate(diffs(n * (n - 1) / 2))
    pos = 0
    do i = 1, n - 1
      do j = i + 1, n
        pos = pos + 1
        diffs(pos) = abs(x(i) - x(j))
      end do
    end do
    call sort_real(diffs)
    h = n / 2 + 1
    k = h * (h - 1) / 2
    if (present(corr_fact)) then
      factor = corr_fact
      if (n <= 9) then
        factor = factor * qn_finite_factor(n) / qn_asymptotic
      else if (mod(n, 2) == 1) then
        factor = factor * real(n, dp) / (real(n, dp) + 1.4_dp)
      else
        factor = factor * real(n, dp) / (real(n, dp) + 3.8_dp)
      end if
    else
      factor = qn_finite_factor(n)
    end if
    value = factor * diffs(k)
  end function qn_scale

  pure real(dp) function robust_scale(x, method) result(value)
    real(dp), intent(in) :: x(:) !! Projected observations whose scatter is estimated.
    integer, intent(in) :: method !! Scale selector: 0=SD, 1=MAD, 2=Qn, 5=root second moment.

    select case (method)
    case (0)
      value = sample_sd(x)
    case (1)
      value = mad_scale(x)
    case (2)
      value = qn_scale(x)
    case (5)
      if (size(x) == 0) then
        value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
        value = sqrt(sum(x**2) / real(size(x), dp))
      end if
    case default
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    end select
  end function robust_scale

  pure subroutine column_center(x, method, center)
    real(dp), intent(in) :: x(:,:) !! Data matrix with observations in rows and variables in columns.
    integer, intent(in) :: method !! Center selector: 0=arithmetic mean, 1=coordinatewise median.
    real(dp), intent(out) :: center(size(x, 2)) !! Per-variable center estimates.
    integer :: j

    do j = 1, size(x, 2)
      select case (method)
      case (0)
        center(j) = mean_vec(x(:, j))
      case (1)
        center(j) = median_vec(x(:, j))
      case default
        center(j) = 0.0_dp
      end select
    end do
  end subroutine column_center

  pure subroutine column_scale(x, method, scales)
    real(dp), intent(in) :: x(:,:) !! Data matrix with observations in rows and variables in columns.
    integer, intent(in) :: method !! Scale selector: -1=none, 0=SD, 1=MAD, 2=Qn.
    real(dp), intent(out) :: scales(size(x, 2)) !! Per-variable scale estimates.
    integer :: j

    if (method < 0) then
      scales = 1.0_dp
      return
    end if
    do j = 1, size(x, 2)
      scales(j) = robust_scale(x(:, j), method)
      if (abs(scales(j)) <= tiny(1.0_dp)) scales(j) = 1.0_dp
    end do
  end subroutine column_scale

  pure real(dp) function cor_fk_vector(x, y) result(value)
    real(dp), intent(in) :: x(:) !! First numeric vector for Kendall tau-b correlation.
    real(dp), intent(in) :: y(:) !! Second numeric vector, with the same length as x.
    integer :: i, j
    integer(kind=8) :: concord, discord, ties_x, ties_y, npair
    real(dp) :: d1, d2

    if (size(x) /= size(y) .or. size(x) < 2) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    concord = 0_8
    discord = 0_8
    ties_x = 0_8
    ties_y = 0_8
    do i = 1, size(x) - 1
      do j = i + 1, size(x)
        if (x(i) <= x(j) .and. x(i) >= x(j)) ties_x = ties_x + 1_8
        if (y(i) <= y(j) .and. y(i) >= y(j)) ties_y = ties_y + 1_8
        if ((x(i) <= x(j) .and. x(i) >= x(j)) .or. &
            (y(i) <= y(j) .and. y(i) >= y(j))) cycle
        if ((x(i) - x(j)) * (y(i) - y(j)) > 0.0_dp) then
          concord = concord + 1_8
        else
          discord = discord + 1_8
        end if
      end do
    end do
    npair = int(size(x), 8) * int(size(x) - 1, 8) / 2_8
    d1 = real(npair - ties_x, dp)
    d2 = real(npair - ties_y, dp)
    if (d1 <= 0.0_dp .or. d2 <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = real(concord - discord, dp) / sqrt(d1 * d2)
    end if
  end function cor_fk_vector

  pure subroutine cor_fk_matrix(x, corr)
    real(dp), intent(in) :: x(:,:) !! Data matrix with observations in rows and variables in columns.
    real(dp), intent(out) :: corr(size(x, 2), size(x, 2)) !! Symmetric Kendall tau-b correlation matrix.
    integer :: i, j

    corr = 0.0_dp
    do i = 1, size(x, 2)
      corr(i, i) = 1.0_dp
      do j = i + 1, size(x, 2)
        corr(i, j) = cor_fk_vector(x(:, i), x(:, j))
        corr(j, i) = corr(i, j)
      end do
    end do
  end subroutine cor_fk_matrix

  pure real(dp) function spatial_objective(x, m) result(value)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for the spatial-median objective.
    real(dp), intent(in) :: m(:) !! Candidate spatial median, one coordinate per variable.
    integer :: i

    value = 0.0_dp
    do i = 1, size(x, 1)
      value = value + sqrt(sum((x(i, :) - m)**2))
    end do
  end function spatial_objective

  pure subroutine spatial_gradient(x, m, zero_tol, gradient)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for the spatial-median objective.
    real(dp), intent(in) :: m(:) !! Candidate spatial median, one coordinate per variable.
    real(dp), intent(in) :: zero_tol !! Distance threshold below which a coincident observation is skipped.
    real(dp), intent(out) :: gradient(size(m)) !! Gradient of the sum-of-Euclidean-distances objective.
    real(dp) :: r
    integer :: i

    gradient = 0.0_dp
    do i = 1, size(x, 1)
      r = sqrt(sum((m - x(i, :))**2))
      if (r > zero_tol) gradient = gradient + (m - x(i, :)) / r
    end do
  end subroutine spatial_gradient

  pure subroutine spatial_hessian(x, m, zero_tol, hessian)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for the spatial-median objective.
    real(dp), intent(in) :: m(:) !! Candidate spatial median, one coordinate per variable.
    real(dp), intent(in) :: zero_tol !! Distance threshold below which the singular Hessian term is skipped.
    real(dp), intent(out) :: hessian(size(m), size(m)) !! Hessian approximation for the spatial objective.
    real(dp), allocatable :: d(:)
    real(dp) :: r
    integer :: i, j

    allocate(d(size(m)))
    hessian = 0.0_dp
    do i = 1, size(x, 1)
      d = m - x(i, :)
      r = sqrt(sum(d**2))
      if (r <= zero_tol) cycle
      do j = 1, size(m)
        hessian(j, j) = hessian(j, j) + 1.0_dp / r
      end do
      hessian = hessian - spread(d, 2, size(m)) * spread(d, 1, size(m)) / r**3
    end do
  end subroutine spatial_hessian

  pure subroutine symmetric_eigen_jacobi(a, values, vectors, tol, max_sweeps)
    real(dp), intent(in) :: a(:,:) !! Real symmetric matrix to diagonalize.
    real(dp), intent(out) :: values(size(a, 1)) !! Eigenvalues, sorted in decreasing order on return.
    real(dp), intent(out) :: vectors(size(a, 1), size(a, 1)) !! Corresponding orthonormal eigenvectors by column.
    real(dp), intent(in), optional :: tol !! Optional convergence threshold for the largest off-diagonal entry.
    integer, intent(in), optional :: max_sweeps !! Optional maximum number of Jacobi rotations.
    real(dp), allocatable :: work(:,:)
    real(dp) :: app, aqq, apq, c, s, tau, t, threshold, max_off, vip, viq
    integer :: i, j, p, q, iter, max_iter, n

    n = size(a, 1)
    work = a
    vectors = 0.0_dp
    do i = 1, n
      vectors(i, i) = 1.0_dp
    end do
    threshold = sqrt(epsilon(1.0_dp))
    if (present(tol)) threshold = tol
    max_iter = max(50, 50 * n * n)
    if (present(max_sweeps)) max_iter = max_sweeps

    do iter = 1, max_iter
      max_off = 0.0_dp
      p = 1
      q = min(2, n)
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(work(i, j)) > max_off) then
            max_off = abs(work(i, j))
            p = i
            q = j
          end if
        end do
      end do
      if (n <= 1 .or. max_off <= threshold) exit
      app = work(p, p)
      aqq = work(q, q)
      apq = work(p, q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t * t)
      s = t * c

      do j = 1, n
        if (j == p .or. j == q) cycle
        vip = work(j, p)
        viq = work(j, q)
        work(j, p) = c * vip - s * viq
        work(p, j) = work(j, p)
        work(j, q) = s * vip + c * viq
        work(q, j) = work(j, q)
      end do
      work(p, p) = c * c * app - 2.0_dp * c * s * apq + s * s * aqq
      work(q, q) = s * s * app + 2.0_dp * c * s * apq + c * c * aqq
      work(p, q) = 0.0_dp
      work(q, p) = 0.0_dp

      do j = 1, n
        vip = vectors(j, p)
        viq = vectors(j, q)
        vectors(j, p) = c * vip - s * viq
        vectors(j, q) = s * vip + c * viq
      end do
    end do

    do i = 1, n
      values(i) = work(i, i)
    end do
    call sort_eigen_desc(values, vectors)
  end subroutine symmetric_eigen_jacobi

  pure subroutine sort_eigen_desc(values, vectors)
    real(dp), intent(inout) :: values(:) !! Eigenvalues reordered from largest to smallest.
    real(dp), intent(inout) :: vectors(:,:) !! Eigenvector columns reordered consistently with values.
    real(dp), allocatable :: temp(:)
    real(dp) :: tv
    integer :: i, j, best

    allocate(temp(size(vectors, 1)))
    do i = 1, size(values) - 1
      best = i
      do j = i + 1, size(values)
        if (values(j) > values(best)) best = j
      end do
      if (best /= i) then
        tv = values(i)
        values(i) = values(best)
        values(best) = tv
        temp = vectors(:, i)
        vectors(:, i) = vectors(:, best)
        vectors(:, best) = temp
      end if
    end do
  end subroutine sort_eigen_desc

  pure subroutine orient_columns(a)
    real(dp), intent(inout) :: a(:,:) !! Matrix whose columns are sign-normalized by their largest absolute element.
    integer :: i, j, idx

    do j = 1, size(a, 2)
      idx = 1
      do i = 2, size(a, 1)
        if (abs(a(i, j)) > abs(a(idx, j))) idx = i
      end do
      if (a(idx, j) < 0.0_dp) a(:, j) = -a(:, j)
    end do
  end subroutine orient_columns

  subroutine random_normal_vector(x)
    real(dp), intent(out) :: x(:) !! Independent standard-normal variates generated with Box-Muller transforms.
    real(dp) :: u1, u2, r, theta
    integer :: i

    i = 1
    do while (i <= size(x))
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      r = sqrt(-2.0_dp * log(u1))
      theta = 2.0_dp * acos(-1.0_dp) * u2
      x(i) = r * cos(theta)
      if (i + 1 <= size(x)) x(i + 1) = r * sin(theta)
      i = i + 2
    end do
  end subroutine random_normal_vector

  subroutine random_normal_matrix(x)
    real(dp), intent(out) :: x(:,:) !! Matrix filled column-major with independent standard-normal variates.
    real(dp), allocatable :: work(:)

    allocate(work(size(x)))
    call random_normal_vector(work)
    x = reshape(work, shape(x))
  end subroutine random_normal_matrix

end module pcapp_stats
