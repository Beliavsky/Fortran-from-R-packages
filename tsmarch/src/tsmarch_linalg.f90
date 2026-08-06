! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_linalg
  use ghyp_kinds, only : dp, i8
  use ghyp_linalg, only : cholesky_lower, inverse_spd, logdet_spd, symmetrize
  implicit none
  private
  public :: column_mean, sample_covariance, covariance_to_correlation
  public :: sample_correlation, nearest_correlation, is_positive_definite
  public :: jacobi_eigen, symmetric_inverse_sqrt, symmetric_sqrt
  public :: outer_product, standardize_columns, matrix_power_symmetric
  public :: set_random_seed, random_normal_matrix, type7_quantile
  public :: matrix_inverse_general, determinant_general, rank_columns

  real(dp), parameter :: pi = acos(-1.0_dp)

contains

  function column_mean(x) result(mu)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable :: mu(:)
    integer :: j
    allocate(mu(size(x, 2)))
    do j = 1, size(x, 2)
      mu(j) = sum(x(:, j)) / real(size(x, 1), dp)
    end do
  end function column_mean

  function sample_covariance(x) result(cov)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable :: cov(:, :)
    real(dp), allocatable :: mu(:), centered(:, :)
    integer :: n
    n = size(x, 1)
    allocate(cov(size(x, 2), size(x, 2)))
    if (n <= 1) then
      cov = 0.0_dp
      return
    end if
    mu = column_mean(x)
    allocate(centered(n, size(x, 2)))
    centered = x - spread(mu, 1, n)
    cov = matmul(transpose(centered), centered) / real(n - 1, dp)
    call symmetrize(cov)
  end function sample_covariance

  function covariance_to_correlation(cov) result(cor)
    real(dp), intent(in) :: cov(:, :)
    real(dp), allocatable :: cor(:, :)
    real(dp), allocatable :: sd(:)
    integer :: i, j, n
    n = size(cov, 1)
    allocate(cor(n, n), sd(n))
    do i = 1, n
      sd(i) = sqrt(max(cov(i, i), tiny(1.0_dp)))
    end do
    do j = 1, n
      do i = 1, n
        cor(i, j) = cov(i, j) / (sd(i) * sd(j))
      end do
    end do
    do i = 1, n
      cor(i, i) = 1.0_dp
    end do
    call symmetrize(cor)
  end function covariance_to_correlation

  function sample_correlation(x, method) result(cor)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in), optional :: method
    real(dp), allocatable :: cor(:, :)
    real(dp), allocatable :: work(:, :)
    integer :: mth, j
    mth = 1
    if (present(method)) mth = method
    if (mth == 1) then
      cor = covariance_to_correlation(sample_covariance(x))
    else
      allocate(work(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
        work(:, j) = rank_columns(x(:, j))
      end do
      cor = covariance_to_correlation(sample_covariance(work))
      if (mth == 3) then
        cor = 2.0_dp * sin(pi * cor / 6.0_dp)
        do j = 1, size(cor, 1)
          cor(j, j) = 1.0_dp
        end do
      end if
    end if
  end function sample_correlation

  function rank_columns(x) result(ranks)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: ranks(:)
    integer, allocatable :: idx(:)
    integer :: i, j, k, n
    real(dp) :: avg
    n = size(x)
    allocate(ranks(n), idx(n))
    idx = [(i, i = 1, n)]
    do i = 2, n
      k = idx(i)
      j = i - 1
      do while (j >= 1)
        if (x(idx(j)) <= x(k)) exit
        idx(j + 1) = idx(j)
        j = j - 1
      end do
      idx(j + 1) = k
    end do
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (abs(x(idx(j + 1)) - x(idx(i))) > epsilon(1.0_dp) * &
          max(1.0_dp, abs(x(idx(j + 1))), abs(x(idx(i))))) exit
        j = j + 1
      end do
      avg = 0.5_dp * real(i + j, dp)
      ranks(idx(i:j)) = avg
      i = j + 1
    end do
  end function rank_columns

  function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), allocatable :: a(:, :)
    allocate(a(size(x), size(y)))
    a = spread(x, 2, size(y)) * spread(y, 1, size(x))
  end function outer_product

  subroutine jacobi_eigen(a, values, vectors, ok, tolerance, max_iterations)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: b(:, :)
    real(dp) :: tol, app, aqq, apq, tau, t, c, s, bip, biq, vip, viq
    integer :: n, i, p, q, iter, maxit, k, order_i
    real(dp), allocatable :: tmpv(:)
    n = size(a, 1)
    allocate(values(n), vectors(n, n), b(n, n), tmpv(n))
    if (size(a, 2) /= n) then
      ok = .false.
      values = 0.0_dp
      vectors = 0.0_dp
      return
    end if
    b = a
    call symmetrize(b)
    vectors = 0.0_dp
    do i = 1, n
      vectors(i, i) = 1.0_dp
    end do
    tol = 1.0e-12_dp
    if (present(tolerance)) tol = tolerance
    maxit = max(100, 50 * n * n)
    if (present(max_iterations)) maxit = max_iterations
    ok = .false.
    do iter = 1, maxit
      p = 1
      q = min(2, n)
      apq = 0.0_dp
      do i = 1, n - 1
        do k = i + 1, n
          if (abs(b(i, k)) > abs(apq)) then
            apq = b(i, k)
            p = i
            q = k
          end if
        end do
      end do
      if (abs(apq) <= tol * max(1.0_dp, maxval(abs(b)))) then
        ok = .true.
        exit
      end if
      app = b(p, p)
      aqq = b(q, q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t * t)
      s = t * c
      do i = 1, n
        if (i /= p .and. i /= q) then
          bip = b(i, p)
          biq = b(i, q)
          b(i, p) = c * bip - s * biq
          b(p, i) = b(i, p)
          b(i, q) = s * bip + c * biq
          b(q, i) = b(i, q)
        end if
        vip = vectors(i, p)
        viq = vectors(i, q)
        vectors(i, p) = c * vip - s * viq
        vectors(i, q) = s * vip + c * viq
      end do
      b(p, p) = c * c * app - 2.0_dp * s * c * apq + s * s * aqq
      b(q, q) = s * s * app + 2.0_dp * s * c * apq + c * c * aqq
      b(p, q) = 0.0_dp
      b(q, p) = 0.0_dp
    end do
    values = [(b(i, i), i = 1, n)]
    do i = 1, n - 1
      order_i = i
      do k = i + 1, n
        if (values(k) > values(order_i)) order_i = k
      end do
      if (order_i /= i) then
        t = values(i)
        values(i) = values(order_i)
        values(order_i) = t
        tmpv = vectors(:, i)
        vectors(:, i) = vectors(:, order_i)
        vectors(:, order_i) = tmpv
      end if
    end do
  end subroutine jacobi_eigen

  function symmetric_inverse_sqrt(a, ok, floor_value) result(out)
    real(dp), intent(in) :: a(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: floor_value
    real(dp), allocatable :: out(:, :), values(:), vectors(:, :), d(:, :)
    real(dp) :: floorv
    integer :: i, n
    call jacobi_eigen(a, values, vectors, ok)
    n = size(a, 1)
    allocate(out(n, n), d(n, n))
    out = 0.0_dp
    d = 0.0_dp
    if (.not. ok) return
    floorv = 1.0e-10_dp
    if (present(floor_value)) floorv = floor_value
    do i = 1, n
      d(i, i) = 1.0_dp / sqrt(max(values(i), floorv))
    end do
    out = matmul(vectors, matmul(d, transpose(vectors)))
    call symmetrize(out)
  end function symmetric_inverse_sqrt

  function symmetric_sqrt(a, ok, floor_value) result(out)
    real(dp), intent(in) :: a(:, :)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: floor_value
    real(dp), allocatable :: out(:, :), values(:), vectors(:, :), d(:, :)
    real(dp) :: floorv
    integer :: i, n
    call jacobi_eigen(a, values, vectors, ok)
    n = size(a, 1)
    allocate(out(n, n), d(n, n))
    out = 0.0_dp
    d = 0.0_dp
    if (.not. ok) return
    floorv = 0.0_dp
    if (present(floor_value)) floorv = floor_value
    do i = 1, n
      d(i, i) = sqrt(max(values(i), floorv))
    end do
    out = matmul(vectors, matmul(d, transpose(vectors)))
    call symmetrize(out)
  end function symmetric_sqrt

  function matrix_power_symmetric(a, power, ok) result(out)
    real(dp), intent(in) :: a(:, :), power
    logical, intent(out) :: ok
    real(dp), allocatable :: out(:, :), values(:), vectors(:, :), d(:, :)
    integer :: i, n
    call jacobi_eigen(a, values, vectors, ok)
    n = size(a, 1)
    allocate(out(n, n), d(n, n))
    out = 0.0_dp
    d = 0.0_dp
    if (.not. ok) return
    do i = 1, n
      if (values(i) <= 0.0_dp .and. power < 0.0_dp) then
        ok = .false.
        return
      end if
      d(i, i) = max(values(i), 0.0_dp) ** power
    end do
    out = matmul(vectors, matmul(d, transpose(vectors)))
    call symmetrize(out)
  end function matrix_power_symmetric

  logical function is_positive_definite(a)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable :: l(:, :)
    logical :: ok
    call cholesky_lower(a, l, ok)
    is_positive_definite = ok
  end function is_positive_definite

  function nearest_correlation(a, eigen_floor) result(cor)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: eigen_floor
    real(dp), allocatable :: cor(:, :), values(:), vectors(:, :), d(:, :)
    real(dp) :: floorv
    logical :: ok
    integer :: i, j, n
    n = size(a, 1)
    floorv = 1.0e-8_dp
    if (present(eigen_floor)) floorv = eigen_floor
    call jacobi_eigen(0.5_dp * (a + transpose(a)), values, vectors, ok)
    allocate(cor(n, n), d(n, n))
    if (.not. ok) then
      cor = 0.0_dp
      do i = 1, n
        cor(i, i) = 1.0_dp
      end do
      return
    end if
    d = 0.0_dp
    do i = 1, n
      d(i, i) = max(values(i), floorv)
    end do
    cor = matmul(vectors, matmul(d, transpose(vectors)))
    do j = 1, n
      do i = 1, n
        cor(i, j) = cor(i, j) / sqrt(max(cor(i, i), floorv) * max(cor(j, j), floorv))
      end do
    end do
    do i = 1, n
      cor(i, i) = 1.0_dp
    end do
    call symmetrize(cor)
  end function nearest_correlation

  subroutine standardize_columns(x, z, mu, sigma)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: z(:, :), mu(:), sigma(:)
    integer :: j, n
    n = size(x, 1)
    allocate(z(n, size(x, 2)), mu(size(x, 2)), sigma(size(x, 2)))
    mu = column_mean(x)
    do j = 1, size(x, 2)
      sigma(j) = sqrt(sum((x(:, j) - mu(j)) ** 2) / real(max(1, n - 1), dp))
      sigma(j) = max(sigma(j), sqrt(tiny(1.0_dp)))
      z(:, j) = (x(:, j) - mu(j)) / sigma(j)
    end do
  end subroutine standardize_columns

  subroutine set_random_seed(seed)
    integer(i8), intent(in) :: seed
    integer, allocatable :: put(:)
    integer :: n, i
    integer(i8) :: state
    call random_seed(size=n)
    allocate(put(n))
    state = max(1_i8, abs(seed))
    do i = 1, n
      state = modulo(6364136223846793005_i8 * state + 1442695040888963407_i8, huge(1_i8))
      put(i) = int(modulo(state, int(huge(1), i8) - 1_i8) + 1_i8)
    end do
    call random_seed(put=put)
  end subroutine set_random_seed

  function random_normal_matrix(n, m, seed) result(z)
    integer, intent(in) :: n, m
    integer(i8), intent(in), optional :: seed
    real(dp), allocatable :: z(:, :)
    real(dp) :: u1, u2
    integer :: i, j
    if (present(seed)) call set_random_seed(seed)
    allocate(z(n, m))
    do j = 1, m
      i = 1
      do while (i <= n)
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        z(i, j) = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
        if (i + 1 <= n) z(i + 1, j) = sqrt(-2.0_dp * log(u1)) * sin(2.0_dp * pi * u2)
        i = i + 2
      end do
    end do
  end function random_normal_matrix

  function type7_quantile(x, probability) result(q)
    real(dp), intent(in) :: x(:), probability
    real(dp) :: q, h, frac, key
    real(dp), allocatable :: work(:)
    integer :: i, j, k, n
    n = size(x)
    if (n == 0) then
      q = huge(1.0_dp)
      return
    end if
    allocate(work(n))
    work = x
    do i = 2, n
      key = work(i)
      j = i - 1
      do while (j >= 1)
        if (work(j) <= key) exit
        work(j + 1) = work(j)
        j = j - 1
      end do
      work(j + 1) = key
    end do
    if (probability <= 0.0_dp) then
      q = work(1)
    else if (probability >= 1.0_dp) then
      q = work(n)
    else
      h = 1.0_dp + real(n - 1, dp) * probability
      k = int(floor(h))
      frac = h - real(k, dp)
      if (k >= n) then
        q = work(n)
      else
        q = (1.0_dp - frac) * work(k) + frac * work(k + 1)
      end if
    end if
  end function type7_quantile

  subroutine matrix_inverse_general(a, ainv, ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: ainv(:, :)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:, :), tmp(:)
    real(dp) :: pivot, factor
    integer :: n, i, k, p
    n = size(a, 1)
    allocate(ainv(n, n), aug(n, 2 * n), tmp(2 * n))
    if (size(a, 2) /= n) then
      ok = .false.
      ainv = 0.0_dp
      return
    end if
    aug = 0.0_dp
    aug(:, 1:n) = a
    do i = 1, n
      aug(i, n + i) = 1.0_dp
    end do
    ok = .true.
    do k = 1, n
      p = k
      do i = k + 1, n
        if (abs(aug(i, k)) > abs(aug(p, k))) p = i
      end do
      if (abs(aug(p, k)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
        ok = .false.
        ainv = 0.0_dp
        return
      end if
      if (p /= k) then
        tmp = aug(k, :)
        aug(k, :) = aug(p, :)
        aug(p, :) = tmp
      end if
      pivot = aug(k, k)
      aug(k, :) = aug(k, :) / pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i, k)
        aug(i, :) = aug(i, :) - factor * aug(k, :)
      end do
    end do
    ainv = aug(:, n + 1:2 * n)
  end subroutine matrix_inverse_general

  function determinant_general(a, ok) result(det)
    real(dp), intent(in) :: a(:, :)
    logical, intent(out) :: ok
    real(dp) :: det, factor, temp
    real(dp), allocatable :: b(:, :)
    integer :: n, i, k, p, sign_det
    n = size(a, 1)
    if (size(a, 2) /= n) then
      ok = .false.
      det = 0.0_dp
      return
    end if
    allocate(b(n, n))
    b = a
    sign_det = 1
    det = 1.0_dp
    ok = .true.
    do k = 1, n
      p = k
      do i = k + 1, n
        if (abs(b(i, k)) > abs(b(p, k))) p = i
      end do
      if (abs(b(p, k)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
        ok = .false.
        det = 0.0_dp
        return
      end if
      if (p /= k) then
        do i = k, n
          temp = b(k, i)
          b(k, i) = b(p, i)
          b(p, i) = temp
        end do
        sign_det = -sign_det
      end if
      det = det * b(k, k)
      do i = k + 1, n
        factor = b(i, k) / b(k, k)
        b(i, k:n) = b(i, k:n) - factor * b(k, k:n)
      end do
    end do
    det = real(sign_det, dp) * det
  end function determinant_general

end module tsmarch_linalg
