! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_linalg
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use intrinsicfrp_kinds, only: dp, status_ok, status_invalid, status_singular
  implicit none
  private
  public :: column_means, center_columns, covariance_matrix, cross_covariance
  public :: solve_linear, inverse_matrix, symmetric_eigen, symmetric_pinv
  public :: singular_values, thin_svd, cholesky_lower, matrix_rank
  public :: diag_vector, outer_product, identity_matrix, median_value
  public :: all_finite_matrix, all_finite_vector, correlation_matrix, kron_matrix

contains

  pure function column_means(x) result(mu)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: mu(size(x, 2))
    if (size(x, 1) > 0) then
      mu = sum(x, dim=1) / real(size(x, 1), dp)
    else
      mu = 0.0_dp
    end if
  end function column_means

  pure function center_columns(x) result(y)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: y(size(x, 1), size(x, 2))
    integer :: j
    real(dp) :: mu(size(x, 2))
    mu = column_means(x)
    do j = 1, size(x, 2)
      y(:, j) = x(:, j) - mu(j)
    end do
  end function center_columns

  pure function covariance_matrix(x) result(cov)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: cov(size(x, 2), size(x, 2))
    real(dp) :: xc(size(x, 1), size(x, 2))
    integer :: n
    n = size(x, 1)
    if (n <= 1) then
      cov = 0.0_dp
      return
    end if
    xc = center_columns(x)
    cov = matmul(transpose(xc), xc) / real(n - 1, dp)
    cov = 0.5_dp * (cov + transpose(cov))
  end function covariance_matrix

  pure function cross_covariance(x, y) result(cov)
    real(dp), intent(in) :: x(:, :), y(:, :)
    real(dp) :: cov(size(x, 2), size(y, 2))
    real(dp) :: xc(size(x, 1), size(x, 2))
    real(dp) :: yc(size(y, 1), size(y, 2))
    integer :: n
    n = size(x, 1)
    if (n <= 1 .or. size(y, 1) /= n) then
      cov = 0.0_dp
      return
    end if
    xc = center_columns(x)
    yc = center_columns(y)
    cov = matmul(transpose(xc), yc) / real(n - 1, dp)
  end function cross_covariance

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function identity_matrix

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: i
    do i = 1, size(x)
      a(i, :) = x(i) * y
    end do
  end function outer_product

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), allocatable, intent(out) :: x(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:, :), rowtmp(:)
    real(dp) :: pivot, factor, scale, tol
    integer :: n, m, i, k, p

    n = size(a, 1)
    m = size(b, 2)
    status = status_invalid
    allocate(x(max(0, n), max(0, m)))
    x = 0.0_dp
    if (size(a, 2) /= n .or. size(b, 1) /= n .or. n == 0) return

    allocate(aug(n, n + m), rowtmp(n + m))
    aug(:, 1:n) = a
    aug(:, n + 1:n + m) = b
    scale = max(1.0_dp, maxval(abs(a)))
    tol = 100.0_dp * epsilon(1.0_dp) * scale * real(n, dp)

    do k = 1, n
      p = k - 1 + maxloc(abs(aug(k:n, k)), dim=1)
      pivot = abs(aug(p, k))
      if (pivot <= tol) then
        call solve_pseudoinverse(a, b, x, status)
        return
      end if
      if (p /= k) then
        rowtmp = aug(k, :)
        aug(k, :) = aug(p, :)
        aug(p, :) = rowtmp
      end if
      aug(k, :) = aug(k, :) / aug(k, k)
      do i = 1, n
        if (i == k) cycle
        factor = aug(i, k)
        if (abs(factor) > 0.0_dp) aug(i, :) = aug(i, :) - factor * aug(k, :)
      end do
    end do
    x = aug(:, n + 1:n + m)
    status = status_ok
  end subroutine solve_linear

  subroutine solve_pseudoinverse(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), allocatable, intent(inout) :: x(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: pinv(:, :)
    integer :: st
    call symmetric_pinv(0.5_dp * (a + transpose(a)), pinv, st)
    if (st /= status_ok) then
      status = status_singular
      x = 0.0_dp
      return
    end if
    x = matmul(pinv, b)
    status = status_ok
  end subroutine solve_pseudoinverse

  subroutine inverse_matrix(a, inva, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: inva(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: eye(:, :)
    integer :: n
    n = size(a, 1)
    allocate(eye(n, n))
    eye = identity_matrix(n)
    call solve_linear(a, eye, inva, status)
  end subroutine inverse_matrix

  subroutine symmetric_eigen(a, values, vectors, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: work(:, :)
    real(dp) :: app, aqq, apq, tau, t, c, s, aik, akq, vip, viq
    real(dp) :: off, tol
    integer :: n, p, q, i, sweep, max_sweeps, j

    n = size(a, 1)
    allocate(values(n), vectors(n, n), work(n, n))
    values = 0.0_dp
    vectors = identity_matrix(n)
    status = status_invalid
    if (size(a, 2) /= n .or. n == 0) return
    work = 0.5_dp * (a + transpose(a))
    tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(work)))
    max_sweeps = max(50, 20 * n * n)

    do sweep = 1, max_sweeps
      off = 0.0_dp
      p = 1
      q = min(2, n)
      do i = 1, n - 1
        j = i + maxloc(abs(work(i, i + 1:n)), dim=1)
        if (abs(work(i, j)) > off) then
          off = abs(work(i, j))
          p = i
          q = j
        end if
      end do
      if (off <= tol .or. n == 1) exit
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

      do i = 1, n
        if (i /= p .and. i /= q) then
          aik = work(i, p)
          akq = work(i, q)
          work(i, p) = c * aik - s * akq
          work(p, i) = work(i, p)
          work(i, q) = s * aik + c * akq
          work(q, i) = work(i, q)
        end if
      end do
      work(p, p) = c * c * app - 2.0_dp * s * c * apq + s * s * aqq
      work(q, q) = s * s * app + 2.0_dp * s * c * apq + c * c * aqq
      work(p, q) = 0.0_dp
      work(q, p) = 0.0_dp
      do i = 1, n
        vip = vectors(i, p)
        viq = vectors(i, q)
        vectors(i, p) = c * vip - s * viq
        vectors(i, q) = s * vip + c * viq
      end do
    end do

    do i = 1, n
      values(i) = work(i, i)
    end do
    call sort_eigenpairs(values, vectors)
    status = status_ok
  end subroutine symmetric_eigen

  subroutine sort_eigenpairs(values, vectors)
    real(dp), intent(inout) :: values(:), vectors(:, :)
    real(dp) :: tv
    real(dp), allocatable :: col(:)
    integer :: i, j, k, n
    n = size(values)
    allocate(col(size(vectors, 1)))
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if (values(j) < values(k)) k = j
      end do
      if (k /= i) then
        tv = values(i)
        values(i) = values(k)
        values(k) = tv
        col = vectors(:, i)
        vectors(:, i) = vectors(:, k)
        vectors(:, k) = col
      end if
    end do
  end subroutine sort_eigenpairs

  subroutine symmetric_pinv(a, pinv, status, rcond)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: pinv(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: rcond
    real(dp), allocatable :: eval(:), evec(:, :)
    real(dp) :: tol, rc
    integer :: i, n, st
    n = size(a, 1)
    allocate(pinv(n, n))
    pinv = 0.0_dp
    call symmetric_eigen(a, eval, evec, st)
    if (st /= status_ok) then
      status = st
      return
    end if
    rc = sqrt(epsilon(1.0_dp))
    if (present(rcond)) rc = max(0.0_dp, rcond)
    tol = rc * max(1.0_dp, maxval(abs(eval)))
    do i = 1, n
      if (abs(eval(i)) > tol) then
        pinv = pinv + outer_product(evec(:, i), evec(:, i)) / eval(i)
      end if
    end do
    status = status_ok
  end subroutine symmetric_pinv

  subroutine singular_values(a, values, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:)
    integer, intent(out) :: status
    real(dp), allocatable :: eval(:), evec(:, :), gram(:, :)
    integer :: nsv, i, st
    if (size(a, 1) >= size(a, 2)) then
      gram = matmul(transpose(a), a)
    else
      gram = matmul(a, transpose(a))
    end if
    call symmetric_eigen(gram, eval, evec, st)
    nsv = size(eval)
    allocate(values(nsv))
    do i = 1, nsv
      values(i) = sqrt(max(0.0_dp, eval(nsv - i + 1)))
    end do
    status = st
  end subroutine singular_values

  subroutine thin_svd(a, u, s, v, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: u(:, :), s(:), v(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: eval(:), evec(:, :), gram(:, :)
    real(dp) :: tol, normc
    integer :: m, n, r, i, idx, st
    m = size(a, 1)
    n = size(a, 2)
    r = min(m, n)
    allocate(u(m, r), s(r), v(n, r))
    u = 0.0_dp
    v = 0.0_dp
    s = 0.0_dp
    if (m >= n) then
      gram = matmul(transpose(a), a)
      call symmetric_eigen(gram, eval, evec, st)
      do i = 1, r
        idx = n - i + 1
        s(i) = sqrt(max(0.0_dp, eval(idx)))
        v(:, i) = evec(:, idx)
      end do
      tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(s))
      do i = 1, r
        if (s(i) > tol) u(:, i) = matmul(a, v(:, i)) / s(i)
      end do
    else
      gram = matmul(a, transpose(a))
      call symmetric_eigen(gram, eval, evec, st)
      do i = 1, r
        idx = m - i + 1
        s(i) = sqrt(max(0.0_dp, eval(idx)))
        u(:, i) = evec(:, idx)
      end do
      tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(s))
      do i = 1, r
        if (s(i) > tol) v(:, i) = matmul(transpose(a), u(:, i)) / s(i)
        normc = sqrt(max(0.0_dp, dot_product(v(:, i), v(:, i))))
        if (normc > 0.0_dp) v(:, i) = v(:, i) / normc
      end do
    end if
    status = st
  end subroutine thin_svd

  subroutine cholesky_lower(a, l, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: l(:, :)
    integer, intent(out) :: status
    real(dp) :: temp, tol
    integer :: n, i, j, k
    n = size(a, 1)
    allocate(l(n, n))
    l = 0.0_dp
    status = status_invalid
    if (size(a, 2) /= n) return
    tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))
    do i = 1, n
      do j = 1, i
        temp = a(i, j)
        do k = 1, j - 1
          temp = temp - l(i, k) * l(j, k)
        end do
        if (i == j) then
          if (temp <= tol) then
            status = status_singular
            return
          end if
          l(i, j) = sqrt(temp)
        else
          l(i, j) = temp / l(j, j)
        end if
      end do
    end do
    status = status_ok
  end subroutine cholesky_lower

  integer function matrix_rank(a, tol_in) result(r)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tol_in
    real(dp), allocatable :: sv(:)
    real(dp) :: tol
    integer :: st
    call singular_values(a, sv, st)
    if (size(sv) == 0) then
      r = 0
      return
    end if
    tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(sv))
    if (present(tol_in)) tol = tol_in
    r = count(sv > tol)
  end function matrix_rank

  pure function diag_vector(a) result(d)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: d(min(size(a, 1), size(a, 2)))
    integer :: i
    do i = 1, size(d)
      d(i) = a(i, i)
    end do
  end function diag_vector

  function median_value(x) result(med)
    real(dp), intent(in) :: x(:)
    real(dp) :: med
    real(dp), allocatable :: y(:)
    real(dp) :: t
    integer :: i, j, n
    n = size(x)
    if (n == 0) then
      med = 0.0_dp
      return
    end if
    y = x
    do i = 2, n
      t = y(i)
      j = i - 1
      do while (j >= 1)
        if (y(j) <= t) exit
        y(j + 1) = y(j)
        j = j - 1
      end do
      y(j + 1) = t
    end do
    if (mod(n, 2) == 1) then
      med = y((n + 1) / 2)
    else
      med = 0.5_dp * (y(n / 2) + y(n / 2 + 1))
    end if
  end function median_value

  pure logical function all_finite_matrix(x)
    real(dp), intent(in) :: x(:, :)
    all_finite_matrix = all(ieee_is_finite(x))
  end function all_finite_matrix

  pure logical function all_finite_vector(x)
    real(dp), intent(in) :: x(:)
    all_finite_vector = all(ieee_is_finite(x))
  end function all_finite_vector

  function correlation_matrix(x, y) result(cor)
    real(dp), intent(in) :: x(:, :), y(:, :)
    real(dp) :: cor(size(x, 2), size(y, 2))
    real(dp) :: cov(size(x, 2), size(y, 2))
    real(dp) :: vx(size(x, 2)), vy(size(y, 2)), den
    integer :: i, j
    cov = cross_covariance(x, y)
    vx = diag_vector(covariance_matrix(x))
    vy = diag_vector(covariance_matrix(y))
    do i = 1, size(x, 2)
      do j = 1, size(y, 2)
        den = sqrt(max(0.0_dp, vx(i) * vy(j)))
        if (den > 0.0_dp) then
          cor(i, j) = cov(i, j) / den
        else
          cor(i, j) = 0.0_dp
        end if
      end do
    end do
  end function correlation_matrix

  function kron_matrix(a, b) result(k)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp) :: k(size(a, 1) * size(b, 1), size(a, 2) * size(b, 2))
    integer :: i, j, r1, r2, c1, c2
    r1 = size(b, 1)
    c1 = size(b, 2)
    do i = 1, size(a, 1)
      do j = 1, size(a, 2)
        r2 = (i - 1) * r1
        c2 = (j - 1) * c1
        k(r2 + 1:r2 + r1, c2 + 1:c2 + c1) = a(i, j) * b
      end do
    end do
  end function kron_matrix

end module intrinsicfrp_linalg
