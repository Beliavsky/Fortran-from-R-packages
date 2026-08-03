! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_linalg
  use multiatsm_kinds, only : dp
  implicit none
  private

  public :: solve_linear, inverse_matrix, least_squares, covariance_rows
  public :: symmetric_eigen, singular_value_decomposition, cholesky_lower
  public :: logdet_spd, spectral_radius, matrix_power, pseudo_inverse
  public :: block_diag, outer_product, frobenius_norm, eye

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import dp
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    end subroutine dgesv

    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      import dp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
    end subroutine dgels

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *), work(*)
      real(dp), intent(out) :: w(*)
    end subroutine dsyev

    subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobu, jobvt
      integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *), work(*)
      real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *)
    end subroutine dgesvd

    subroutine dpotrf(uplo, n, a, lda, info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *)
    end subroutine dpotrf

    subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobvl, jobvr
      integer, intent(in) :: n, lda, ldvl, ldvr, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *), work(*)
      real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *)
    end subroutine dgeev
  end interface

contains

  pure function eye(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function eye

  subroutine solve_linear(a, b, x, info, ridge)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), allocatable, intent(out) :: x(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: ridge
    real(dp), allocatable :: ac(:, :), bc(:, :)
    integer, allocatable :: ipiv(:)
    integer :: n, nrhs, i, attempt
    real(dp) :: reg

    n = size(a, 1)
    nrhs = size(b, 2)
    if (size(a, 2) /= n .or. size(b, 1) /= n) then
      info = -1
      allocate(x(0, 0))
      return
    end if

    reg = 0.0_dp
    if (present(ridge)) reg = max(ridge, 0.0_dp)
    allocate(ac(n, n), bc(n, nrhs), ipiv(n), x(n, nrhs))

    do attempt = 0, 6
      ac = a
      if (reg > 0.0_dp) then
        do i = 1, n
          ac(i, i) = ac(i, i) + reg
        end do
      end if
      bc = b
      call dgesv(n, nrhs, ac, n, ipiv, bc, n, info)
      if (info == 0) then
        x = bc
        return
      end if
      if (reg <= tiny(1.0_dp)) then
        reg = 1.0e-12_dp * max(1.0_dp, maxval(abs(a)))
      else
        reg = 10.0_dp * reg
      end if
    end do
    x = 0.0_dp
  end subroutine solve_linear

  subroutine inverse_matrix(a, ainv, info, ridge)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: ainv(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: ridge
    real(dp), allocatable :: rhs(:, :)
    integer :: n

    n = size(a, 1)
    if (size(a, 2) /= n) then
      info = -1
      allocate(ainv(0, 0))
      return
    end if
    allocate(rhs(n, n))
    rhs = eye(n)
    if (present(ridge)) then
      call solve_linear(a, rhs, ainv, info, ridge)
    else
      call solve_linear(a, rhs, ainv, info)
    end if
  end subroutine inverse_matrix

  subroutine least_squares(a, b, x, info)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), allocatable, intent(out) :: x(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: ac(:, :), bc(:, :), work(:)
    real(dp) :: workq(1)
    integer :: m, n, nrhs, ldb, lwork

    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    if (size(b, 1) /= m) then
      info = -1
      allocate(x(0, 0))
      return
    end if
    ldb = max(m, n)
    allocate(ac(m, n), bc(ldb, nrhs))
    ac = a
    bc = 0.0_dp
    bc(1:m, :) = b
    call dgels('N', m, n, nrhs, ac, m, bc, ldb, workq, -1, info)
    if (info /= 0) then
      allocate(x(0, 0))
      return
    end if
    lwork = max(1, int(workq(1)))
    allocate(work(lwork))
    ac = a
    bc = 0.0_dp
    bc(1:m, :) = b
    call dgels('N', m, n, nrhs, ac, m, bc, ldb, work, lwork, info)
    allocate(x(n, nrhs))
    if (info == 0) then
      x = bc(1:n, :)
    else
      x = 0.0_dp
    end if
  end subroutine least_squares

  function covariance_rows(x, demean) result(cov)
    real(dp), intent(in) :: x(:, :)
    logical, intent(in), optional :: demean
    real(dp) :: cov(size(x, 1), size(x, 1))
    real(dp), allocatable :: xc(:, :), mu(:)
    integer :: t
    logical :: center

    t = size(x, 2)
    center = .true.
    if (present(demean)) center = demean
    allocate(xc(size(x, 1), t))
    xc = x
    if (center .and. t > 0) then
      allocate(mu(size(x, 1)))
      mu = sum(x, dim=2) / real(t, dp)
      xc = x - spread(mu, 2, t)
    end if
    if (t > 1) then
      cov = matmul(xc, transpose(xc)) / real(t - 1, dp)
    else
      cov = 0.0_dp
    end if
  end function covariance_rows

  subroutine symmetric_eigen(a, values, vectors, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: ac(:, :), work(:)
    real(dp) :: workq(1)
    integer :: n, lwork

    n = size(a, 1)
    if (size(a, 2) /= n) then
      info = -1
      allocate(values(0), vectors(0, 0))
      return
    end if
    allocate(ac(n, n), values(n), vectors(n, n))
    ac = 0.5_dp * (a + transpose(a))
    call dsyev('V', 'U', n, ac, n, values, workq, -1, info)
    if (info /= 0) return
    lwork = max(1, int(workq(1)))
    allocate(work(lwork))
    ac = 0.5_dp * (a + transpose(a))
    call dsyev('V', 'U', n, ac, n, values, work, lwork, info)
    vectors = ac
  end subroutine symmetric_eigen

  subroutine singular_value_decomposition(a, u, s, vt, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: u(:, :), s(:), vt(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: ac(:, :), work(:)
    real(dp) :: workq(1)
    integer :: m, n, k, lwork

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)
    allocate(ac(m, n), u(m, m), s(k), vt(n, n))
    ac = a
    call dgesvd('A', 'A', m, n, ac, m, s, u, m, vt, n, workq, -1, info)
    if (info /= 0) return
    lwork = max(1, int(workq(1)))
    allocate(work(lwork))
    ac = a
    call dgesvd('A', 'A', m, n, ac, m, s, u, m, vt, n, work, lwork, info)
  end subroutine singular_value_decomposition

  subroutine pseudo_inverse(a, ap, info, tolerance)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: ap(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: u(:, :), s(:), vt(:, :), sinv(:, :)
    real(dp) :: tol
    integer :: i, m, n, k

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)
    call singular_value_decomposition(a, u, s, vt, info)
    if (info /= 0) then
      allocate(ap(0, 0))
      return
    end if
    tol = epsilon(1.0_dp) * real(max(m, n), dp) * maxval(s)
    if (present(tolerance)) tol = tolerance
    allocate(sinv(n, m))
    sinv = 0.0_dp
    do i = 1, k
      if (s(i) > tol) sinv(i, i) = 1.0_dp / s(i)
    end do
    allocate(ap(n, m))
    ap = matmul(transpose(vt), matmul(sinv, transpose(u)))
  end subroutine pseudo_inverse

  subroutine cholesky_lower(a, l, info, jitter)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: l(:, :)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: jitter
    integer :: n, i, attempt
    real(dp) :: eps

    n = size(a, 1)
    if (size(a, 2) /= n) then
      info = -1
      allocate(l(0, 0))
      return
    end if
    allocate(l(n, n))
    eps = 0.0_dp
    if (present(jitter)) eps = max(jitter, 0.0_dp)
    do attempt = 0, 7
      l = 0.5_dp * (a + transpose(a))
      do i = 1, n
        l(i, i) = l(i, i) + eps
      end do
      call dpotrf('L', n, l, n, info)
      if (info == 0) then
        do i = 1, n
          if (i < n) l(i, i+1:n) = 0.0_dp
        end do
        return
      end if
      if (eps <= tiny(1.0_dp)) then
        eps = 1.0e-12_dp * max(1.0_dp, maxval(abs(a)))
      else
        eps = 10.0_dp * eps
      end if
    end do
    l = 0.0_dp
  end subroutine cholesky_lower

  function logdet_spd(a, info) result(value)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out) :: info
    real(dp) :: value
    real(dp), allocatable :: l(:, :)
    integer :: i

    call cholesky_lower(a, l, info)
    if (info /= 0) then
      value = huge(1.0_dp)
      return
    end if
    value = 0.0_dp
    do i = 1, size(l, 1)
      value = value + 2.0_dp * log(l(i, i))
    end do
  end function logdet_spd

  function spectral_radius(a, info) result(rho)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out) :: info
    real(dp) :: rho
    real(dp), allocatable :: ac(:, :), wr(:), wi(:), work(:), vl(:, :), vr(:, :)
    real(dp) :: workq(1)
    integer :: n, lwork

    n = size(a, 1)
    if (size(a, 2) /= n) then
      info = -1
      rho = huge(1.0_dp)
      return
    end if
    allocate(ac(n, n), wr(n), wi(n), vl(1, 1), vr(1, 1))
    ac = a
    call dgeev('N', 'N', n, ac, n, wr, wi, vl, 1, vr, 1, workq, -1, info)
    if (info /= 0) then
      rho = huge(1.0_dp)
      return
    end if
    lwork = max(1, int(workq(1)))
    allocate(work(lwork))
    ac = a
    call dgeev('N', 'N', n, ac, n, wr, wi, vl, 1, vr, 1, work, lwork, info)
    if (info == 0) then
      rho = maxval(sqrt(wr * wr + wi * wi))
    else
      rho = huge(1.0_dp)
    end if
  end function spectral_radius

  function matrix_power(a, p) result(out)
    real(dp), intent(in) :: a(:, :)
    integer, intent(in) :: p
    real(dp) :: out(size(a, 1), size(a, 2))
    real(dp) :: base(size(a, 1), size(a, 2))
    integer :: q

    out = eye(size(a, 1))
    if (p <= 0) return
    base = a
    q = p
    do while (q > 0)
      if (mod(q, 2) == 1) out = matmul(out, base)
      q = q / 2
      if (q > 0) base = matmul(base, base)
    end do
  end function matrix_power

  function block_diag(a, b) result(c)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp) :: c(size(a, 1) + size(b, 1), size(a, 2) + size(b, 2))
    c = 0.0_dp
    c(1:size(a, 1), 1:size(a, 2)) = a
    c(size(a, 1)+1:, size(a, 2)+1:) = b
  end function block_diag

  pure function outer_product(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))
    c = spread(a, 2, size(b)) * spread(b, 1, size(a))
  end function outer_product

  pure function frobenius_norm(a) result(value)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: value
    value = sqrt(sum(a * a))
  end function frobenius_norm

end module multiatsm_linalg
