! Computational translation of the R package irlba 2.3.7.
! Upstream package and native core: GPL-3 / GPL-3-or-later.
! See ../LICENSE and ../UPSTREAM.md for provenance and copyright details.
module irlba_linalg
  use irlba_kinds, only : dp
  use irlba_lapack, only : dgesdd, dgeqrf, dorgqr, dsyev, zgesdd
  implicit none
  private
  public :: vec_norm2, orthogonalize, thin_qr, svd_real, svd_complex, symmetric_eigen
  public :: fill_normal, sort_real_ascending

contains

  pure real(dp) function vec_norm2(x) result(r)
    real(dp), intent(in) :: x(:)
    r = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vec_norm2

  subroutine orthogonalize(x, y, ncol)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(inout) :: y(:)
    integer, intent(in) :: ncol
    real(dp), allocatable :: t(:)

    if (ncol <= 0) return
    allocate(t(ncol))
    t = matmul(transpose(x(:, 1:ncol)), y)
    y = y - matmul(x(:, 1:ncol), t)
  end subroutine orthogonalize

  subroutine fill_normal(x)
    real(dp), intent(out) :: x(:)
    real(dp) :: u1, u2, r, theta
    integer :: i
    real(dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp

    i = 1
    do while (i <= size(x))
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      r = sqrt(-2.0_dp * log(u1))
      theta = twopi * u2
      x(i) = r * cos(theta)
      if (i + 1 <= size(x)) x(i + 1) = r * sin(theta)
      i = i + 2
    end do
  end subroutine fill_normal

  subroutine thin_qr(a, q, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: q(:, :)
    integer, intent(out), optional :: info
    real(dp), allocatable :: work(:), tau(:), tmp(:, :)
    real(dp) :: workq(1)
    integer :: m, n, k, lwork, ierr

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)
    allocate(tmp(m, n), tau(k))
    tmp = a
    lwork = -1
    call dgeqrf(m, n, tmp, m, tau, workq, lwork, ierr)
    if (ierr /= 0) then
      if (present(info)) info = ierr
      allocate(q(m, k)); q = 0.0_dp
      return
    end if
    lwork = max(1, int(workq(1)))
    allocate(work(lwork))
    call dgeqrf(m, n, tmp, m, tau, work, lwork, ierr)
    if (ierr /= 0) then
      if (present(info)) info = ierr
      allocate(q(m, k)); q = 0.0_dp
      return
    end if
    deallocate(work)
    lwork = -1
    call dorgqr(m, k, k, tmp, m, tau, workq, lwork, ierr)
    lwork = max(1, int(workq(1)))
    allocate(work(lwork))
    call dorgqr(m, k, k, tmp, m, tau, work, lwork, ierr)
    allocate(q(m, k))
    if (ierr == 0) then
      q = tmp(:, 1:k)
    else
      q = 0.0_dp
    end if
    if (present(info)) info = ierr
  end subroutine thin_qr

  subroutine svd_real(a, s, u, v, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: s(:), u(:, :), v(:, :)
    integer, intent(out), optional :: info
    real(dp), allocatable :: ac(:, :), vt(:, :), work(:)
    real(dp) :: workq(1)
    integer, allocatable :: iwork(:)
    integer :: m, n, r, lwork, ierr

    m = size(a, 1)
    n = size(a, 2)
    r = min(m, n)
    allocate(ac(m, n), s(r), u(m, r), vt(r, n), iwork(8 * max(1, r)))
    ac = a
    lwork = -1
    call dgesdd('S', m, n, ac, m, s, u, m, vt, r, workq, lwork, iwork, ierr)
    if (ierr == 0) then
      lwork = max(1, int(workq(1)))
      allocate(work(lwork))
      ac = a
      call dgesdd('S', m, n, ac, m, s, u, m, vt, r, work, lwork, iwork, ierr)
    end if
    allocate(v(n, r))
    if (ierr == 0) then
      v = transpose(vt)
    else
      s = 0.0_dp
      u = 0.0_dp
      v = 0.0_dp
    end if
    if (present(info)) info = ierr
  end subroutine svd_real

  subroutine svd_complex(a, s, u, v, info)
    complex(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: s(:)
    complex(dp), allocatable, intent(out) :: u(:, :), v(:, :)
    integer, intent(out), optional :: info
    complex(dp), allocatable :: ac(:, :), vt(:, :), work(:)
    complex(dp) :: workq(1)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    integer :: m, n, r, lwork, ierr, rwsize

    m = size(a, 1)
    n = size(a, 2)
    r = min(m, n)
    rwsize = max(1, 5 * r * r + 7 * r)
    allocate(ac(m, n), s(r), u(m, r), vt(r, n), rwork(rwsize), iwork(8 * max(1, r)))
    ac = a
    lwork = -1
    call zgesdd('S', m, n, ac, m, s, u, m, vt, r, workq, lwork, rwork, iwork, ierr)
    if (ierr == 0) then
      lwork = max(1, int(real(workq(1), dp)))
      allocate(work(lwork))
      ac = a
      call zgesdd('S', m, n, ac, m, s, u, m, vt, r, work, lwork, rwork, iwork, ierr)
    end if
    allocate(v(n, r))
    if (ierr == 0) then
      v = transpose(conjg(vt))
    else
      s = 0.0_dp
      u = (0.0_dp, 0.0_dp)
      v = (0.0_dp, 0.0_dp)
    end if
    if (present(info)) info = ierr
  end subroutine svd_complex

  subroutine symmetric_eigen(a, values, vectors, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out), optional :: info
    real(dp), allocatable :: ac(:, :), work(:)
    real(dp) :: workq(1)
    integer :: n, lwork, ierr

    if (size(a, 1) /= size(a, 2)) error stop "symmetric_eigen: matrix must be square"
    n = size(a, 1)
    allocate(ac(n, n), values(n))
    ac = a
    lwork = -1
    call dsyev('V', 'U', n, ac, n, values, workq, lwork, ierr)
    if (ierr == 0) then
      lwork = max(1, int(workq(1)))
      allocate(work(lwork))
      ac = a
      call dsyev('V', 'U', n, ac, n, values, work, lwork, ierr)
    end if
    allocate(vectors(n, n))
    if (ierr == 0) then
      vectors = ac
      values = values(n:1:-1)
      vectors = vectors(:, n:1:-1)
    else
      values = 0.0_dp
      vectors = 0.0_dp
    end if
    if (present(info)) info = ierr
  end subroutine symmetric_eigen

  subroutine sort_real_ascending(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: t
    do i = 2, size(x)
      t = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= t) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = t
    end do
  end subroutine sort_real_ascending

end module irlba_linalg
