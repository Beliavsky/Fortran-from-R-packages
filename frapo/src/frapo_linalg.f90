! SPDX-License-Identifier: GPL-3.0-or-later
module frapo_linalg
  use frapo_kinds, only : dp
  use frapo_types, only : frapo_ok, frapo_singular, frapo_invalid_input
  implicit none
  private

  public :: solve_linear, symmetric_eigen, symmetric_matrix_sqrt, sqrm

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import dp
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    end subroutine dgesv

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*)
      real(dp), intent(inout) :: work(*)
    end subroutine dsyev
  end interface

contains

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: acopy(:, :), rhs(:, :)
    integer, allocatable :: ipiv(:)
    integer :: n, info

    n = size(b)
    allocate(x(n))
    if (size(a, 1) /= n .or. size(a, 2) /= n) then
      x = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    allocate(acopy(n, n), rhs(n, 1), ipiv(n))
    acopy = a
    rhs(:, 1) = b
    call dgesv(n, 1, acopy, n, ipiv, rhs, n, info)
    if (info == 0) then
      x = rhs(:, 1)
      if (present(status)) status = frapo_ok
    else
      x = 0.0_dp
      if (present(status)) status = frapo_singular
    end if
  end subroutine solve_linear

  subroutine symmetric_eigen(a, values, vectors, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: work(:)
    real(dp) :: work_query(1)
    integer :: n, lwork, info

    n = size(a, 1)
    allocate(values(n), vectors(n, n))
    if (size(a, 2) /= n) then
      values = 0.0_dp
      vectors = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    vectors = a
    call dsyev('V', 'U', n, vectors, n, values, work_query, -1, info)
    lwork = max(1, int(work_query(1)))
    allocate(work(lwork))
    call dsyev('V', 'U', n, vectors, n, values, work, lwork, info)
    if (present(status)) then
      if (info == 0) then
        status = frapo_ok
      else
        status = frapo_singular
      end if
    end if
  end subroutine symmetric_eigen

  subroutine symmetric_matrix_sqrt(a, root, status, clip_negative)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: root(:, :)
    integer, intent(out), optional :: status
    logical, intent(in), optional :: clip_negative
    real(dp), allocatable :: values(:), vectors(:, :), scaled(:, :)
    logical :: clip
    integer :: i, n, istat

    clip = .false.
    if (present(clip_negative)) clip = clip_negative
    call symmetric_eigen(a, values, vectors, istat)
    n = size(values)
    allocate(root(n, n), scaled(n, n))
    if (istat /= frapo_ok) then
      root = 0.0_dp
      if (present(status)) status = istat
      return
    end if
    if (any(values < 0.0_dp) .and. .not. clip) then
      root = 0.0_dp
      if (present(status)) status = frapo_invalid_input
      return
    end if
    do i = 1, n
      values(i) = sqrt(max(values(i), 0.0_dp))
      scaled(:, i) = vectors(:, i) * values(i)
    end do
    root = matmul(scaled, transpose(vectors))
    if (present(status)) status = frapo_ok
  end subroutine symmetric_matrix_sqrt

  function sqrm(a, status) result(root)
    real(dp), intent(in) :: a(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: root(:, :)
    integer :: istat

    call symmetric_matrix_sqrt(a, root, istat)
    if (present(status)) status = istat
  end function sqrm
end module frapo_linalg
