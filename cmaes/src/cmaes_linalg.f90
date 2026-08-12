! cmaes-fortran - GPL-2.0-only
module cmaes_linalg
  use cmaes_kinds, only : dp
  implicit none
  private
  public :: symmetric_eigen_descending

  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*)
      real(dp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface
contains
  subroutine symmetric_eigen_descending(a, values, vectors, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: acopy(:, :), wasc(:), work(:)
    integer :: n, lwork, j

    n = size(a, 1)
    if (size(a, 2) /= n) then
      info = -1
      allocate(values(0), vectors(0, 0))
      return
    end if

    allocate(acopy(n, n), wasc(n))
    acopy = 0.5_dp * (a + transpose(a))
    lwork = max(1, 3 * n - 1)
    allocate(work(lwork))
    call dsyev('V', 'U', n, acopy, n, wasc, work, lwork, info)
    if (info /= 0) then
      allocate(values(0), vectors(0, 0))
      return
    end if

    allocate(values(n), vectors(n, n))
    do j = 1, n
      values(j) = wasc(n - j + 1)
      vectors(:, j) = acopy(:, n - j + 1)
    end do
  end subroutine symmetric_eigen_descending
end module cmaes_linalg
