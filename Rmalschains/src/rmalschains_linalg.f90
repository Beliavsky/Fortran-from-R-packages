module rmalschains_linalg
  use rmalschains_kinds, only : dp
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
    real(dp), allocatable :: work(:), tmpv(:)
    real(dp) :: query(1), tv
    integer :: n, lwork, i, j
    n = size(a, 1)
    if (size(a, 2) /= n) error stop "symmetric_eigen_descending: matrix must be square"
    allocate(values(n), vectors(n, n))
    vectors = 0.5_dp * (a + transpose(a))
    call dsyev('V', 'U', n, vectors, n, values, query, -1, info)
    if (info /= 0) return
    lwork = max(1, int(query(1)))
    allocate(work(lwork))
    call dsyev('V', 'U', n, vectors, n, values, work, lwork, info)
    if (info /= 0) return
    allocate(tmpv(n))
    do i = 1, n / 2
      j = n - i + 1
      tv = values(i); values(i) = values(j); values(j) = tv
      tmpv = vectors(:, i); vectors(:, i) = vectors(:, j); vectors(:, j) = tmpv
    end do
  end subroutine symmetric_eigen_descending
end module rmalschains_linalg
