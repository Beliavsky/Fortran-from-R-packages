module irlba_lapack
  use irlba_kinds, only : dp
  implicit none
  private
  public :: dgesdd, dgeqrf, dorgqr, dsyev, zgesdd

  interface
    subroutine dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
      import :: dp
      character(len=1), intent(in) :: jobz
      integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: s(*)
      real(dp), intent(out) :: u(ldu, *)
      real(dp), intent(out) :: vt(ldvt, *)
      real(dp), intent(inout) :: work(*)
      integer, intent(out) :: iwork(*)
      integer, intent(out) :: info
    end subroutine dgesdd

    subroutine zgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
      import :: dp
      character(len=1), intent(in) :: jobz
      integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
      complex(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: s(*)
      complex(dp), intent(out) :: u(ldu, *)
      complex(dp), intent(out) :: vt(ldvt, *)
      complex(dp), intent(inout) :: work(*)
      real(dp), intent(inout) :: rwork(*)
      integer, intent(out) :: iwork(*)
      integer, intent(out) :: info
    end subroutine zgesdd

    subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
      import :: dp
      integer, intent(in) :: m, n, lda, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: tau(*)
      real(dp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dgeqrf

    subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
      import :: dp
      integer, intent(in) :: m, n, k, lda, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(in) :: tau(*)
      real(dp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dorgqr

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
end module irlba_lapack
