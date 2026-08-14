module mixsqp_lapack
  use mixsqp_kinds, only : dp
  implicit none
  private
  public :: dpotrf, dpotrs, dgesdd
  interface
    subroutine dpotrf(uplo,n,a,lda,info)
      import :: dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine dpotrf
    subroutine dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
      import :: dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, nrhs, lda, ldb
      real(dp), intent(in) :: a(lda,*)
      real(dp), intent(inout) :: b(ldb,*)
      integer, intent(out) :: info
    end subroutine dpotrs
    subroutine dgesdd(jobz,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,iwork,info)
      import :: dp
      character(len=1), intent(in) :: jobz
      integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: s(*)
      real(dp), intent(out) :: u(ldu,*), vt(ldvt,*), work(*)
      integer, intent(out) :: iwork(*), info
    end subroutine dgesdd
  end interface
end module mixsqp_lapack
