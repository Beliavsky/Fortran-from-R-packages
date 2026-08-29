! Compatibility entry points for the legacy BLAS and LAPACK calls retained
! from KFAS. The implementation is supplied by the fortran-lapack package.

real(dp) function ddot(n, x, incx, y, incy)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_ddot
    implicit none
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: x(*), y(*)

    ddot = la_ddot(n, x, incx, y, incy)
end function ddot

subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dgemm
    implicit none
    character, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
    real(dp), intent(inout) :: c(ldc, *)

    call la_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine dgemm

subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dgemv
    implicit none
    character, intent(in) :: trans
    integer, intent(in) :: m, n, lda, incx, incy
    real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
    real(dp), intent(inout) :: y(*)

    call la_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine dgemv

subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dger
    implicit none
    integer, intent(in) :: m, n, incx, incy, lda
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda, *)

    call la_dger(m, n, alpha, x, incx, y, incy, a, lda)
end subroutine dger

subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dsymm
    implicit none
    character, intent(in) :: side, uplo
    integer, intent(in) :: m, n, lda, ldb, ldc
    real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
    real(dp), intent(inout) :: c(ldc, *)

    call la_dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine dsymm

subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dsymv
    implicit none
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda, incx, incy
    real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
    real(dp), intent(inout) :: y(*)

    call la_dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine dsymv

subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dsyr
    implicit none
    character, intent(in) :: uplo
    integer, intent(in) :: n, incx, lda
    real(dp), intent(in) :: alpha, x(*)
    real(dp), intent(inout) :: a(lda, *)

    call la_dsyr(uplo, n, alpha, x, incx, a, lda)
end subroutine dsyr

subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dsyr2
    implicit none
    character, intent(in) :: uplo
    integer, intent(in) :: n, incx, incy, lda
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda, *)

    call la_dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
end subroutine dsyr2

subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dsyrk
    implicit none
    character, intent(in) :: uplo, trans
    integer, intent(in) :: n, k, lda, ldc
    real(dp), intent(in) :: alpha, beta, a(lda, *)
    real(dp), intent(inout) :: c(ldc, *)

    call la_dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
end subroutine dsyrk

subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dtrmm
    implicit none
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    real(dp), intent(in) :: alpha, a(lda, *)
    real(dp), intent(inout) :: b(ldb, *)

    call la_dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
end subroutine dtrmm

subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
    use kfas_kinds, only: dp
    use la_blas_d, only: la_dtrmv
    implicit none
    character, intent(in) :: uplo, trans, diag
    integer, intent(in) :: n, lda, incx
    real(dp), intent(in) :: a(lda, *)
    real(dp), intent(inout) :: x(*)

    call la_dtrmv(uplo, trans, diag, n, a, lda, x, incx)
end subroutine dtrmv

subroutine dpotrf(uplo, n, a, lda, info)
    use kfas_kinds, only: dp
    use la_lapack_d, only: la_dpotrf
    implicit none
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda
    integer, intent(out) :: info
    real(dp), intent(inout) :: a(lda, *)

    call la_dpotrf(uplo, n, a, lda, info)
end subroutine dpotrf

subroutine dtrtri(uplo, diag, n, a, lda, info)
    use kfas_kinds, only: dp
    use la_lapack_d, only: la_dtrtri
    implicit none
    character, intent(in) :: uplo, diag
    integer, intent(in) :: n, lda
    integer, intent(out) :: info
    real(dp), intent(inout) :: a(lda, *)

    call la_dtrtri(uplo, diag, n, a, lda, info)
end subroutine dtrtri
