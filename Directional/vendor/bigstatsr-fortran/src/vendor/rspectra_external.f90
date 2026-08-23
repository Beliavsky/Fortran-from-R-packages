! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra_external
    use rspectra_kinds, only: dp
    implicit none
    private
    public :: dsaupd, dseupd, dnaupd, dneupd
    public :: dgetrf, dgetrs, zgetrf, zgetrs, dsyev, dgeev, dgesvd

    interface
        subroutine dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                          iparam, ipntr, workd, workl, lworkl, info)
            import dp
            integer :: ido, n, nev, ncv, ldv, lworkl, info
            character(len=1) :: bmat
            character(len=2) :: which
            real(dp) :: tol
            real(dp) :: resid(*), v(ldv,*), workd(*), workl(*)
            integer :: iparam(*), ipntr(*)
        end subroutine dsaupd

        subroutine dseupd(rvec, howmny, select, d, z, ldz, sigma, bmat, n, &
                          which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, &
                          workd, workl, lworkl, info)
            import dp
            logical :: rvec, select(*)
            character(len=1) :: howmny, bmat
            character(len=2) :: which
            integer :: ldz, n, nev, ncv, ldv, lworkl, info
            real(dp) :: d(*), z(ldz,*), sigma, tol, resid(*), v(ldv,*), workd(*), workl(*)
            integer :: iparam(*), ipntr(*)
        end subroutine dseupd

        subroutine dnaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                          iparam, ipntr, workd, workl, lworkl, info)
            import dp
            integer :: ido, n, nev, ncv, ldv, lworkl, info
            character(len=1) :: bmat
            character(len=2) :: which
            real(dp) :: tol
            real(dp) :: resid(*), v(ldv,*), workd(*), workl(*)
            integer :: iparam(*), ipntr(*)
        end subroutine dnaupd

        subroutine dneupd(rvec, howmny, select, dr, di, z, ldz, sigmar, sigmai, &
                          workev, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                          iparam, ipntr, workd, workl, lworkl, info)
            import dp
            logical :: rvec, select(*)
            character(len=1) :: howmny, bmat
            character(len=2) :: which
            integer :: ldz, n, nev, ncv, ldv, lworkl, info
            real(dp) :: dr(*), di(*), z(ldz,*), sigmar, sigmai, workev(*)
            real(dp) :: tol, resid(*), v(ldv,*), workd(*), workl(*)
            integer :: iparam(*), ipntr(*)
        end subroutine dneupd

        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import dp
            integer :: m, n, lda, ipiv(*), info
            real(dp) :: a(lda,*)
        end subroutine dgetrf

        subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            character(len=1) :: trans
            integer :: n, nrhs, lda, ipiv(*), ldb, info
            real(dp) :: a(lda,*), b(ldb,*)
        end subroutine dgetrs

        subroutine zgetrf(m, n, a, lda, ipiv, info)
            import dp
            integer :: m, n, lda, ipiv(*), info
            complex(dp) :: a(lda,*)
        end subroutine zgetrf

        subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            character(len=1) :: trans
            integer :: n, nrhs, lda, ipiv(*), ldb, info
            complex(dp) :: a(lda,*), b(ldb,*)
        end subroutine zgetrs

        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import dp
            character(len=1) :: jobz, uplo
            integer :: n, lda, lwork, info
            real(dp) :: a(lda,*), w(*), work(*)
        end subroutine dsyev

        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
                         work, lwork, info)
            import dp
            character(len=1) :: jobvl, jobvr
            integer :: n, lda, ldvl, ldvr, lwork, info
            real(dp) :: a(lda,*), wr(*), wi(*), vl(ldvl,*), vr(ldvr,*), work(*)
        end subroutine dgeev

        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import dp
            character(len=1) :: jobu, jobvt
            integer :: m, n, lda, ldu, ldvt, lwork, info
            real(dp) :: a(lda,*), s(*), u(ldu,*), vt(ldvt,*), work(*)
        end subroutine dgesvd
    end interface
end module rspectra_external
