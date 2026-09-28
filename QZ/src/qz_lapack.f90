module qz_lapack
    use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, ieee_value
    use qz_kinds, only: dp
    use qz_types, only: qz_eigen_real_result, qz_eigen_complex_result, &
        qz_geigen_real_result, qz_geigen_complex_result, &
        qz_schur_real_result, qz_schur_complex_result, &
        qz_gschur_real_result, qz_gschur_complex_result, &
        qz_reorder_real_result, qz_reorder_complex_result, &
        qz_greorder_real_result, qz_greorder_complex_result
    use la_lapack, only: gees, geev, gges, ggev, tgsen, trsen
    implicit none
    private

    public :: qz_dgees, qz_dgeev, qz_dgges, qz_dggev
    public :: qz_dtgsen, qz_dtrsen
    public :: qz_zgees, qz_zgeev, qz_zgges, qz_zggev
    public :: qz_ztgsen, qz_ztrsen

contains

    subroutine qz_dgeev(a, result, vl, vr, lwork_in)
        real(dp), intent(in) :: a(:,:) !! Real square matrix whose eigensystem is required.
        type(qz_eigen_real_result), intent(out) :: result !! Eigenvalues, optional eigenvectors, workspace hint, and LAPACK status.
        logical, intent(in), optional :: vl !! If true, compute left eigenvectors; default is true.
        logical, intent(in), optional :: vr !! If true, compute right eigenvectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvl, ldvr
        logical :: want_vl, want_vr
        character :: jobvl, jobvr
        real(dp), allocatable :: acopy(:,:), work(:), vl_raw(:,:), vr_raw(:,:)

        n = size(a, 1)
        if (size(a, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vl = .true.
        want_vr = .true.
        if (present(vl)) want_vl = vl
        if (present(vr)) want_vr = vr
        jobvl = merge('V', 'N', want_vl)
        jobvr = merge('V', 'N', want_vr)
        ldvl = merge(max(1, n), 1, want_vl)
        ldvr = merge(max(1, n), 1, want_vr)
        lwork = max(1, 4 * max(1, n))
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(acopy(n, n), result%wr(n), result%wi(n), work(lwork))
        allocate(vl_raw(ldvl, max(1, n)), vr_raw(ldvr, max(1, n)))
        acopy = a
        vl_raw = 0.0_dp
        vr_raw = 0.0_dp
        call geev(jobvl, jobvr, n, acopy, max(1, n), result%wr, result%wi, &
            vl_raw, ldvl, vr_raw, ldvr, work, lwork, result%info)
        result%work_opt = work(1)
        allocate(result%w(n))
        result%w = cmplx(result%wr, result%wi, kind=dp)
        if (want_vl) then
            allocate(result%vl(n, n))
            result%vl = vl_raw(1:n, 1:n)
            call real_eigenvectors_to_complex(result%vl, result%wi, result%u)
        end if
        if (want_vr) then
            allocate(result%vr(n, n))
            result%vr = vr_raw(1:n, 1:n)
            call real_eigenvectors_to_complex(result%vr, result%wi, result%v)
        end if
    end subroutine qz_dgeev

    subroutine qz_zgeev(a, result, vl, vr, lwork_in)
        complex(dp), intent(in) :: a(:,:) !! Complex square matrix whose eigensystem is required.
        type(qz_eigen_complex_result), intent(out) :: result !! Eigenvalues, optional vectors, workspace hint, and LAPACK status.
        logical, intent(in), optional :: vl !! If true, compute left eigenvectors; default is true.
        logical, intent(in), optional :: vr !! If true, compute right eigenvectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvl, ldvr
        logical :: want_vl, want_vr
        character :: jobvl, jobvr
        complex(dp), allocatable :: acopy(:,:), work(:), vl_raw(:,:), vr_raw(:,:)
        real(dp), allocatable :: rwork(:)

        n = size(a, 1)
        if (size(a, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vl = .true.
        want_vr = .true.
        if (present(vl)) want_vl = vl
        if (present(vr)) want_vr = vr
        jobvl = merge('V', 'N', want_vl)
        jobvr = merge('V', 'N', want_vr)
        ldvl = merge(max(1, n), 1, want_vl)
        ldvr = merge(max(1, n), 1, want_vr)
        lwork = max(1, 2 * max(1, n))
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(acopy(n, n), result%w(n), work(lwork), rwork(max(1, 2 * n)))
        allocate(vl_raw(ldvl, max(1, n)), vr_raw(ldvr, max(1, n)))
        acopy = a
        vl_raw = cmplx(0.0_dp, 0.0_dp, kind=dp)
        vr_raw = cmplx(0.0_dp, 0.0_dp, kind=dp)
        call geev(jobvl, jobvr, n, acopy, max(1, n), result%w, vl_raw, ldvl, &
            vr_raw, ldvr, work, lwork, rwork, result%info)
        result%work_opt = work(1)
        if (want_vl) then
            allocate(result%vl(n, n), result%u(n, n))
            result%vl = vl_raw(1:n, 1:n)
            result%u = result%vl
        end if
        if (want_vr) then
            allocate(result%vr(n, n), result%v(n, n))
            result%vr = vr_raw(1:n, 1:n)
            result%v = result%vr
        end if
    end subroutine qz_zgeev

    subroutine qz_dggev(a, b, result, vl, vr, lwork_in)
        real(dp), intent(in) :: a(:,:) !! First real matrix of the generalized eigenproblem A*v=lambda*B*v.
        real(dp), intent(in) :: b(:,:) !! Second real matrix of the generalized eigenproblem, with the same square shape as A.
        type(qz_geigen_real_result), intent(out) :: result !! Generalized eigenvalue factors, vectors, workspace hint, and status.
        logical, intent(in), optional :: vl !! If true, compute left generalized eigenvectors; default is true.
        logical, intent(in), optional :: vr !! If true, compute right generalized eigenvectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvl, ldvr
        logical :: want_vl, want_vr
        character :: jobvl, jobvr
        real(dp), allocatable :: acopy(:,:), bcopy(:,:), work(:), vl_raw(:,:), vr_raw(:,:)

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vl = .true.
        want_vr = .true.
        if (present(vl)) want_vl = vl
        if (present(vr)) want_vr = vr
        jobvl = merge('V', 'N', want_vl)
        jobvr = merge('V', 'N', want_vr)
        ldvl = merge(max(1, n), 1, want_vl)
        ldvr = merge(max(1, n), 1, want_vr)
        lwork = max(1, 8 * max(1, n))
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(acopy(n, n), bcopy(n, n), result%alphar(n), result%alphai(n), result%beta(n))
        allocate(work(lwork), vl_raw(ldvl, max(1, n)), vr_raw(ldvr, max(1, n)))
        acopy = a
        bcopy = b
        vl_raw = 0.0_dp
        vr_raw = 0.0_dp
        call ggev(jobvl, jobvr, n, acopy, max(1, n), bcopy, max(1, n), &
            result%alphar, result%alphai, result%beta, vl_raw, ldvl, vr_raw, ldvr, &
            work, lwork, result%info)
        result%work_opt = work(1)
        allocate(result%alpha(n), result%lambda(n))
        result%alpha = cmplx(result%alphar, result%alphai, kind=dp)
        call generalized_ratio(result%alpha, cmplx(result%beta, 0.0_dp, kind=dp), result%lambda)
        if (want_vl) then
            allocate(result%vl(n, n))
            result%vl = vl_raw(1:n, 1:n)
            call real_eigenvectors_to_complex(result%vl, result%alphai, result%u)
        end if
        if (want_vr) then
            allocate(result%vr(n, n))
            result%vr = vr_raw(1:n, 1:n)
            call real_eigenvectors_to_complex(result%vr, result%alphai, result%v)
        end if
    end subroutine qz_dggev

    subroutine qz_zggev(a, b, result, vl, vr, lwork_in)
        complex(dp), intent(in) :: a(:,:) !! First complex matrix of the generalized eigenproblem A*v=lambda*B*v.
        complex(dp), intent(in) :: b(:,:) !! Second complex matrix of the generalized eigenproblem, with the same square shape as A.
        type(qz_geigen_complex_result), intent(out) :: result !! Generalized eigenpairs, vectors, workspace hint, and status.
        logical, intent(in), optional :: vl !! If true, compute left generalized eigenvectors; default is true.
        logical, intent(in), optional :: vr !! If true, compute right generalized eigenvectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvl, ldvr
        logical :: want_vl, want_vr
        character :: jobvl, jobvr
        complex(dp), allocatable :: acopy(:,:), bcopy(:,:), work(:), vl_raw(:,:), vr_raw(:,:)
        real(dp), allocatable :: rwork(:)

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vl = .true.
        want_vr = .true.
        if (present(vl)) want_vl = vl
        if (present(vr)) want_vr = vr
        jobvl = merge('V', 'N', want_vl)
        jobvr = merge('V', 'N', want_vr)
        ldvl = merge(max(1, n), 1, want_vl)
        ldvr = merge(max(1, n), 1, want_vr)
        lwork = max(1, 2 * max(1, n))
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(acopy(n, n), bcopy(n, n), result%alpha(n), result%beta(n), result%lambda(n))
        allocate(work(lwork), rwork(max(1, 8 * n)))
        allocate(vl_raw(ldvl, max(1, n)), vr_raw(ldvr, max(1, n)))
        acopy = a
        bcopy = b
        vl_raw = cmplx(0.0_dp, 0.0_dp, kind=dp)
        vr_raw = cmplx(0.0_dp, 0.0_dp, kind=dp)
        call ggev(jobvl, jobvr, n, acopy, max(1, n), bcopy, max(1, n), result%alpha, &
            result%beta, vl_raw, ldvl, vr_raw, ldvr, work, lwork, rwork, result%info)
        result%work_opt = work(1)
        call generalized_ratio(result%alpha, result%beta, result%lambda)
        if (want_vl) then
            allocate(result%vl(n, n), result%u(n, n))
            result%vl = vl_raw(1:n, 1:n)
            result%u = result%vl
        end if
        if (want_vr) then
            allocate(result%vr(n, n), result%v(n, n))
            result%vr = vr_raw(1:n, 1:n)
            result%v = result%vr
        end if
    end subroutine qz_zggev

    subroutine qz_dgees(a, result, vs, lwork_in)
        real(dp), intent(in) :: a(:,:) !! Real square matrix to reduce to real Schur form.
        type(qz_schur_real_result), intent(out) :: result !! Schur form, eigenvalues, optional Schur vectors, and status.
        logical, intent(in), optional :: vs !! If true, compute Schur vectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvs
        logical :: want_vs
        character :: jobvs
        real(dp), allocatable :: work(:), vs_raw(:,:)
        logical, allocatable :: bwork(:)

        n = size(a, 1)
        if (size(a, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vs = .true.
        if (present(vs)) want_vs = vs
        jobvs = merge('V', 'N', want_vs)
        ldvs = merge(max(1, n), 1, want_vs)
        lwork = max(1, 3 * max(1, n))
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(result%t(n, n), result%wr(n), result%wi(n), result%w(n))
        allocate(work(lwork), vs_raw(ldvs, max(1, n)), bwork(max(1, n)))
        result%t = a
        vs_raw = 0.0_dp
        bwork = .false.
        call gees(jobvs, 'N', select_none_real, n, result%t, max(1, n), result%sdim, &
            result%wr, result%wi, vs_raw, ldvs, work, lwork, bwork, result%info)
        result%work_opt = work(1)
        result%w = cmplx(result%wr, result%wi, kind=dp)
        if (want_vs) then
            allocate(result%q(n, n))
            result%q = vs_raw(1:n, 1:n)
        end if
    end subroutine qz_dgees

    subroutine qz_zgees(a, result, vs, lwork_in)
        complex(dp), intent(in) :: a(:,:) !! Complex square matrix to reduce to complex Schur form.
        type(qz_schur_complex_result), intent(out) :: result !! Schur form, eigenvalues, optional Schur vectors, and status.
        logical, intent(in), optional :: vs !! If true, compute Schur vectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvs
        logical :: want_vs
        character :: jobvs
        complex(dp), allocatable :: work(:), vs_raw(:,:)
        real(dp), allocatable :: rwork(:)
        logical, allocatable :: bwork(:)

        n = size(a, 1)
        if (size(a, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vs = .true.
        if (present(vs)) want_vs = vs
        jobvs = merge('V', 'N', want_vs)
        ldvs = merge(max(1, n), 1, want_vs)
        lwork = max(1, 2 * max(1, n))
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(result%t(n, n), result%w(n), work(lwork), rwork(max(1, n)), bwork(max(1, n)))
        allocate(vs_raw(ldvs, max(1, n)))
        result%t = a
        vs_raw = cmplx(0.0_dp, 0.0_dp, kind=dp)
        bwork = .false.
        call gees(jobvs, 'N', select_none_complex, n, result%t, max(1, n), result%sdim, &
            result%w, vs_raw, ldvs, work, lwork, rwork, bwork, result%info)
        result%work_opt = work(1)
        if (want_vs) then
            allocate(result%q(n, n))
            result%q = vs_raw(1:n, 1:n)
        end if
    end subroutine qz_zgees

    subroutine qz_dgges(a, b, result, vsl, vsr, lwork_in)
        real(dp), intent(in) :: a(:,:) !! First real matrix in the generalized Schur decomposition.
        real(dp), intent(in) :: b(:,:) !! Second real matrix in the generalized Schur decomposition.
        type(qz_gschur_real_result), intent(out) :: result !! Generalized Schur factors, eigenvalue factors, vectors, and status.
        logical, intent(in), optional :: vsl !! If true, compute left generalized Schur vectors; default is true.
        logical, intent(in), optional :: vsr !! If true, compute right generalized Schur vectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvsl, ldvsr
        logical :: want_vsl, want_vsr
        character :: jobvsl, jobvsr
        real(dp), allocatable :: work(:), vsl_raw(:,:), vsr_raw(:,:)
        logical, allocatable :: bwork(:)

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vsl = .true.
        want_vsr = .true.
        if (present(vsl)) want_vsl = vsl
        if (present(vsr)) want_vsr = vsr
        jobvsl = merge('V', 'N', want_vsl)
        jobvsr = merge('V', 'N', want_vsr)
        ldvsl = merge(max(1, n), 1, want_vsl)
        ldvsr = merge(max(1, n), 1, want_vsr)
        lwork = max(1, 8 * max(1, n) + 16)
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(result%s(n, n), result%t(n, n), result%alphar(n), result%alphai(n), result%beta(n))
        allocate(result%alpha(n), result%lambda(n), work(lwork), bwork(max(1, n)))
        allocate(vsl_raw(ldvsl, max(1, n)), vsr_raw(ldvsr, max(1, n)))
        result%s = a
        result%t = b
        vsl_raw = 0.0_dp
        vsr_raw = 0.0_dp
        bwork = .false.
        call gges(jobvsl, jobvsr, 'N', select_none_generalized_real, n, result%s, max(1, n), &
            result%t, max(1, n), result%sdim, result%alphar, result%alphai, result%beta, &
            vsl_raw, ldvsl, vsr_raw, ldvsr, work, lwork, bwork, result%info)
        result%work_opt = work(1)
        result%alpha = cmplx(result%alphar, result%alphai, kind=dp)
        call generalized_ratio(result%alpha, cmplx(result%beta, 0.0_dp, kind=dp), result%lambda)
        if (want_vsl) then
            allocate(result%q(n, n))
            result%q = vsl_raw(1:n, 1:n)
        end if
        if (want_vsr) then
            allocate(result%z(n, n))
            result%z = vsr_raw(1:n, 1:n)
        end if
    end subroutine qz_dgges

    subroutine qz_zgges(a, b, result, vsl, vsr, lwork_in)
        complex(dp), intent(in) :: a(:,:) !! First complex matrix in the generalized Schur decomposition.
        complex(dp), intent(in) :: b(:,:) !! Second complex matrix in the generalized Schur decomposition.
        type(qz_gschur_complex_result), intent(out) :: result !! Generalized Schur factors, eigenvalue factors, vectors, and status.
        logical, intent(in), optional :: vsl !! If true, compute left generalized Schur vectors; default is true.
        logical, intent(in), optional :: vsr !! If true, compute right generalized Schur vectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested LAPACK work length; a safe minimum is enforced.
        integer :: n, lwork, ldvsl, ldvsr
        logical :: want_vsl, want_vsr
        character :: jobvsl, jobvsr
        complex(dp), allocatable :: work(:), vsl_raw(:,:), vsr_raw(:,:)
        real(dp), allocatable :: rwork(:)
        logical, allocatable :: bwork(:)

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(b, 2) /= n) then
            result%info = -1001
            return
        end if
        want_vsl = .true.
        want_vsr = .true.
        if (present(vsl)) want_vsl = vsl
        if (present(vsr)) want_vsr = vsr
        jobvsl = merge('V', 'N', want_vsl)
        jobvsr = merge('V', 'N', want_vsr)
        ldvsl = merge(max(1, n), 1, want_vsl)
        ldvsr = merge(max(1, n), 1, want_vsr)
        lwork = max(1, 2 * max(1, n))
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(result%s(n, n), result%t(n, n), result%alpha(n), result%beta(n), result%lambda(n))
        allocate(work(lwork), rwork(max(1, 8 * n)), bwork(max(1, n)))
        allocate(vsl_raw(ldvsl, max(1, n)), vsr_raw(ldvsr, max(1, n)))
        result%s = a
        result%t = b
        vsl_raw = cmplx(0.0_dp, 0.0_dp, kind=dp)
        vsr_raw = cmplx(0.0_dp, 0.0_dp, kind=dp)
        bwork = .false.
        call gges(jobvsl, jobvsr, 'N', select_none_generalized_complex, n, result%s, max(1, n), &
            result%t, max(1, n), result%sdim, result%alpha, result%beta, vsl_raw, ldvsl, &
            vsr_raw, ldvsr, work, lwork, rwork, bwork, result%info)
        result%work_opt = work(1)
        call generalized_ratio(result%alpha, result%beta, result%lambda)
        if (want_vsl) then
            allocate(result%q(n, n))
            result%q = vsl_raw(1:n, 1:n)
        end if
        if (want_vsr) then
            allocate(result%z(n, n))
            result%z = vsr_raw(1:n, 1:n)
        end if
    end subroutine qz_zgges

    subroutine qz_dtrsen(t, q, select, result, job, want_q, lwork_in, liwork_in)
        real(dp), intent(in) :: t(:,:) !! Real Schur-form matrix to reorder.
        real(dp), intent(in) :: q(:,:) !! Real Schur-vector matrix associated with T.
        logical, intent(in) :: select(:) !! Selection mask for eigenvalues or real 2-by-2 Schur blocks.
        type(qz_reorder_real_result), intent(out) :: result !! Reordered Schur form, eigenvalues, condition estimates, and status.
        character(len=1), intent(in), optional :: job !! LAPACK TRSEN condition-estimate job: B, V, E, or N; default B.
        logical, intent(in), optional :: want_q !! If true, update Schur vectors Q; default is true.
        integer, intent(in), optional :: lwork_in !! Requested real work-array length; a conservative safe minimum is enforced.
        integer, intent(in), optional :: liwork_in !! Requested integer work-array length; a conservative safe minimum is enforced.
        integer :: n, lwork, liwork
        logical :: update_q
        character :: job_use, compq
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)

        n = size(t, 1)
        if (size(t, 2) /= n .or. size(q, 1) /= n .or. size(q, 2) /= n .or. size(select) /= n) then
            result%info = -1001
            return
        end if
        job_use = 'B'
        if (present(job)) job_use = job
        update_q = .true.
        if (present(want_q)) update_q = want_q
        compq = merge('V', 'N', update_q)
        lwork = max(1, 4 * n * n + 32 * n + 64)
        liwork = max(1, 4 * n * n + 32 * n + 64)
        if (present(lwork_in)) lwork = max(lwork, lwork_in)
        if (present(liwork_in)) liwork = max(liwork, liwork_in)

        allocate(result%t(n, n), result%q(n, n), result%wr(n), result%wi(n), result%w(n))
        allocate(work(lwork), iwork(liwork))
        result%t = t
        result%q = q
        call trsen(job_use, compq, select, n, result%t, max(1, n), result%q, max(1, n), &
            result%wr, result%wi, result%m, result%s, result%sep, work, lwork, iwork, liwork, result%info)
        result%work_opt = work(1)
        result%iwork_opt = iwork(1)
        result%w = cmplx(result%wr, result%wi, kind=dp)
        if (.not. update_q) deallocate(result%q)
    end subroutine qz_dtrsen

    subroutine qz_ztrsen(t, q, select, result, job, want_q, lwork_in)
        complex(dp), intent(in) :: t(:,:) !! Complex Schur-form matrix to reorder.
        complex(dp), intent(in) :: q(:,:) !! Complex Schur-vector matrix associated with T.
        logical, intent(in) :: select(:) !! Selection mask for eigenvalues to move to the leading Schur block.
        type(qz_reorder_complex_result), intent(out) :: result !! Reordered Schur form, eigenvalues, estimates, and status.
        character(len=1), intent(in), optional :: job !! LAPACK TRSEN condition-estimate job: B, V, E, or N; default B.
        logical, intent(in), optional :: want_q !! If true, update Schur vectors Q; default is true.
        integer, intent(in), optional :: lwork_in !! Requested complex work-array length; a conservative safe minimum is enforced.
        integer :: n, lwork
        logical :: update_q
        character :: job_use, compq
        complex(dp), allocatable :: work(:)

        n = size(t, 1)
        if (size(t, 2) /= n .or. size(q, 1) /= n .or. size(q, 2) /= n .or. size(select) /= n) then
            result%info = -1001
            return
        end if
        job_use = 'B'
        if (present(job)) job_use = job
        update_q = .true.
        if (present(want_q)) update_q = want_q
        compq = merge('V', 'N', update_q)
        lwork = max(1, 4 * n * n + 32 * n + 64)
        if (present(lwork_in)) lwork = max(lwork, lwork_in)

        allocate(result%t(n, n), result%q(n, n), result%w(n), work(lwork))
        result%t = t
        result%q = q
        call trsen(job_use, compq, select, n, result%t, max(1, n), result%q, max(1, n), &
            result%w, result%m, result%s, result%sep, work, lwork, result%info)
        result%work_opt = work(1)
        if (.not. update_q) deallocate(result%q)
    end subroutine qz_ztrsen

    subroutine qz_dtgsen(s, t, q, z, select, result, ijob, want_q, want_z, lwork_in, liwork_in)
        real(dp), intent(in) :: s(:,:) !! First real generalized Schur-form matrix.
        real(dp), intent(in) :: t(:,:) !! Second real generalized Schur-form matrix.
        real(dp), intent(in) :: q(:,:) !! Left generalized Schur-vector matrix.
        real(dp), intent(in) :: z(:,:) !! Right generalized Schur-vector matrix.
        logical, intent(in) :: select(:) !! Selection mask for generalized eigenvalue blocks to move to the leading block.
        type(qz_greorder_real_result), intent(out) :: result !! Reordered generalized Schur form, condition estimates, and status.
        integer, intent(in), optional :: ijob !! LAPACK TGSEN condition-estimate selector from 0 through 5; default 4.
        logical, intent(in), optional :: want_q !! If true, update left Schur vectors; default is true.
        logical, intent(in), optional :: want_z !! If true, update right Schur vectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested real work-array length; a conservative safe minimum is enforced.
        integer, intent(in), optional :: liwork_in !! Requested integer work-array length; a conservative safe minimum is enforced.
        integer :: n, lwork, liwork, ijob_use
        logical :: update_q, update_z
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)

        n = size(s, 1)
        if (size(s, 2) /= n .or. size(t, 1) /= n .or. size(t, 2) /= n .or. &
            size(q, 1) /= n .or. size(q, 2) /= n .or. size(z, 1) /= n .or. &
            size(z, 2) /= n .or. size(select) /= n) then
            result%info = -1001
            return
        end if
        ijob_use = 4
        if (present(ijob)) ijob_use = ijob
        update_q = .true.
        update_z = .true.
        if (present(want_q)) update_q = want_q
        if (present(want_z)) update_z = want_z
        lwork = max(1, 4 * n * n + 32 * n + 64)
        liwork = max(1, 4 * n * n + 32 * n + 64)
        if (present(lwork_in)) lwork = max(lwork, lwork_in)
        if (present(liwork_in)) liwork = max(liwork, liwork_in)

        allocate(result%s(n, n), result%t(n, n), result%q(n, n), result%z(n, n))
        allocate(result%alphar(n), result%alphai(n), result%beta(n), result%alpha(n), result%lambda(n))
        allocate(work(lwork), iwork(liwork))
        result%s = s
        result%t = t
        result%q = q
        result%z = z
        call tgsen(ijob_use, update_q, update_z, select, n, result%s, max(1, n), &
            result%t, max(1, n), result%alphar, result%alphai, result%beta, &
            result%q, max(1, n), result%z, max(1, n), result%m, result%pl, result%pr, &
            result%dif, work, lwork, iwork, liwork, result%info)
        result%work_opt = work(1)
        result%iwork_opt = iwork(1)
        result%alpha = cmplx(result%alphar, result%alphai, kind=dp)
        call generalized_ratio(result%alpha, cmplx(result%beta, 0.0_dp, kind=dp), result%lambda)
        if (.not. update_q) deallocate(result%q)
        if (.not. update_z) deallocate(result%z)
    end subroutine qz_dtgsen

    subroutine qz_ztgsen(s, t, q, z, select, result, ijob, want_q, want_z, lwork_in, liwork_in)
        complex(dp), intent(in) :: s(:,:) !! First complex generalized Schur-form matrix.
        complex(dp), intent(in) :: t(:,:) !! Second complex generalized Schur-form matrix.
        complex(dp), intent(in) :: q(:,:) !! Left generalized Schur-vector matrix.
        complex(dp), intent(in) :: z(:,:) !! Right generalized Schur-vector matrix.
        logical, intent(in) :: select(:) !! Selection mask for generalized eigenvalues to move to the leading block.
        type(qz_greorder_complex_result), intent(out) :: result !! Reordered generalized Schur form, estimates, and status.
        integer, intent(in), optional :: ijob !! LAPACK TGSEN condition-estimate selector from 0 through 5; default 4.
        logical, intent(in), optional :: want_q !! If true, update left Schur vectors; default is true.
        logical, intent(in), optional :: want_z !! If true, update right Schur vectors; default is true.
        integer, intent(in), optional :: lwork_in !! Requested complex work-array length; a conservative safe minimum is enforced.
        integer, intent(in), optional :: liwork_in !! Requested integer work-array length; a conservative safe minimum is enforced.
        integer :: n, lwork, liwork, ijob_use
        logical :: update_q, update_z
        complex(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)

        n = size(s, 1)
        if (size(s, 2) /= n .or. size(t, 1) /= n .or. size(t, 2) /= n .or. &
            size(q, 1) /= n .or. size(q, 2) /= n .or. size(z, 1) /= n .or. &
            size(z, 2) /= n .or. size(select) /= n) then
            result%info = -1001
            return
        end if
        ijob_use = 4
        if (present(ijob)) ijob_use = ijob
        update_q = .true.
        update_z = .true.
        if (present(want_q)) update_q = want_q
        if (present(want_z)) update_z = want_z
        lwork = max(1, 4 * n * n + 32 * n + 64)
        liwork = max(1, 4 * n * n + 32 * n + 64)
        if (present(lwork_in)) lwork = max(lwork, lwork_in)
        if (present(liwork_in)) liwork = max(liwork, liwork_in)

        allocate(result%s(n, n), result%t(n, n), result%q(n, n), result%z(n, n))
        allocate(result%alpha(n), result%beta(n), result%lambda(n), work(lwork), iwork(liwork))
        result%s = s
        result%t = t
        result%q = q
        result%z = z
        call tgsen(ijob_use, update_q, update_z, select, n, result%s, max(1, n), &
            result%t, max(1, n), result%alpha, result%beta, result%q, max(1, n), &
            result%z, max(1, n), result%m, result%pl, result%pr, result%dif, &
            work, lwork, iwork, liwork, result%info)
        result%work_opt = work(1)
        result%iwork_opt = iwork(1)
        call generalized_ratio(result%alpha, result%beta, result%lambda)
        if (.not. update_q) deallocate(result%q)
        if (.not. update_z) deallocate(result%z)
    end subroutine qz_ztgsen

    subroutine real_eigenvectors_to_complex(raw, imag_part, vectors)
        real(dp), intent(in) :: raw(:,:) !! LAPACK real eigenvector storage, with imaginary parts in adjacent columns.
        real(dp), intent(in) :: imag_part(:) !! Eigenvalue imaginary parts identifying real and conjugate-pair columns.
        complex(dp), allocatable, intent(out) :: vectors(:,:) !! Reconstructed complex eigenvectors with one vector per column.
        integer :: n, j

        n = size(raw, 1)
        allocate(vectors(n, size(raw, 2)))
        vectors = cmplx(raw, 0.0_dp, kind=dp)
        j = 1
        do while (j <= size(raw, 2))
            if (imag_part(j) > 0.0_dp .and. j < size(raw, 2)) then
                vectors(:, j) = cmplx(raw(:, j), raw(:, j + 1), kind=dp)
                vectors(:, j + 1) = conjg(vectors(:, j))
                j = j + 2
            else
                j = j + 1
            end if
        end do
    end subroutine real_eigenvectors_to_complex

    subroutine generalized_ratio(alpha, beta, lambda)
        complex(dp), intent(in) :: alpha(:) !! Generalized-eigenvalue numerator factors.
        complex(dp), intent(in) :: beta(:) !! Generalized-eigenvalue denominator factors.
        complex(dp), intent(out) :: lambda(:) !! Ratios alpha/beta, using IEEE nonfinite values when beta is effectively zero.
        integer :: i
        real(dp) :: inf_value, nan_value, real_part, imag_part, scale

        inf_value = ieee_value(0.0_dp, ieee_positive_inf)
        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        do i = 1, size(alpha)
            scale = abs(beta(i))
            if (scale > tiny(1.0_dp)) then
                lambda(i) = alpha(i) / beta(i)
            else if (abs(alpha(i)) <= tiny(1.0_dp)) then
                lambda(i) = cmplx(nan_value, nan_value, kind=dp)
            else
                real_part = 0.0_dp
                imag_part = 0.0_dp
                if (abs(real(alpha(i), dp)) > tiny(1.0_dp)) then
                    real_part = sign(inf_value, real(alpha(i), dp))
                end if
                if (abs(aimag(alpha(i))) > tiny(1.0_dp)) then
                    imag_part = sign(inf_value, aimag(alpha(i)))
                end if
                lambda(i) = cmplx(real_part, imag_part, kind=dp)
            end if
        end do
    end subroutine generalized_ratio

    logical function select_none_real(wr, wi)
        real(dp), intent(in) :: wr !! Candidate real part supplied by GEES; unused because sorting is disabled.
        real(dp), intent(in) :: wi !! Candidate imaginary part supplied by GEES; unused because sorting is disabled.
        select_none_real = abs(wr) < 0.0_dp .and. abs(wi) < 0.0_dp
    end function select_none_real

    logical function select_none_complex(w)
        complex(dp), intent(in) :: w !! Candidate complex eigenvalue supplied by GEES; unused because sorting is disabled.
        select_none_complex = abs(w) < 0.0_dp
    end function select_none_complex

    logical function select_none_generalized_real(alphar, alphai, beta)
        real(dp), intent(in) :: alphar !! Candidate generalized-eigenvalue real numerator; unused because sorting is disabled.
        real(dp), intent(in) :: alphai !! Candidate generalized-eigenvalue imaginary numerator; unused because sorting is disabled.
        real(dp), intent(in) :: beta !! Candidate generalized-eigenvalue denominator; unused because sorting is disabled.
        select_none_generalized_real = abs(alphar) + abs(alphai) + abs(beta) < 0.0_dp
    end function select_none_generalized_real

    logical function select_none_generalized_complex(alpha, beta)
        complex(dp), intent(in) :: alpha !! Candidate generalized-eigenvalue numerator; unused because sorting is disabled.
        complex(dp), intent(in) :: beta !! Candidate generalized-eigenvalue denominator; unused because sorting is disabled.
        select_none_generalized_complex = abs(alpha) + abs(beta) < 0.0_dp
    end function select_none_generalized_complex
end module qz_lapack
