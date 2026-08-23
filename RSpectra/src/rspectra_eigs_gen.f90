! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra_eigs_gen_mod
    use rspectra_kinds, only: dp
    use rspectra_types, only: eigs_opts, eigs_result
    use rspectra_external, only: dnaupd, dneupd, dgetrf, dgetrs, zgetrf, zgetrs, dgeev
    use rspectra_operators, only: linear_operator, dense_operator, csr_operator, &
        make_dense_operator, csr_to_dense
    use rspectra_sort, only: sort_complex_pairs
    implicit none
    private
    public :: eigs_operator, eigs_dense, eigs_csr

contains

    function eigs_operator(op, k, which, opts) result(res)
        class(linear_operator), intent(inout) :: op
        integer, intent(in) :: k
        character(len=*), intent(in), optional :: which
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_result) :: res
        type(eigs_opts) :: o
        character(len=2) :: w
        integer :: n, ncv, ido, info, ierr, lworkl, nconv, ix, iy
        integer :: iparam(11), ipntr(14)
        real(dp), allocatable :: resid(:), v(:,:), workd(:), workl(:), dr(:), di(:), z(:,:), workev(:)
        logical, allocatable :: select(:)

        n = op%nrow
        if (op%nrow /= op%ncol .or. n < 3 .or. k < 1 .or. k >= n - 1) then
            res%info = -100
            return
        end if
        o = eigs_opts()
        if (present(opts)) o = opts
        w = 'LM'
        if (present(which)) w = adjustl(which)
        if (.not. valid_gen_which(w)) then
            res%info = -101
            return
        end if
        ncv = o%ncv
        if (ncv == 0) ncv = min(n, max(2 * k + 1, 20))
        if (ncv < k + 2 .or. ncv > n) then
            res%info = -102
            return
        end if
        lworkl = 3 * ncv * ncv + 6 * ncv
        allocate(resid(n), v(n,ncv), workd(3*n), workl(lworkl))
        allocate(dr(k+1), di(k+1), z(n,k+1), workev(3*ncv), select(ncv))
        resid = 0.0_dp
        v = 0.0_dp
        workd = 0.0_dp
        workl = 0.0_dp
        iparam = 0
        ipntr = 0
        iparam(1) = 1
        iparam(3) = o%maxitr
        iparam(7) = 1
        ido = 0
        if (allocated(o%initvec)) then
            if (size(o%initvec) /= n) then
                res%info = -103
                return
            end if
            resid = o%initvec
            info = 1
        else
            info = 0
        end if
        do
            call dnaupd(ido, 'I', n, w, k, o%tol, resid, ncv, v, n, iparam, ipntr, &
                workd, workl, lworkl, info)
            if (ido == -1 .or. ido == 1) then
                ix = ipntr(1)
                iy = ipntr(2)
                call op%prod(workd(ix:ix+n-1), workd(iy:iy+n-1))
            else if (ido == 2) then
                ix = ipntr(1)
                iy = ipntr(2)
                workd(iy:iy+n-1) = workd(ix:ix+n-1)
            else
                exit
            end if
        end do
        res%niter = iparam(3)
        res%nops = iparam(9)
        nconv = max(0, min(k, iparam(5)))
        res%nconv = nconv
        if (info < 0) then
            res%info = info
            return
        end if
        if (nconv == 0) then
            res%info = info
            allocate(res%values(0), res%vectors(n,0))
            return
        end if
        select = .false.
        ierr = 0
        call dneupd(o%retvec, 'A', select, dr, di, z, n, 0.0_dp, 0.0_dp, workev, &
            'I', n, w, k, o%tol, resid, ncv, v, n, iparam, ipntr, workd, workl, lworkl, ierr)
        if (ierr /= 0) then
            res%info = ierr
            return
        end if
        call pack_general_result(dr, di, z, nconv, o%retvec, res)
        if (o%retvec) then
            call sort_complex_pairs(res%values, res%vectors, w)
        else
            call sort_complex_pairs(res%values, which=w)
        end if
        res%info = info
    end function eigs_operator

    function eigs_dense(a, k, which, sigma, opts) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        character(len=*), intent(in), optional :: which
        complex(dp), intent(in), optional :: sigma
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_result) :: res
        type(dense_operator) :: op
        character(len=2) :: w
        integer :: n
        n = size(a,1)
        if (size(a,2) /= n .or. n < 1 .or. k < 1 .or. k > n) then
            res%info = -100
            return
        end if
        w = 'LM'
        if (present(which)) w = adjustl(which)
        if (.not. valid_gen_which(w)) then
            res%info = -101
            return
        end if
        if (k >= n - 1) then
            res = eigs_full(a, k, w, opts)
            return
        end if
        if (present(sigma)) then
            res = eigs_shift_full(a, k, w, sigma, opts)
        else
            op = make_dense_operator(a)
            res = eigs_operator(op, k, w, opts)
        end if
    end function eigs_dense

    function eigs_csr(op, k, which, sigma, opts) result(res)
        type(csr_operator), intent(inout) :: op
        integer, intent(in) :: k
        character(len=*), intent(in), optional :: which
        complex(dp), intent(in), optional :: sigma
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_result) :: res
        real(dp), allocatable :: a(:,:)
        character(len=2) :: w
        w = 'LM'
        if (present(which)) w = adjustl(which)
        if (present(sigma)) then
            a = csr_to_dense(op)
            res = eigs_dense(a, k, w, sigma, opts)
        else
            res = eigs_operator(op, k, w, opts)
        end if
    end function eigs_csr

    function eigs_shift_full(a, k, w, sigma, opts) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        character(len=2), intent(in) :: w
        complex(dp), intent(in) :: sigma
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_result) :: res
        type(eigs_result) :: full
        type(eigs_opts) :: o
        integer :: n
        n = size(a,1)
        o = eigs_opts()
        if (present(opts)) o = opts
        full = eigs_full(a, n, w, o)
        if (full%info /= 0) then
            res = full
            return
        end if
        if (o%retvec) then
            call sort_complex_pairs(full%values, full%vectors, w, sigma, .true.)
        else
            call sort_complex_pairs(full%values, which=w, sigma=sigma, shift_mode=.true.)
        end if
        allocate(res%values(k))
        res%values = full%values(1:k)
        if (o%retvec) then
            allocate(res%vectors(n,k))
            res%vectors = full%vectors(:,1:k)
        else
            allocate(res%vectors(n,0))
        end if
        res%nconv = k
        res%niter = 0
        res%nops = 0
        res%info = 0
    end function eigs_shift_full


    function eigs_full(a, k, w, opts) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        character(len=2), intent(in) :: w
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_result) :: res
        type(eigs_opts) :: o
        real(dp), allocatable :: ac(:,:), wr(:), wi(:), vl(:,:), vr(:,:), work(:)
        complex(dp), allocatable :: vals(:), vecs(:,:)
        integer :: n, lwork, info, nr
        n = size(a,1)
        o = eigs_opts()
        if (present(opts)) o = opts
        allocate(ac(n,n), wr(n), wi(n), vl(1,1), vr(n,n), work(1))
        ac = a
        lwork = -1
        call dgeev('N', merge('V','N',o%retvec), n, ac, n, wr, wi, vl, 1, vr, n, work, lwork, info)
        lwork = max(1, int(work(1)))
        deallocate(work)
        allocate(work(lwork))
        ac = a
        call dgeev('N', merge('V','N',o%retvec), n, ac, n, wr, wi, vl, 1, vr, n, work, lwork, info)
        if (info /= 0) then
            res%info = info
            return
        end if
        nr = n
        allocate(vals(nr), vecs(n,nr))
        call pack_dgeev(wr, wi, vr, o%retvec, vals, vecs)
        if (o%retvec) then
            call sort_complex_pairs(vals, vecs, w)
        else
            call sort_complex_pairs(vals, which=w)
        end if
        allocate(res%values(k))
        res%values = vals(1:k)
        if (o%retvec) then
            allocate(res%vectors(n,k))
            res%vectors = vecs(:,1:k)
        else
            allocate(res%vectors(n,0))
        end if
        res%nconv = k
        res%niter = 0
        res%nops = 0
        res%info = 0
    end function eigs_full

    subroutine pack_general_result(dr, di, z, nconv, retvec, res)
        real(dp), intent(in) :: dr(:), di(:), z(:,:)
        integer, intent(in) :: nconv
        logical, intent(in) :: retvec
        type(eigs_result), intent(inout) :: res
        integer :: i, n
        n = size(z,1)
        allocate(res%values(nconv))
        if (retvec) allocate(res%vectors(n,nconv))
        if (.not. retvec) allocate(res%vectors(n,0))
        i = 1
        do while (i <= nconv)
            res%values(i) = cmplx(dr(i), di(i), dp)
            if (retvec) then
                if (abs(di(i)) <= 100.0_dp * epsilon(1.0_dp)) then
                    res%vectors(:,i) = cmplx(z(:,i), 0.0_dp, dp)
                else if (di(i) > 0.0_dp .and. i < nconv) then
                    res%vectors(:,i) = cmplx(z(:,i), z(:,i+1), dp)
                    res%values(i+1) = cmplx(dr(i+1), di(i+1), dp)
                    res%vectors(:,i+1) = conjg(res%vectors(:,i))
                    i = i + 1
                else if (i > 1) then
                    res%vectors(:,i) = conjg(res%vectors(:,i-1))
                end if
            end if
            i = i + 1
        end do
    end subroutine pack_general_result

    subroutine pack_dgeev(wr, wi, vr, retvec, vals, vecs)
        real(dp), intent(in) :: wr(:), wi(:), vr(:,:)
        logical, intent(in) :: retvec
        complex(dp), intent(out) :: vals(:), vecs(:,:)
        integer :: i, n
        n = size(wr)
        vals = cmplx(wr, wi, dp)
        vecs = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. retvec) return
        i = 1
        do while (i <= n)
            if (abs(wi(i)) <= 100.0_dp * epsilon(1.0_dp)) then
                vecs(:,i) = cmplx(vr(:,i), 0.0_dp, dp)
            else if (wi(i) > 0.0_dp .and. i < n) then
                vecs(:,i) = cmplx(vr(:,i), vr(:,i+1), dp)
                vecs(:,i+1) = conjg(vecs(:,i))
                i = i + 1
            else if (i > 1) then
                vecs(:,i) = conjg(vecs(:,i-1))
            end if
            i = i + 1
        end do
    end subroutine pack_dgeev

    logical function valid_gen_which(w) result(ans)
        character(len=*), intent(in) :: w
        select case (trim(w))
        case ('LM','SM','LR','SR','LI','SI')
            ans = .true.
        case default
            ans = .false.
        end select
    end function valid_gen_which

end module rspectra_eigs_gen_mod
