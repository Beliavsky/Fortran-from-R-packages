! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra_eigs_sym_mod
    use rspectra_kinds, only: dp
    use rspectra_types, only: eigs_opts, eigs_sym_result
    use rspectra_external, only: dsaupd, dseupd, dgetrf, dgetrs, dsyev
    use rspectra_operators, only: linear_operator, dense_operator, csr_operator, &
        make_dense_operator, csr_to_dense
    use rspectra_sort, only: sort_real_pairs
    implicit none
    private
    public :: eigs_sym_operator, eigs_sym_dense, eigs_sym_csr

contains

    function eigs_sym_operator(op, k, which, opts) result(res)
        class(linear_operator), intent(inout) :: op
        integer, intent(in) :: k
        character(len=*), intent(in), optional :: which
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_sym_result) :: res
        type(eigs_opts) :: o
        character(len=2) :: w
        integer :: n, ncv, ido, info, ierr, lworkl
        integer :: iparam(11), ipntr(11)
        real(dp), allocatable :: resid(:), v(:,:), workd(:), workl(:), d(:), z(:,:)
        logical, allocatable :: select(:)
        integer :: ix, iy, nconv

        n = op%nrow
        if (op%nrow /= op%ncol .or. n < 1 .or. k < 1 .or. k >= n) then
            res%info = -100
            return
        end if
        o = eigs_opts()
        if (present(opts)) o = opts
        w = 'LM'
        if (present(which)) w = adjustl(which)
        if (.not. valid_sym_which(w)) then
            res%info = -101
            return
        end if
        ncv = o%ncv
        if (ncv == 0) ncv = min(n, max(2 * k + 1, 20))
        if (ncv <= k .or. ncv > n) then
            res%info = -102
            return
        end if
        lworkl = ncv * ncv + 8 * ncv
        allocate(resid(n), v(n,ncv), workd(3*n), workl(lworkl))
        allocate(d(k), z(n,k), select(ncv))
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
            call dsaupd(ido, 'I', n, w, k, o%tol, resid, ncv, v, n, iparam, ipntr, &
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
        call dseupd(o%retvec, 'A', select, d, z, n, 0.0_dp, 'I', n, w, k, o%tol, &
            resid, ncv, v, n, iparam, ipntr, workd, workl, lworkl, ierr)
        if (ierr /= 0) then
            res%info = ierr
            return
        end if
        allocate(res%values(nconv))
        res%values = d(1:nconv)
        if (o%retvec) then
            allocate(res%vectors(n,nconv))
            res%vectors = z(:,1:nconv)
            call sort_real_pairs(res%values, res%vectors, w)
        else
            allocate(res%vectors(n,0))
            call sort_real_pairs(res%values, which=w)
        end if
        res%info = info
    end function eigs_sym_operator

    function eigs_sym_dense(a, k, which, sigma, opts) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        character(len=*), intent(in), optional :: which
        real(dp), intent(in), optional :: sigma
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_sym_result) :: res
        type(dense_operator) :: op
        character(len=2) :: w
        integer :: n

        n = size(a,1)
        if (size(a,2) /= n .or. k < 1 .or. k > n) then
            res%info = -100
            return
        end if
        w = 'LM'
        if (present(which)) w = adjustl(which)
        if (.not. valid_sym_which(w)) then
            res%info = -101
            return
        end if
        if (k == n) then
            res = eigs_sym_full(a, k, w, opts)
            return
        end if
        if (present(sigma)) then
            res = eigs_sym_shift_dense(a, k, w, sigma, opts)
        else
            op = make_dense_operator(a)
            res = eigs_sym_operator(op, k, w, opts)
        end if
    end function eigs_sym_dense

    function eigs_sym_csr(op, k, which, sigma, opts) result(res)
        type(csr_operator), intent(inout) :: op
        integer, intent(in) :: k
        character(len=*), intent(in), optional :: which
        real(dp), intent(in), optional :: sigma
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_sym_result) :: res
        real(dp), allocatable :: a(:,:)
        character(len=2) :: w
        w = 'LM'
        if (present(which)) w = adjustl(which)
        if (present(sigma)) then
            a = csr_to_dense(op)
            res = eigs_sym_dense(a, k, w, sigma, opts)
        else
            res = eigs_sym_operator(op, k, w, opts)
        end if
    end function eigs_sym_csr

    function eigs_sym_shift_dense(a, k, w, sigma, opts) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        character(len=2), intent(in) :: w
        real(dp), intent(in) :: sigma
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_sym_result) :: res
        type(eigs_opts) :: o
        integer :: n, ncv, ido, info, ierr, lworkl, i
        integer :: iparam(11), ipntr(11)
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: lu(:,:), resid(:), v(:,:), workd(:), workl(:), d(:), z(:,:), rhs(:,:)
        logical, allocatable :: select(:)
        integer :: ix, iy, nconv

        n = size(a,1)
        o = eigs_opts()
        if (present(opts)) o = opts
        ncv = o%ncv
        if (ncv == 0) ncv = min(n, max(2 * k + 1, 20))
        if (ncv <= k .or. ncv > n) then
            res%info = -102
            return
        end if
        allocate(lu(n,n), ipiv(n), rhs(n,1))
        lu = a
        do i = 1, n
            lu(i,i) = lu(i,i) - sigma
        end do
        call dgetrf(n, n, lu, n, ipiv, ierr)
        if (ierr /= 0) then
            res%info = -200 - ierr
            return
        end if
        lworkl = ncv * ncv + 8 * ncv
        allocate(resid(n), v(n,ncv), workd(3*n), workl(lworkl))
        allocate(d(k), z(n,k), select(ncv))
        resid = 0.0_dp
        v = 0.0_dp
        workd = 0.0_dp
        workl = 0.0_dp
        iparam = 0
        ipntr = 0
        iparam(1) = 1
        iparam(3) = o%maxitr
        iparam(7) = 3
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
            call dsaupd(ido, 'I', n, w, k, o%tol, resid, ncv, v, n, iparam, ipntr, &
                workd, workl, lworkl, info)
            if (ido == -1 .or. ido == 1) then
                ix = ipntr(1)
                iy = ipntr(2)
                rhs(:,1) = workd(ix:ix+n-1)
                call dgetrs('N', n, 1, lu, n, ipiv, rhs, n, ierr)
                if (ierr /= 0) then
                    res%info = -300 - ierr
                    return
                end if
                workd(iy:iy+n-1) = rhs(:,1)
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
        call dseupd(o%retvec, 'A', select, d, z, n, sigma, 'I', n, w, k, o%tol, &
            resid, ncv, v, n, iparam, ipntr, workd, workl, lworkl, ierr)
        if (ierr /= 0) then
            res%info = ierr
            return
        end if
        allocate(res%values(nconv))
        res%values = d(1:nconv)
        if (o%retvec) then
            allocate(res%vectors(n,nconv))
            res%vectors = z(:,1:nconv)
            call sort_real_pairs(res%values, res%vectors, w, sigma, .true.)
        else
            allocate(res%vectors(n,0))
            call sort_real_pairs(res%values, which=w, sigma=sigma, shift_mode=.true.)
        end if
        res%info = info
    end function eigs_sym_shift_dense

    function eigs_sym_full(a, k, w, opts) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        character(len=2), intent(in) :: w
        type(eigs_opts), intent(in), optional :: opts
        type(eigs_sym_result) :: res
        type(eigs_opts) :: o
        real(dp), allocatable :: ac(:,:), eval(:), work(:)
        integer :: n, lwork, info
        n = size(a,1)
        o = eigs_opts()
        if (present(opts)) o = opts
        allocate(ac(n,n), eval(n), work(1))
        ac = a
        lwork = -1
        call dsyev(merge('V','N',o%retvec), 'L', n, ac, n, eval, work, lwork, info)
        lwork = max(1, int(work(1)))
        deallocate(work)
        allocate(work(lwork))
        call dsyev(merge('V','N',o%retvec), 'L', n, ac, n, eval, work, lwork, info)
        if (info /= 0) then
            res%info = info
            return
        end if
        allocate(res%values(n))
        res%values = eval
        if (o%retvec) then
            allocate(res%vectors(n,n))
            res%vectors = ac
            call sort_real_pairs(res%values, res%vectors, w)
        else
            allocate(res%vectors(n,0))
            call sort_real_pairs(res%values, which=w)
        end if
        if (k < n) then
            res%values = res%values(1:k)
            if (o%retvec) res%vectors = res%vectors(:,1:k)
        end if
        res%nconv = k
        res%niter = 0
        res%nops = 0
        res%info = 0
    end function eigs_sym_full

    logical function valid_sym_which(w) result(ans)
        character(len=*), intent(in) :: w
        select case (trim(w))
        case ('LM','SM','LA','SA','BE')
            ans = .true.
        case default
            ans = .false.
        end select
    end function valid_sym_which

end module rspectra_eigs_sym_mod
