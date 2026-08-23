! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra_svds_mod
    use rspectra_kinds, only: dp
    use rspectra_types, only: eigs_opts, eigs_sym_result, svds_opts, svds_result
    use rspectra_external, only: dgesvd
    use rspectra_operators, only: linear_operator, dense_operator, csr_operator, make_dense_operator
    use rspectra_eigs_sym_mod, only: eigs_sym_operator
    implicit none
    private
    public :: svds_operator, svds_dense, svds_csr

    type, extends(linear_operator) :: normal_operator
        class(linear_operator), pointer :: base => null()
        logical :: tall = .true.
        real(dp), allocatable :: ctr(:), scl(:)
        real(dp), allocatable :: wm(:), wn(:)
    contains
        procedure :: prod => normal_prod
        procedure :: tprod => normal_tprod
        procedure :: apply_b
        procedure :: apply_bt
    end type normal_operator

contains

    function svds_operator(op, k, nu, nv, opts, center, scale) result(res)
        class(linear_operator), target, intent(inout) :: op
        integer, intent(in) :: k
        integer, intent(in), optional :: nu, nv
        type(svds_opts), intent(in), optional :: opts
        real(dp), intent(in), optional :: center(:), scale(:)
        type(svds_result) :: res
        type(svds_opts) :: so
        type(eigs_opts) :: eo
        type(eigs_sym_result) :: er
        type(normal_operator) :: nop
        integer :: m, n, wd, nru, nrv, i
        real(dp) :: tiny_d
        real(dp), allocatable :: tmpm(:), tmpn(:)

        m = op%nrow
        n = op%ncol
        wd = min(m,n)
        if (wd < 3 .or. k < 1 .or. k >= wd) then
            res%info = -100
            return
        end if
        nru = k
        nrv = k
        if (present(nu)) nru = nu
        if (present(nv)) nrv = nv
        if (nru < 0 .or. nrv < 0 .or. nru > k .or. nrv > k) then
            res%info = -101
            return
        end if
        so = svds_opts()
        if (present(opts)) so = opts
        eo%ncv = so%ncv
        eo%maxitr = so%maxitr
        eo%tol = so%tol
        eo%retvec = .true.

        nop%base => op
        nop%tall = m > n
        nop%nrow = wd
        nop%ncol = wd
        allocate(nop%ctr(n), nop%scl(n), nop%wm(m), nop%wn(n))
        nop%ctr = 0.0_dp
        nop%scl = 1.0_dp
        if (present(center)) then
            if (size(center) /= n) then
                res%info = -102
                return
            end if
            nop%ctr = center
        end if
        if (present(scale)) then
            if (size(scale) /= n .or. any(scale <= 0.0_dp)) then
                res%info = -103
                return
            end if
            nop%scl = scale
        end if

        er = eigs_sym_operator(nop, k, 'LA', eo)
        res%info = er%info
        res%nconv = er%nconv
        res%niter = er%niter
        res%nops = 2 * er%nops
        if (er%nconv <= 0) then
            allocate(res%d(0), res%u(m,0), res%v(n,0))
            return
        end if
        allocate(res%d(er%nconv))
        res%d = sqrt(max(er%values, 0.0_dp))
        nru = min(nru, er%nconv)
        nrv = min(nrv, er%nconv)
        allocate(res%u(m,nru), res%v(n,nrv))
        tiny_d = sqrt(epsilon(1.0_dp))
        allocate(tmpm(m), tmpn(n))
        if (nop%tall) then
            if (nrv > 0) res%v = er%vectors(:,1:nrv)
            do i = 1, nru
                if (res%d(i) > tiny_d) then
                    call nop%apply_b(er%vectors(:,i), tmpm)
                    res%u(:,i) = tmpm / res%d(i)
                else
                    res%u(:,i) = 0.0_dp
                end if
                res%nops = res%nops + 1
            end do
        else
            if (nru > 0) res%u = er%vectors(:,1:nru)
            do i = 1, nrv
                if (res%d(i) > tiny_d) then
                    call nop%apply_bt(er%vectors(:,i), tmpn)
                    res%v(:,i) = tmpn / res%d(i)
                else
                    res%v(:,i) = 0.0_dp
                end if
                res%nops = res%nops + 1
            end do
        end if
    end function svds_operator

    function svds_dense(a, k, nu, nv, opts, center, scale, center_vec, scale_vec) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        integer, intent(in), optional :: nu, nv
        type(svds_opts), intent(in), optional :: opts
        logical, intent(in), optional :: center, scale
        real(dp), intent(in), optional :: center_vec(:), scale_vec(:)
        type(svds_result) :: res
        type(dense_operator) :: op
        integer :: m, n, wd, j
        real(dp), allocatable :: ctr(:), scl(:)
        logical :: do_center, do_scale

        m = size(a,1)
        n = size(a,2)
        wd = min(m,n)
        if (wd < 1 .or. k < 1 .or. k > wd) then
            res%info = -100
            return
        end if
        allocate(ctr(n), scl(n))
        ctr = 0.0_dp
        scl = 1.0_dp
        do_center = .false.
        do_scale = .false.
        if (present(center)) do_center = center
        if (present(scale)) do_scale = scale
        if (present(center_vec)) then
            if (size(center_vec) /= n) then
                res%info = -102
                return
            end if
            ctr = center_vec
            do_center = .true.
        else if (do_center) then
            ctr = sum(a, dim=1) / real(m,dp)
        end if
        if (present(scale_vec)) then
            if (size(scale_vec) /= n .or. any(scale_vec <= 0.0_dp)) then
                res%info = -103
                return
            end if
            scl = scale_vec
            do_scale = .true.
        else if (do_scale) then
            do j = 1, n
                scl(j) = sqrt(sum((a(:,j) - ctr(j))**2))
            end do
            if (any(scl <= 0.0_dp)) then
                res%info = -104
                return
            end if
        end if
        if (.not. do_center) ctr = 0.0_dp
        if (.not. do_scale) scl = 1.0_dp
        if (k == wd) then
            res = svds_full_dense(a, k, nu, nv, ctr, scl)
            return
        end if
        op = make_dense_operator(a)
        res = svds_operator(op, k, nu, nv, opts, ctr, scl)
    end function svds_dense

    function svds_csr(op, k, nu, nv, opts, center, scale) result(res)
        type(csr_operator), target, intent(inout) :: op
        integer, intent(in) :: k
        integer, intent(in), optional :: nu, nv
        type(svds_opts), intent(in), optional :: opts
        real(dp), intent(in), optional :: center(:), scale(:)
        type(svds_result) :: res
        res = svds_operator(op, k, nu, nv, opts, center, scale)
    end function svds_csr

    subroutine normal_prod(self, x, y)
        class(normal_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        if (self%tall) then
            call self%apply_b(x, self%wm)
            call self%apply_bt(self%wm, y)
        else
            call self%apply_bt(x, self%wn)
            call self%apply_b(self%wn, y)
        end if
    end subroutine normal_prod

    subroutine normal_tprod(self, x, y)
        class(normal_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        call self%prod(x, y)
    end subroutine normal_tprod

    subroutine apply_b(self, x, y)
        class(normal_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        self%wn = x / self%scl
        call self%base%prod(self%wn, y)
        y = y - dot_product(self%ctr, self%wn)
    end subroutine apply_b

    subroutine apply_bt(self, x, y)
        class(normal_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        call self%base%tprod(x, y)
        y = (y - self%ctr * sum(x)) / self%scl
    end subroutine apply_bt

    function svds_full_dense(a, k, nu, nv, ctr, scl) result(res)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        integer, intent(in), optional :: nu, nv
        real(dp), intent(in) :: ctr(:), scl(:)
        type(svds_result) :: res
        real(dp), allocatable :: b(:,:), s(:), ufull(:,:), vt(:,:), work(:)
        integer :: m, n, wd, lwork, info, nru, nrv, j
        m = size(a,1)
        n = size(a,2)
        wd = min(m,n)
        nru = k
        nrv = k
        if (present(nu)) nru = nu
        if (present(nv)) nrv = nv
        if (nru < 0 .or. nrv < 0 .or. nru > k .or. nrv > k) then
            res%info = -101
            return
        end if
        allocate(b(m,n), s(wd), ufull(m,wd), vt(wd,n), work(1))
        do j = 1, n
            b(:,j) = (a(:,j) - ctr(j)) / scl(j)
        end do
        lwork = -1
        call dgesvd('S','S',m,n,b,m,s,ufull,m,vt,wd,work,lwork,info)
        lwork = max(1,int(work(1)))
        deallocate(work)
        allocate(work(lwork))
        do j = 1, n
            b(:,j) = (a(:,j) - ctr(j)) / scl(j)
        end do
        call dgesvd('S','S',m,n,b,m,s,ufull,m,vt,wd,work,lwork,info)
        if (info /= 0) then
            res%info = info
            return
        end if
        allocate(res%d(k), res%u(m,nru), res%v(n,nrv))
        res%d = s(1:k)
        if (nru > 0) res%u = ufull(:,1:nru)
        if (nrv > 0) res%v = transpose(vt(1:nrv,:))
        res%nconv = k
        res%niter = 0
        res%nops = 0
        res%info = 0
    end function svds_full_dense

end module rspectra_svds_mod
