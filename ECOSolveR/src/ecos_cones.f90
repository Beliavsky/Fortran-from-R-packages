! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_cones
    use ecos_types, only : dp, ecos_problem
    use ecos_linalg, only : vecnorm2
    use ecos_sparse_cones, only : sparse_cone_slack, sparse_cone_values
    implicit none
    private
    public :: cone_scalar_eval, cone_slack, cone_dual_from_scalar, cone_violation

contains

    subroutine cone_slack(prob, x, s)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: s(:)
        if (size(s) == 0) return
        if (prob%sparse_storage) then
            call sparse_cone_slack(prob,x,s)
        else
            s = prob%h - matmul(prob%gmat,x)
        end if
    end subroutine cone_slack

    subroutine cone_scalar_eval(prob, x, lambda, g, jac, hlag)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: lambda(:)
        real(dp), intent(out) :: g(:), jac(:,:), hlag(:,:)
        real(dp), allocatable :: s(:), u(:), hu(:,:), qv(:), qs(:,:), tmp(:,:)
        real(dp) :: nr, a0, bb, cc, bs, cs, lg, epss
        integer :: row, ir, iq, k, qdim, idx, n, mscalar
        n = size(x)
        mscalar = prob%dims%scalar_inequalities()
        if (prob%sparse_storage) error stop 'cone_scalar_eval: dense routine called for sparse problem'
        g = 0.0_dp
        jac = 0.0_dp
        hlag = 0.0_dp
        if (mscalar == 0) return
        allocate(s(prob%ncone()))
        call cone_slack(prob,x,s)
        epss = 1.0e-12_dp
        idx = 0
        row = 0
        ! Positive orthant: s_i >= 0 -> -s_i <= 0.
        do ir = 1, prob%dims%l
            row = row + 1
            idx = idx + 1
            g(idx) = -s(row)
            jac(idx,:) = prob%gmat(row,:)
        end do
        ! SOC blocks: ||u|| - t <= 0.
        if (allocated(prob%dims%q)) then
            do iq = 1, size(prob%dims%q)
                qdim = prob%dims%q(iq)
                idx = idx + 1
                if (qdim <= 0) cycle
                a0 = s(row+1)
                if (qdim == 1) then
                    g(idx) = -a0
                    jac(idx,:) = prob%gmat(row+1,:)
                    row = row + 1
                else
                    allocate(u(qdim-1), hu(qdim-1,n))
                    u = s(row+2:row+qdim)
                    hu = prob%gmat(row+2:row+qdim,:)
                    nr = max(vecnorm2(u),epss)
                    g(idx) = nr - a0
                    jac(idx,:) = prob%gmat(row+1,:) - matmul(u,hu)/nr
                    if (lambda(idx) > 0.0_dp) then
                        allocate(qs(qdim-1,qdim-1), tmp(qdim-1,n))
                        qs = 0.0_dp
                        do k = 1, qdim-1
                            qs(k,k) = 1.0_dp/nr
                        end do
                        qs = qs - spread(u,2,qdim-1)*spread(u,1,qdim-1)/(nr**3)
                        tmp = matmul(qs,hu)
                        hlag = hlag + lambda(idx)*matmul(transpose(hu),tmp)
                        deallocate(qs,tmp)
                    end if
                    row = row + qdim
                    deallocate(u,hu)
                end if
            end do
        end if
        ! Exponential cones: b>0, c>0, a-c*log(b/c)<=0.
        do ir = 1, prob%dims%e
            a0 = s(row+1)
            bb = s(row+2)
            cc = s(row+3)
            bs = max(bb,epss)
            cs = max(cc,epss)
            idx = idx + 1
            g(idx) = -bb
            jac(idx,:) = prob%gmat(row+2,:)
            idx = idx + 1
            g(idx) = -cc
            jac(idx,:) = prob%gmat(row+3,:)
            idx = idx + 1
            lg = log(bs/cs)
            g(idx) = a0 - cc*lg
            allocate(qv(3))
            qv(1) = 1.0_dp
            qv(2) = -cc/bs
            qv(3) = 1.0_dp - lg
            jac(idx,:) = -matmul(qv,prob%gmat(row+1:row+3,:))
            if (lambda(idx) > 0.0_dp) then
                allocate(qs(3,3), tmp(3,n))
                qs = 0.0_dp
                qs(2,2) = max(cc,epss)/(bs*bs)
                qs(2,3) = -1.0_dp/bs
                qs(3,2) = qs(2,3)
                qs(3,3) = 1.0_dp/cs
                tmp = matmul(qs,prob%gmat(row+1:row+3,:))
                hlag = hlag + lambda(idx)*matmul(transpose(prob%gmat(row+1:row+3,:)),tmp)
                deallocate(qs,tmp)
            end if
            deallocate(qv)
            row = row + 3
        end do
    end subroutine cone_scalar_eval

    subroutine cone_dual_from_scalar(prob, x, lambda, z)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:), lambda(:)
        real(dp), intent(out) :: z(:)
        real(dp), allocatable :: s(:), u(:)
        real(dp) :: nr, bb, cc, bs, cs, lg
        integer :: row, idx, ir, iq, qdim
        z = 0.0_dp
        if (size(z) == 0) return
        allocate(s(prob%ncone()))
        call cone_slack(prob,x,s)
        row = 0
        idx = 0
        do ir = 1, prob%dims%l
            row = row+1
            idx = idx+1
            z(row) = lambda(idx)
        end do
        if (allocated(prob%dims%q)) then
            do iq = 1, size(prob%dims%q)
                qdim = prob%dims%q(iq)
                idx = idx+1
                if (qdim == 1) then
                    z(row+1) = lambda(idx)
                    row = row+1
                else
                    allocate(u(qdim-1))
                    u = s(row+2:row+qdim)
                    nr = max(vecnorm2(u),1.0e-14_dp)
                    z(row+1) = lambda(idx)
                    z(row+2:row+qdim) = -lambda(idx)*u/nr
                    row = row+qdim
                    deallocate(u)
                end if
            end do
        end if
        do ir = 1, prob%dims%e
            bb = s(row+2)
            cc = s(row+3)
            bs = max(bb,1.0e-14_dp)
            cs = max(cc,1.0e-14_dp)
            lg = log(bs/cs)
            z(row+1) = -lambda(idx+3)
            z(row+2) = lambda(idx+1) + lambda(idx+3)*cc/bs
            z(row+3) = lambda(idx+2) + lambda(idx+3)*(lg-1.0_dp)
            idx = idx+3
            row = row+3
        end do
    end subroutine cone_dual_from_scalar

    real(dp) function cone_violation(prob, x) result(v)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: g(:), jac(:,:), hess(:,:), lam(:)
        integer :: ni, n
        ni = prob%dims%scalar_inequalities()
        n = size(x)
        if (ni == 0) then
            v = 0.0_dp
            return
        end if
        allocate(g(ni))
        if (prob%sparse_storage) then
            call sparse_cone_values(prob,x,g)
        else
            allocate(jac(ni,n),hess(n,n),lam(ni))
            lam = 0.0_dp
            call cone_scalar_eval(prob,x,lam,g,jac,hess)
        end if
        v = max(0.0_dp,maxval(g))
    end function cone_violation

end module ecos_cones
