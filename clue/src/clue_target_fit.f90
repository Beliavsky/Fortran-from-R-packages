! SPDX-License-Identifier: GPL-2.0-only
module clue_target_fit
    use clue_kinds, only: dp
    use clue_pava, only: pava_mean, weighted_median
    implicit none
    private
    public :: fit_ultrametric_target_l2, fit_ultrametric_target_l1
contains
    recursive subroutine get_leaves(code, merge, leaves, nleaf)
        integer, intent(in) :: code
        integer, intent(in) :: merge(:,:)
        integer, intent(out) :: leaves(:)
        integer, intent(out) :: nleaf
        integer :: a(size(leaves)), b(size(leaves)), na, nb
        if (code < 0) then
            nleaf = 1
            leaves(1) = -code
        else
            call get_leaves(merge(code,1), merge, a, na)
            call get_leaves(merge(code,2), merge, b, nb)
            nleaf = na + nb
            leaves(1:na) = a(1:na)
            leaves(na+1:nleaf) = b(1:nb)
        end if
    end subroutine

    function fit_ultrametric_target_l2(d, merge, weights) result(u)
        real(dp), intent(in) :: d(:,:)
        integer, intent(in) :: merge(:,:)
        real(dp), intent(in), optional :: weights(:,:)
        real(dp), allocatable :: u(:,:), block_mean(:), block_w(:), fitted(:)
        integer, allocatable :: left(:), right(:)
        integer :: n, nm, i, a, b, nl, nr
        real(dp) :: sw, sx, wij
        n = size(d,1); nm = size(merge,1)
        allocate(u(n,n), block_mean(nm), block_w(nm), left(n), right(n))
        u = 0.0_dp
        do i = 1, nm
            call get_leaves(merge(i,1), merge, left, nl)
            call get_leaves(merge(i,2), merge, right, nr)
            sw = 0.0_dp; sx = 0.0_dp
            do a = 1, nl
                do b = 1, nr
                    wij = 1.0_dp
                    if (present(weights)) wij = weights(left(a),right(b))
                    sw = sw + wij
                    sx = sx + wij * d(left(a),right(b))
                end do
            end do
            block_w(i) = sw
            if (sw > 0.0_dp) then
                block_mean(i) = sx / sw
            else
                block_mean(i) = 0.0_dp
            end if
        end do
        fitted = pava_mean(block_mean, block_w)
        do i = 1, nm
            call get_leaves(merge(i,1), merge, left, nl)
            call get_leaves(merge(i,2), merge, right, nr)
            do a = 1, nl
                do b = 1, nr
                    u(left(a),right(b)) = fitted(i)
                    u(right(b),left(a)) = fitted(i)
                end do
            end do
        end do
    end function

    function fit_ultrametric_target_l1(d, merge, weights) result(u)
        real(dp), intent(in) :: d(:,:)
        integer, intent(in) :: merge(:,:)
        real(dp), intent(in), optional :: weights(:,:)
        real(dp), allocatable :: u(:,:), obs(:), wt(:), val(:), fitted(:)
        integer, allocatable :: left(:), right(:), bs(:), be(:), blo(:), bhi(:)
        integer, allocatable :: pi(:), pj(:), pblock(:)
        integer :: n, nm, np, i, a, b, nl, nr, pos, m, j
        n = size(d,1); nm = size(merge,1); np = n*(n-1)/2
        allocate(u(n,n),obs(np),wt(np),pi(np),pj(np),pblock(np))
        allocate(left(n),right(n),bs(nm),be(nm),val(nm),blo(nm),bhi(nm),fitted(nm))
        u=0.0_dp;pos=0
        do i=1,nm
            bs(i)=pos+1
            call get_leaves(merge(i,1),merge,left,nl)
            call get_leaves(merge(i,2),merge,right,nr)
            do a=1,nl
                do b=1,nr
                    pos=pos+1;pi(pos)=left(a);pj(pos)=right(b);pblock(pos)=i
                    obs(pos)=d(left(a),right(b));wt(pos)=1.0_dp
                    if(present(weights))wt(pos)=weights(left(a),right(b))
                end do
            end do
            be(i)=pos
            val(i)=weighted_median(obs(bs(i):be(i)),wt(bs(i):be(i)))
            blo(i)=i;bhi(i)=i
        end do
        m=nm;i=1
        do while(i<m)
            if(val(i)<=val(i+1))then
                i=i+1
            else
                bhi(i)=bhi(i+1)
                val(i)=weighted_median(obs(bs(blo(i)):be(bhi(i))),wt(bs(blo(i)):be(bhi(i))))
                do j=i+1,m-1
                    val(j)=val(j+1);blo(j)=blo(j+1);bhi(j)=bhi(j+1)
                end do
                m=m-1
                if(i>1)i=i-1
            end if
        end do
        do i=1,m
            fitted(blo(i):bhi(i))=val(i)
        end do
        do pos=1,np
            u(pi(pos),pj(pos))=fitted(pblock(pos))
            u(pj(pos),pi(pos))=fitted(pblock(pos))
        end do
    end function
end module clue_target_fit
