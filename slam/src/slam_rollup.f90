! SPDX-License-Identifier: GPL-2.0-only
!
! Translation of slam's grouped rollup computations for numeric payloads.
module slam_rollup
    use iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use slam_kinds, only : dp
    use slam_utils, only : argsort_int64, linear_indices
    use slam_stm, only : simple_triplet_matrix, make_stm, simple_triplet_zero_matrix
    use slam_ssa, only : simple_sparse_array, make_ssa, simple_sparse_zero_array
    implicit none
    private

    public :: rollup_stm_sum
    public :: rollup_ssa_sum
    public :: factorize_groups

contains

    subroutine factorize_groups(groups, codes, nlevels)
        integer, intent(in) :: groups(:)
        integer, allocatable, intent(out) :: codes(:)
        integer, intent(out) :: nlevels
        integer(int64), allocatable :: key(:), uniq(:)
        integer, allocatable :: ord(:)
        integer :: k, m, q

        allocate(codes(size(groups)), key(size(groups)), uniq(size(groups)))
        codes = 0
        key = int(groups, int64)
        if (size(groups) == 0) then
            nlevels = 0
            return
        end if
        call argsort_int64(key, ord)
        m = 0
        do k = 1, size(ord)
            if (groups(ord(k)) == 0) cycle
            if (m == 0) then
                m = 1
                uniq(m) = int(groups(ord(k)), int64)
            else if (int(groups(ord(k)),int64) /= uniq(m)) then
                m = m + 1
                uniq(m) = int(groups(ord(k)), int64)
            end if
        end do
        nlevels = m
        do k = 1, size(groups)
            if (groups(k) == 0) cycle
            do q = 1, m
                if (int(groups(k),int64) == uniq(q)) then
                    codes(k) = q
                    exit
                end if
            end do
        end do
    end subroutine factorize_groups

    function rollup_stm_sum(x, margin, groups, na_rm, reduce_zeros) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        integer, intent(in) :: margin
        integer, intent(in) :: groups(:)
        logical, intent(in), optional :: na_rm, reduce_zeros
        type(simple_triplet_matrix) :: y
        integer, allocatable :: codes(:), ii(:), jj(:), ord(:), oi(:), oj(:)
        real(dp), allocatable :: vv(:), ov(:)
        integer(int64), allocatable :: key(:)
        logical :: rm, reduce
        integer :: ng, k, n, m, q, nr, nc
        integer(int64) :: current
        real(dp) :: value

        if (margin /= 1 .and. margin /= 2) error stop "rollup_stm_sum: margin must be 1 or 2"
        if (margin == 1) then
            if (size(groups) /= x%nrow) error stop "rollup_stm_sum: group length mismatch"
        else
            if (size(groups) /= x%ncol) error stop "rollup_stm_sum: group length mismatch"
        end if
        rm = .false.; reduce = .false.
        if (present(na_rm)) rm = na_rm
        if (present(reduce_zeros)) reduce = reduce_zeros
        call factorize_groups(groups, codes, ng)
        nr = x%nrow; nc = x%ncol
        if (margin == 1) nr = ng
        if (margin == 2) nc = ng
        if (ng == 0 .or. x%nnz() == 0) then
            y = simple_triplet_zero_matrix(nr,nc)
            return
        end if

        allocate(ii(x%nnz()), jj(x%nnz()), vv(x%nnz()), key(x%nnz()), &
                 oi(x%nnz()), oj(x%nnz()), ov(x%nnz()))
        n = 0
        do k = 1, x%nnz()
            if (margin == 1) then
                q = codes(x%i(k))
                if (q == 0) cycle
                n = n + 1; ii(n)=q; jj(n)=x%j(k); vv(n)=x%v(k)
            else
                q = codes(x%j(k))
                if (q == 0) cycle
                n = n + 1; ii(n)=x%i(k); jj(n)=q; vv(n)=x%v(k)
            end if
        end do
        if (n == 0) then
            y = simple_triplet_zero_matrix(nr,nc)
            return
        end if
        key(:n)=int(ii(:n),int64)+int(jj(:n)-1,int64)*int(nr,int64)
        call argsort_int64(key(:n),ord)
        m=0; k=1
        do while(k<=n)
            current=key(ord(k))
            value=0.0_dp
            do while(k<=n)
                if(key(ord(k))/=current) exit
                if (.not. (rm .and. ieee_is_nan(vv(ord(k))))) then
                    value=value+vv(ord(k))
                end if
                k=k+1
            end do
            if (.not. reduce .or. ieee_is_nan(value) .or. value/=0.0_dp) then
                m=m+1
                oi(m)=int(mod(current-1_int64,int(nr,int64)))+1
                oj(m)=int((current-1_int64)/int(nr,int64))+1
                ov(m)=value
            end if
        end do
        y=make_stm(oi(:m),oj(:m),ov(:m),nr,nc)
    end function rollup_stm_sum

    function rollup_ssa_sum(x, margin, groups, na_rm, reduce_zeros) result(y)
        type(simple_sparse_array), intent(in) :: x
        integer, intent(in) :: margin
        integer, intent(in) :: groups(:)
        logical, intent(in), optional :: na_rm, reduce_zeros
        type(simple_sparse_array) :: y
        integer, allocatable :: codes(:), coords(:,:), outcoords(:,:), dims(:), ord(:)
        real(dp), allocatable :: vals(:), outvals(:)
        integer(int64), allocatable :: key(:)
        logical :: rm, reduce
        integer :: ng,k,n,m,q
        integer(int64) :: current
        real(dp) :: value

        if(margin<1 .or. margin>x%rank()) error stop "rollup_ssa_sum: invalid margin"
        if(size(groups)/=x%dim(margin)) error stop "rollup_ssa_sum: group length mismatch"
        rm=.false.; reduce=.false.
        if(present(na_rm)) rm=na_rm
        if(present(reduce_zeros)) reduce=reduce_zeros
        call factorize_groups(groups,codes,ng)
        dims=x%dim; dims(margin)=ng
        if(ng==0 .or. x%nnz()==0) then
            y=simple_sparse_zero_array(dims)
            return
        end if
        allocate(coords(x%nnz(),x%rank()),vals(x%nnz()), &
                 outcoords(x%nnz(),x%rank()),outvals(x%nnz()))
        n=0
        do k=1,x%nnz()
            q=codes(x%i(k,margin))
            if(q==0) cycle
            n=n+1
            coords(n,:)=x%i(k,:)
            coords(n,margin)=q
            vals(n)=x%v(k)
        end do
        if(n==0) then
            y=simple_sparse_zero_array(dims)
            return
        end if
        key=linear_indices(coords(:n,:),dims)
        call argsort_int64(key,ord)
        m=0; k=1
        do while(k<=n)
            current=key(ord(k)); value=0.0_dp
            do while(k<=n)
                if(key(ord(k))/=current) exit
                if(.not.(rm .and. ieee_is_nan(vals(ord(k))))) value=value+vals(ord(k))
                k=k+1
            end do
            if(.not.reduce .or. ieee_is_nan(value) .or. value/=0.0_dp) then
                m=m+1
                outcoords(m,:)=coords(ord(k-1),:)
                outvals(m)=value
            end if
        end do
        y=make_ssa(outcoords(:m,:),outvals(:m),dims)
    end function rollup_ssa_sum

end module slam_rollup
