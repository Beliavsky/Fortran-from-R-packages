! SPDX-License-Identifier: GPL-2.0-only
!
! Numeric translation of slam's simple_sparse_array representation.
module slam_ssa
    use iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use slam_kinds, only : dp
    use slam_utils, only : linear_indices, coordinates_from_linear, product_dims, argsort_int64
    use slam_stm, only : simple_triplet_matrix, make_stm
    implicit none
    private

    type, public :: simple_sparse_array
        integer, allocatable :: i(:,:)  ! (nnz, rank), 1-based coordinates
        real(dp), allocatable :: v(:)
        integer, allocatable :: dim(:)
    contains
        procedure :: valid => ssa_valid
        procedure :: nnz => ssa_nnz
        procedure :: rank => ssa_rank
        procedure :: to_dense_flat => ssa_to_dense_flat
        procedure :: permute => ssa_permute
        procedure :: reshape => ssa_reshape
        procedure :: mean => ssa_mean
    end type simple_sparse_array

    public :: make_ssa
    public :: dense_flat_to_ssa
    public :: simple_sparse_zero_array
    public :: reduce_simple_sparse_array
    public :: drop_simple_sparse_array
    public :: extend_simple_sparse_array
    public :: abind_simple_sparse_array
    public :: ssa_set
    public :: ssa_to_stm
    public :: stm_to_ssa

contains

    function make_ssa(i, v, dim, check_duplicates) result(x)
        integer, intent(in) :: i(:,:)
        real(dp), intent(in) :: v(:)
        integer, intent(in) :: dim(:)
        logical, intent(in), optional :: check_duplicates
        type(simple_sparse_array) :: x
        logical :: check
        integer(int64), allocatable :: key(:)

        if (size(i,1) /= size(v)) error stop "make_ssa: i and v lengths do not conform"
        if (size(i,2) /= size(dim)) error stop "make_ssa: i and dim ranks do not conform"
        if (any(dim < 0)) error stop "make_ssa: negative dimension"
        if (any(dim == 0) .and. size(v) > 0) error stop "make_ssa: nonzero entries with zero dimension"
        if (size(i,1) > 0) then
            key = linear_indices(i, dim)
            if (any(key == 0_int64)) error stop "make_ssa: coordinate out of bounds"
        end if
        allocate(x%i(size(i,1), size(i,2)), x%v(size(v)), x%dim(size(dim)))
        x%i = i
        x%v = v
        x%dim = dim
        check = .true.
        if (present(check_duplicates)) check = check_duplicates
        if (check) then
            if (has_duplicate_coords(x)) error stop "make_ssa: duplicate coordinate rows are not allowed"
        end if
    end function make_ssa

    pure function ssa_nnz(self) result(n)
        class(simple_sparse_array), intent(in) :: self
        integer :: n
        if (allocated(self%v)) then
            n = size(self%v)
        else
            n = 0
        end if
    end function ssa_nnz

    pure function ssa_rank(self) result(n)
        class(simple_sparse_array), intent(in) :: self
        integer :: n
        if (allocated(self%dim)) then
            n = size(self%dim)
        else
            n = 0
        end if
    end function ssa_rank

    function ssa_valid(self) result(ok)
        class(simple_sparse_array), intent(in) :: self
        logical :: ok
        integer(int64), allocatable :: key(:)

        ok = allocated(self%i) .and. allocated(self%v) .and. allocated(self%dim)
        if (.not. ok) return
        ok = size(self%i,1) == size(self%v) .and. size(self%i,2) == size(self%dim)
        if (.not. ok) return
        ok = all(self%dim >= 0)
        if (.not. ok) return
        if (any(self%dim == 0) .and. self%nnz() > 0) then
            ok = .false.
            return
        end if
        if (self%nnz() > 0) then
            key = linear_indices(self%i, self%dim)
            ok = all(key > 0_int64)
            if (.not. ok) return
        end if
        ok = .not. has_duplicate_coords(self)
    end function ssa_valid

    function has_duplicate_coords(x) result(dup)
        class(simple_sparse_array), intent(in) :: x
        logical :: dup
        integer(int64), allocatable :: key(:)
        integer, allocatable :: order(:)
        integer :: k

        dup = .false.
        if (x%nnz() < 2) return
        key = linear_indices(x%i, x%dim)
        call argsort_int64(key, order)
        do k = 2, size(order)
            if (key(order(k)) == key(order(k-1))) then
                dup = .true.
                return
            end if
        end do
    end function has_duplicate_coords

    function ssa_to_dense_flat(self) result(a)
        class(simple_sparse_array), intent(in) :: self
        real(dp), allocatable :: a(:)
        integer(int64), allocatable :: key(:)
        integer(int64) :: n
        integer :: k

        n = product_dims(self%dim)
        if (n < 0_int64 .or. n > int(huge(1), int64)) &
            error stop "ssa_to_dense_flat: dense array too large"
        allocate(a(int(n)))
        a = 0.0_dp
        if (self%nnz() == 0) return
        key = linear_indices(self%i, self%dim)
        do k = 1, self%nnz()
            a(int(key(k))) = self%v(k)
        end do
    end function ssa_to_dense_flat

    function dense_flat_to_ssa(a, dim) result(x)
        real(dp), intent(in) :: a(:)
        integer, intent(in) :: dim(:)
        type(simple_sparse_array) :: x
        integer(int64) :: n
        integer :: k, m
        integer, allocatable :: coords(:,:), c(:)
        real(dp), allocatable :: vv(:)

        n = product_dims(dim)
        if (n /= int(size(a), int64)) error stop "dense_flat_to_ssa: data length does not equal product(dim)"
        m = 0
        do k = 1, size(a)
            if (ieee_is_nan(a(k)) .or. a(k) /= 0.0_dp) m = m + 1
        end do
        allocate(coords(m,size(dim)), vv(m), c(size(dim)))
        m = 0
        do k = 1, size(a)
            if (ieee_is_nan(a(k)) .or. a(k) /= 0.0_dp) then
                m = m + 1
                call coordinates_from_linear(int(k,int64), dim, c)
                coords(m,:) = c
                vv(m) = a(k)
            end if
        end do
        x = make_ssa(coords, vv, dim)
    end function dense_flat_to_ssa

    function simple_sparse_zero_array(dim) result(x)
        integer, intent(in) :: dim(:)
        type(simple_sparse_array) :: x
        integer, allocatable :: ii(:,:)
        real(dp), allocatable :: vv(:)
        allocate(ii(0,size(dim)), vv(0))
        x = make_ssa(ii, vv, dim)
    end function simple_sparse_zero_array

    function reduce_simple_sparse_array(x, strict, order_entries) result(y)
        type(simple_sparse_array), intent(in) :: x
        logical, intent(in), optional :: strict, order_entries
        type(simple_sparse_array) :: y
        logical :: s, ord, same
        integer(int64), allocatable :: key(:)
        integer, allocatable :: order(:), ii(:,:)
        real(dp), allocatable :: vv(:)
        integer :: m, k, q, first
        real(dp) :: value

        s = .false.; ord = .false.
        if (present(strict)) s = strict
        if (present(order_entries)) ord = order_entries
        if (x%nnz() == 0) then
            y = simple_sparse_zero_array(x%dim)
            return
        end if
        key = linear_indices(x%i, x%dim)
        call argsort_int64(key, order)
        allocate(ii(x%nnz(),x%rank()), vv(x%nnz()))
        m = 0
        k = 1
        do while (k <= x%nnz())
            first = order(k)
            value = x%v(first)
            same = .true.
            q = k + 1
            do while (q <= x%nnz())
                if (key(order(q)) /= key(first)) exit
                if (.not. values_equal(value, x%v(order(q)))) same = .false.
                q = q + 1
            end do
            if (q - k > 1 .and. s) error stop "reduce_simple_sparse_array: multiple entries"
            if (.not. same) value = ieee_value(0.0_dp, ieee_quiet_nan)
            if (.not. ieee_is_nan(value) .and. value == 0.0_dp) then
                if (s) error stop "reduce_simple_sparse_array: zero entries"
            else
                m = m + 1
                ii(m,:) = x%i(first,:)
                vv(m) = value
            end if
            k = q
        end do
        if (ord) then
            y = make_ssa(ii(:m,:), vv(:m), x%dim)
        else
            ! R's reduction preserves first-occurrence order unless order=TRUE.
            ! Reconstruct that order using each group's first original index.
            y = make_ssa(ii(:m,:), vv(:m), x%dim)
            call reorder_by_original_first(y, x)
        end if
    end function reduce_simple_sparse_array

    subroutine reorder_by_original_first(y, x)
        type(simple_sparse_array), intent(inout) :: y
        type(simple_sparse_array), intent(in) :: x
        integer(int64), allocatable :: ky(:), kx(:), firstkey(:)
        integer, allocatable :: p(:), ii(:,:)
        real(dp), allocatable :: vv(:)
        integer :: a, b, first

        if (y%nnz() < 2) return
        ky = linear_indices(y%i, y%dim)
        kx = linear_indices(x%i, x%dim)
        allocate(firstkey(y%nnz()))
        do a = 1, y%nnz()
            first = huge(1)
            do b = 1, x%nnz()
                if (kx(b) == ky(a)) then
                    first = b
                    exit
                end if
            end do
            firstkey(a) = int(first, int64)
        end do
        call argsort_int64(firstkey, p)
        allocate(ii(y%nnz(),y%rank()), vv(y%nnz()))
        ii = y%i(p,:)
        vv = y%v(p)
        call move_alloc(ii, y%i)
        call move_alloc(vv, y%v)
    end subroutine reorder_by_original_first

    function drop_simple_sparse_array(x) result(y)
        type(simple_sparse_array), intent(in) :: x
        type(simple_sparse_array) :: y
        integer, allocatable :: keep(:), dims(:), ii(:,:)
        integer :: k, n

        if (any(x%dim == 0)) then
            y = simple_sparse_zero_array([0])
            return
        end if
        n = count(x%dim /= 1)
        if (n == 0) then
            allocate(ii(x%nnz(),1), dims(1))
            dims = 1
            if (x%nnz() > 0) ii(:,1) = 1
            y = make_ssa(ii, x%v, dims)
            return
        end if
        allocate(keep(n), dims(n), ii(x%nnz(),n))
        n = 0
        do k = 1, x%rank()
            if (x%dim(k) /= 1) then
                n = n + 1
                keep(n) = k
                dims(n) = x%dim(k)
                if (x%nnz() > 0) ii(:,n) = x%i(:,k)
            end if
        end do
        y = make_ssa(ii, x%v, dims)
    end function drop_simple_sparse_array

    function ssa_reshape(self, dim) result(y)
        class(simple_sparse_array), intent(in) :: self
        integer, intent(in) :: dim(:)
        type(simple_sparse_array) :: y
        integer(int64), allocatable :: idx(:)
        integer, allocatable :: coords(:,:)
        integer :: k
        logical :: ok

        if (size(dim) == 0) error stop "ssa_reshape: dim must have positive length"
        if (any(dim < 0)) error stop "ssa_reshape: negative dimension"
        if (product_dims(dim) /= product_dims(self%dim)) &
            error stop "ssa_reshape: element count must not change"
        idx = linear_indices(self%i,self%dim)
        allocate(coords(self%nnz(),size(dim)))
        do k=1,self%nnz()
            call coordinates_from_linear(idx(k),dim,coords(k,:),ok)
            if (.not. ok) error stop "ssa_reshape: internal index conversion failed"
        end do
        y=make_ssa(coords,self%v,dim)
    end function ssa_reshape

    function ssa_permute(self, perm) result(y)
        class(simple_sparse_array), intent(in) :: self
        integer, intent(in) :: perm(:)
        type(simple_sparse_array) :: y
        logical, allocatable :: seen(:)
        integer :: k

        if (size(perm) /= self%rank()) error stop "ssa_permute: invalid permutation length"
        allocate(seen(self%rank()))
        seen = .false.
        do k = 1, self%rank()
            if (perm(k) < 1 .or. perm(k) > self%rank()) error stop "ssa_permute: invalid permutation"
            if (seen(perm(k))) error stop "ssa_permute: duplicate permutation entry"
            seen(perm(k)) = .true.
        end do
        y = make_ssa(self%i(:,perm), self%v, self%dim(perm))
    end function ssa_permute

    function extend_simple_sparse_array(x, margin) result(y)
        type(simple_sparse_array), intent(in) :: x
        integer, intent(in) :: margin
        type(simple_sparse_array) :: y
        integer :: before, r
        integer, allocatable :: dims(:), ii(:,:)

        if (margin < 0) then
            before = -margin - 1
        else
            before = margin
        end if
        r = x%rank()
        if (before < 0 .or. before > r) error stop "extend_simple_sparse_array: invalid margin"
        allocate(dims(r+1), ii(x%nnz(),r+1))
        if (before > 0) then
            dims(:before)=x%dim(:before)
            if (x%nnz()>0) ii(:,:before)=x%i(:,:before)
        end if
        dims(before+1)=1
        if (x%nnz()>0) ii(:,before+1)=1
        if (before < r) then
            dims(before+2:)=x%dim(before+1:)
            if (x%nnz()>0) ii(:,before+2:)=x%i(:,before+1:)
        end if
        y = make_ssa(ii, x%v, dims)
    end function extend_simple_sparse_array

    function abind_simple_sparse_array(xin, yin, margin) result(z)
        type(simple_sparse_array), intent(in) :: xin, yin
        integer, intent(in) :: margin
        type(simple_sparse_array) :: z
        type(simple_sparse_array) :: x, y
        integer :: m, nx, ny, k
        integer, allocatable :: dims(:), ii(:,:)
        real(dp), allocatable :: vv(:)

        if (margin == 0) error stop "abind_simple_sparse_array: margin cannot be zero"
        x = xin; y = yin
        if (margin < 0) then
            x = extend_simple_sparse_array(x, margin)
            y = extend_simple_sparse_array(y, margin)
            m = abs(margin)
        else
            m = margin
        end if
        if (x%rank() /= y%rank()) error stop "abind_simple_sparse_array: ranks do not conform"
        if (m < 1 .or. m > x%rank()) error stop "abind_simple_sparse_array: invalid margin"
        do k = 1, x%rank()
            if (k /= m .and. x%dim(k) /= y%dim(k)) &
                error stop "abind_simple_sparse_array: common dimensions do not conform"
        end do
        dims = x%dim
        dims(m) = dims(m) + y%dim(m)
        nx=x%nnz(); ny=y%nnz()
        allocate(ii(nx+ny,x%rank()), vv(nx+ny))
        if (nx>0) then
            ii(:nx,:)=x%i; vv(:nx)=x%v
        end if
        if (ny>0) then
            ii(nx+1:,:)=y%i
            ii(nx+1:,m)=ii(nx+1:,m)+x%dim(m)
            vv(nx+1:)=y%v
        end if
        z=make_ssa(ii,vv,dims)
    end function abind_simple_sparse_array

    subroutine ssa_set(x, coords, values)
        type(simple_sparse_array), intent(inout) :: x
        integer, intent(in) :: coords(:,:)
        real(dp), intent(in) :: values(:)
        integer(int64), allocatable :: oldkey(:), newkey(:)
        integer, allocatable :: ii(:,:)
        real(dp), allocatable :: vv(:)
        integer :: k, p, n, cap
        logical :: found
        real(dp) :: value

        if (size(coords,2) /= x%rank()) error stop "ssa_set: coordinate rank mismatch"
        if (size(values)==0 .and. size(coords,1)>0) error stop "ssa_set: empty replacement"
        newkey=linear_indices(coords,x%dim)
        if (any(newkey==0_int64)) error stop "ssa_set: coordinate out of bounds"
        oldkey=linear_indices(x%i,x%dim)
        cap=x%nnz()+size(coords,1)
        allocate(ii(cap,x%rank()),vv(cap))
        n=x%nnz()
        if(n>0) then
            ii(:n,:)=x%i; vv(:n)=x%v
        end if
        do k=1,size(coords,1)
            value=values(mod(k-1,size(values))+1)
            found=.false.
            do p=1,n
                if (all(ii(p,:)==coords(k,:))) then
                    vv(p)=value; found=.true.; exit
                end if
            end do
            if(.not.found) then
                n=n+1; ii(n,:)=coords(k,:); vv(n)=value
            end if
        end do
        x=make_ssa(ii(:n,:),vv(:n),x%dim)
        call remove_ssa_zeros(x)
    end subroutine ssa_set

    subroutine remove_ssa_zeros(x)
        type(simple_sparse_array), intent(inout) :: x
        integer, allocatable :: ii(:,:)
        real(dp), allocatable :: vv(:)
        integer :: k,n
        n=0
        do k=1,x%nnz()
            if (ieee_is_nan(x%v(k)) .or. x%v(k)/=0.0_dp) n=n+1
        end do
        if(n==x%nnz()) return
        allocate(ii(n,x%rank()),vv(n)); n=0
        do k=1,x%nnz()
            if (ieee_is_nan(x%v(k)) .or. x%v(k)/=0.0_dp) then
                n=n+1; ii(n,:)=x%i(k,:); vv(n)=x%v(k)
            end if
        end do
        call move_alloc(ii,x%i); call move_alloc(vv,x%v)
    end subroutine remove_ssa_zeros

    function ssa_to_stm(x) result(m)
        type(simple_sparse_array), intent(in) :: x
        type(simple_triplet_matrix) :: m
        integer, allocatable :: j(:)

        select case(x%rank())
        case(1)
            allocate(j(x%nnz())); j=1
            m=make_stm(x%i(:,1),j,x%v,x%dim(1),1)
        case(2)
            m=make_stm(x%i(:,1),x%i(:,2),x%v,x%dim(1),x%dim(2))
        case default
            error stop "ssa_to_stm: only rank 1 or 2 is supported"
        end select
    end function ssa_to_stm

    function stm_to_ssa(m) result(x)
        type(simple_triplet_matrix), intent(in) :: m
        type(simple_sparse_array) :: x
        integer, allocatable :: ii(:,:)
        allocate(ii(m%nnz(),2))
        if(m%nnz()>0) then
            ii(:,1)=m%i; ii(:,2)=m%j
        end if
        x=make_ssa(ii,m%v,[m%nrow,m%ncol])
    end function stm_to_ssa

    function ssa_mean(self) result(m)
        class(simple_sparse_array), intent(in) :: self
        real(dp) :: m
        integer(int64) :: n
        n=product_dims(self%dim)
        if(n==0_int64) then
            m=ieee_value(0.0_dp,ieee_quiet_nan)
        else
            m=sum(self%v)/real(n,dp)
        end if
    end function ssa_mean

    pure function values_equal(a,b) result(eq)
        real(dp), intent(in) :: a,b
        logical :: eq
        eq=(a==b) .or. (ieee_is_nan(a) .and. ieee_is_nan(b))
    end function values_equal

end module slam_ssa
