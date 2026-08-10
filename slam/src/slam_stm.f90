! SPDX-License-Identifier: GPL-2.0-only
!
! Numeric modern-Fortran translation of the simple_triplet_matrix core from
! R package slam 0.1-56 (Hornik, Meyer, Buchta).
module slam_stm
    use iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
    use slam_kinds, only : dp
    use slam_utils, only : argsort_int64
    implicit none
    private

    type, public :: simple_triplet_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        integer, allocatable :: i(:)
        integer, allocatable :: j(:)
        real(dp), allocatable :: v(:)
    contains
        procedure :: valid => stm_valid
        procedure :: nnz => stm_nnz
        procedure :: to_dense => stm_to_dense
        procedure :: transpose => stm_transpose
        procedure :: row_sums => stm_row_sums
        procedure :: col_sums => stm_col_sums
        procedure :: row_means => stm_row_means
        procedure :: col_means => stm_col_means
        procedure :: row_norms => stm_row_norms
        procedure :: col_norms => stm_col_norms
        procedure :: mean => stm_mean
        procedure :: extract => stm_extract
        procedure :: reshape => stm_reshape
    end type simple_triplet_matrix

    interface operator(+)
        module procedure stm_add
    end interface
    interface operator(-)
        module procedure stm_subtract
        module procedure stm_negate
    end interface
    interface operator(*)
        module procedure stm_hadamard
        module procedure stm_times_scalar
        module procedure scalar_times_stm
    end interface
    interface operator(/)
        module procedure stm_div_scalar
    end interface
    interface operator(**)
        module procedure stm_power
    end interface

    public :: make_stm
    public :: dense_to_stm
    public :: simple_triplet_zero_matrix
    public :: simple_triplet_diag_matrix
    public :: stm_rbind
    public :: stm_cbind
    public :: stm_row_scale
    public :: stm_col_scale
    public :: stm_abs
    public :: stm_sqrt
    public :: stm_tcrossprod
    public :: stm_crossprod
    public :: stm_matprod
    public :: stm_tcrossprod_dense
    public :: stm_matprod_dense
    public :: dense_matprod_stm
    public :: stm_set
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)

contains

    function make_stm(i, j, v, nrow, ncol, check_duplicates) result(x)
        integer, intent(in) :: i(:), j(:)
        real(dp), intent(in) :: v(:)
        integer, intent(in), optional :: nrow, ncol
        logical, intent(in), optional :: check_duplicates
        type(simple_triplet_matrix) :: x
        logical :: check

        if (size(i) /= size(j) .or. size(i) /= size(v)) &
            error stop "make_stm: i, j, and v must have equal lengths"

        if (present(nrow)) then
            x%nrow = nrow
        else if (size(i) > 0) then
            x%nrow = maxval(i)
        else
            x%nrow = 0
        end if

        if (present(ncol)) then
            x%ncol = ncol
        else if (size(j) > 0) then
            x%ncol = maxval(j)
        else
            x%ncol = 0
        end if

        allocate(x%i(size(i)), x%j(size(j)), x%v(size(v)))
        x%i = i
        x%j = j
        x%v = v

        if (x%nrow < 0 .or. x%ncol < 0) error stop "make_stm: negative dimension"
        if (size(i) > 0) then
            if (any(i < 1) .or. any(i > x%nrow)) error stop "make_stm: row index out of bounds"
            if (any(j < 1) .or. any(j > x%ncol)) error stop "make_stm: column index out of bounds"
        end if
        check = .true.
        if (present(check_duplicates)) check = check_duplicates
        if (check) then
            if (has_duplicate_pairs(x)) error stop "make_stm: duplicate (i,j) pairs are not allowed"
        end if
    end function make_stm

    pure function stm_nnz(self) result(n)
        class(simple_triplet_matrix), intent(in) :: self
        integer :: n
        if (allocated(self%v)) then
            n = size(self%v)
        else
            n = 0
        end if
    end function stm_nnz

    function stm_valid(self) result(ok)
        class(simple_triplet_matrix), intent(in) :: self
        logical :: ok

        ok = self%nrow >= 0 .and. self%ncol >= 0
        ok = ok .and. allocated(self%i) .and. allocated(self%j) .and. allocated(self%v)
        if (.not. ok) return
        ok = size(self%i) == size(self%j) .and. size(self%i) == size(self%v)
        if (.not. ok) return
        if (size(self%i) > 0) then
            ok = all(self%i >= 1) .and. all(self%i <= self%nrow) .and. &
                 all(self%j >= 1) .and. all(self%j <= self%ncol)
            if (.not. ok) return
        end if
        ok = .not. has_duplicate_pairs(self)
    end function stm_valid

    function has_duplicate_pairs(x) result(dup)
        class(simple_triplet_matrix), intent(in) :: x
        logical :: dup
        integer(int64), allocatable :: key(:)
        integer, allocatable :: order(:)
        integer :: k, n

        n = x%nnz()
        dup = .false.
        if (n < 2) return
        allocate(key(n))
        key = int(x%i, int64) + int(x%j - 1, int64) * int(x%nrow, int64)
        call argsort_int64(key, order)
        do k = 2, n
            if (key(order(k)) == key(order(k - 1))) then
                dup = .true.
                return
            end if
        end do
    end function has_duplicate_pairs

    function dense_to_stm(a) result(x)
        real(dp), intent(in) :: a(:,:)
        type(simple_triplet_matrix) :: x
        integer :: r, c, k, n

        n = 0
        do c = 1, size(a, 2)
            do r = 1, size(a, 1)
                if (ieee_is_nan(a(r,c)) .or. a(r,c) /= 0.0_dp) n = n + 1
            end do
        end do
        allocate(x%i(n), x%j(n), x%v(n))
        x%nrow = size(a, 1)
        x%ncol = size(a, 2)
        k = 0
        do c = 1, x%ncol
            do r = 1, x%nrow
                if (ieee_is_nan(a(r,c)) .or. a(r,c) /= 0.0_dp) then
                    k = k + 1
                    x%i(k) = r
                    x%j(k) = c
                    x%v(k) = a(r,c)
                end if
            end do
        end do
    end function dense_to_stm

    function stm_to_dense(self) result(a)
        class(simple_triplet_matrix), intent(in) :: self
        real(dp), allocatable :: a(:,:)
        integer :: k

        allocate(a(self%nrow, self%ncol))
        a = 0.0_dp
        do k = 1, self%nnz()
            a(self%i(k), self%j(k)) = self%v(k)
        end do
    end function stm_to_dense

    function simple_triplet_zero_matrix(nrow, ncol) result(x)
        integer, intent(in) :: nrow
        integer, intent(in), optional :: ncol
        type(simple_triplet_matrix) :: x
        integer :: nc
        integer, allocatable :: empty(:)
        real(dp), allocatable :: z(:)

        nc = nrow
        if (present(ncol)) nc = ncol
        allocate(empty(0), z(0))
        x = make_stm(empty, empty, z, nrow, nc)
    end function simple_triplet_zero_matrix

    function simple_triplet_diag_matrix(v, nrow) result(x)
        real(dp), intent(in) :: v(:)
        integer, intent(in), optional :: nrow
        type(simple_triplet_matrix) :: x
        integer :: n, k
        integer, allocatable :: ind(:)
        real(dp), allocatable :: vv(:)

        if (present(nrow)) then
            n = nrow
        else
            n = size(v)
        end if
        if (n < 0) error stop "simple_triplet_diag_matrix: negative dimension"
        if (n > 0 .and. size(v) == 0) error stop "simple_triplet_diag_matrix: empty v cannot be recycled"
        allocate(ind(n), vv(n))
        do k = 1, n
            ind(k) = k
            vv(k) = v(mod(k - 1, size(v)) + 1)
        end do
        x = make_stm(ind, ind, vv, n, n)
    end function simple_triplet_diag_matrix

    function stm_transpose(self) result(y)
        class(simple_triplet_matrix), intent(in) :: self
        type(simple_triplet_matrix) :: y
        y = make_stm(self%j, self%i, self%v, self%ncol, self%nrow)
    end function stm_transpose

    function stm_reshape(self, nrow, ncol) result(y)
        class(simple_triplet_matrix), intent(in) :: self
        integer, intent(in) :: nrow, ncol
        type(simple_triplet_matrix) :: y
        integer(int64), allocatable :: pos(:)
        integer, allocatable :: ii(:), jj(:)

        if (nrow < 0 .or. ncol < 0) error stop "stm_reshape: negative dimension"
        if (int(nrow,int64) * int(ncol,int64) /= &
            int(self%nrow,int64) * int(self%ncol,int64)) &
            error stop "stm_reshape: element count must not change"
        allocate(pos(self%nnz()), ii(self%nnz()), jj(self%nnz()))
        if (self%nnz() > 0) then
            if (nrow == 0) error stop "stm_reshape: invalid nonzero matrix shape"
            pos = int(self%i - 1,int64) + int(self%j - 1,int64) * int(self%nrow,int64)
            ii = int(mod(pos,int(nrow,int64))) + 1
            jj = int(pos / int(nrow,int64)) + 1
        end if
        y = make_stm(ii,jj,self%v,nrow,ncol)
    end function stm_reshape

    function stm_row_sums(self, na_rm) result(s)
        class(simple_triplet_matrix), intent(in) :: self
        logical, intent(in), optional :: na_rm
        real(dp), allocatable :: s(:)
        logical :: rm
        integer :: k

        rm = .false.
        if (present(na_rm)) rm = na_rm
        allocate(s(self%nrow))
        s = 0.0_dp
        do k = 1, self%nnz()
            if (rm .and. ieee_is_nan(self%v(k))) cycle
            s(self%i(k)) = s(self%i(k)) + self%v(k)
        end do
    end function stm_row_sums

    function stm_col_sums(self, na_rm) result(s)
        class(simple_triplet_matrix), intent(in) :: self
        logical, intent(in), optional :: na_rm
        real(dp), allocatable :: s(:)
        logical :: rm
        integer :: k

        rm = .false.
        if (present(na_rm)) rm = na_rm
        allocate(s(self%ncol))
        s = 0.0_dp
        do k = 1, self%nnz()
            if (rm .and. ieee_is_nan(self%v(k))) cycle
            s(self%j(k)) = s(self%j(k)) + self%v(k)
        end do
    end function stm_col_sums

    function stm_row_means(self, na_rm) result(m)
        class(simple_triplet_matrix), intent(in) :: self
        logical, intent(in), optional :: na_rm
        real(dp), allocatable :: m(:)
        integer, allocatable :: denom(:)
        logical :: rm
        integer :: k

        rm = .false.
        if (present(na_rm)) rm = na_rm
        m = self%row_sums(rm)
        allocate(denom(self%nrow))
        denom = self%ncol
        if (rm) then
            do k = 1, self%nnz()
                if (ieee_is_nan(self%v(k))) denom(self%i(k)) = denom(self%i(k)) - 1
            end do
        end if
        do k = 1, self%nrow
            if (denom(k) > 0) then
                m(k) = m(k) / real(denom(k), dp)
            else
                m(k) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
        end do
    end function stm_row_means

    function stm_col_means(self, na_rm) result(m)
        class(simple_triplet_matrix), intent(in) :: self
        logical, intent(in), optional :: na_rm
        real(dp), allocatable :: m(:)
        integer, allocatable :: denom(:)
        logical :: rm
        integer :: k

        rm = .false.
        if (present(na_rm)) rm = na_rm
        m = self%col_sums(rm)
        allocate(denom(self%ncol))
        denom = self%nrow
        if (rm) then
            do k = 1, self%nnz()
                if (ieee_is_nan(self%v(k))) denom(self%j(k)) = denom(self%j(k)) - 1
            end do
        end if
        do k = 1, self%ncol
            if (denom(k) > 0) then
                m(k) = m(k) / real(denom(k), dp)
            else
                m(k) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
        end do
    end function stm_col_means

    function stm_row_norms(self, p) result(nrm)
        class(simple_triplet_matrix), intent(in) :: self
        real(dp), intent(in), optional :: p
        real(dp), allocatable :: nrm(:)
        real(dp) :: q
        integer :: k

        q = 2.0_dp
        if (present(p)) q = p
        allocate(nrm(self%nrow))
        nrm = 0.0_dp
        if (q <= 0.0_dp) error stop "stm_row_norms: p must be positive"
        if (.not. ieee_is_finite(q)) then
            do k = 1, self%nnz()
                if (ieee_is_nan(self%v(k))) then
                    nrm(self%i(k)) = self%v(k)
                else if (.not. ieee_is_nan(nrm(self%i(k)))) then
                    nrm(self%i(k)) = max(nrm(self%i(k)), abs(self%v(k)))
                end if
            end do
        else if (q == 1.0_dp) then
            do k = 1, self%nnz()
                nrm(self%i(k)) = nrm(self%i(k)) + abs(self%v(k))
            end do
        else if (q == 2.0_dp) then
            do k = 1, self%nnz()
                nrm(self%i(k)) = nrm(self%i(k)) + self%v(k) * self%v(k)
            end do
            nrm = sqrt(nrm)
        else
            do k = 1, self%nnz()
                nrm(self%i(k)) = nrm(self%i(k)) + abs(self%v(k)) ** q
            end do
            nrm = nrm ** (1.0_dp / q)
        end if
    end function stm_row_norms

    function stm_col_norms(self, p) result(nrm)
        class(simple_triplet_matrix), intent(in) :: self
        real(dp), intent(in), optional :: p
        real(dp), allocatable :: nrm(:)
        real(dp) :: q
        integer :: k

        q = 2.0_dp
        if (present(p)) q = p
        allocate(nrm(self%ncol))
        nrm = 0.0_dp
        if (q <= 0.0_dp) error stop "stm_col_norms: p must be positive"
        if (.not. ieee_is_finite(q)) then
            do k = 1, self%nnz()
                if (ieee_is_nan(self%v(k))) then
                    nrm(self%j(k)) = self%v(k)
                else if (.not. ieee_is_nan(nrm(self%j(k)))) then
                    nrm(self%j(k)) = max(nrm(self%j(k)), abs(self%v(k)))
                end if
            end do
        else if (q == 1.0_dp) then
            do k = 1, self%nnz()
                nrm(self%j(k)) = nrm(self%j(k)) + abs(self%v(k))
            end do
        else if (q == 2.0_dp) then
            do k = 1, self%nnz()
                nrm(self%j(k)) = nrm(self%j(k)) + self%v(k) * self%v(k)
            end do
            nrm = sqrt(nrm)
        else
            do k = 1, self%nnz()
                nrm(self%j(k)) = nrm(self%j(k)) + abs(self%v(k)) ** q
            end do
            nrm = nrm ** (1.0_dp / q)
        end if
    end function stm_col_norms

    function stm_mean(self) result(m)
        class(simple_triplet_matrix), intent(in) :: self
        real(dp) :: m
        integer(int64) :: n
        n = int(self%nrow, int64) * int(self%ncol, int64)
        if (n == 0_int64) then
            m = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            m = sum(self%v) / real(n, dp)
        end if
    end function stm_mean

    function stm_abs(x) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        type(simple_triplet_matrix) :: y
        y = x
        y%v = abs(y%v)
    end function stm_abs

    function stm_sqrt(x) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        type(simple_triplet_matrix) :: y
        y = x
        y%v = sqrt(y%v)
    end function stm_sqrt

    function stm_negate(x) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        type(simple_triplet_matrix) :: y
        y = x
        y%v = -y%v
    end function stm_negate

    function stm_add(x, y) result(z)
        type(simple_triplet_matrix), intent(in) :: x, y
        type(simple_triplet_matrix) :: z
        z = combine_add_sub(x, y, 1.0_dp)
    end function stm_add

    function stm_subtract(x, y) result(z)
        type(simple_triplet_matrix), intent(in) :: x, y
        type(simple_triplet_matrix) :: z
        z = combine_add_sub(x, y, -1.0_dp)
    end function stm_subtract

    function combine_add_sub(x, y, yscale) result(z)
        type(simple_triplet_matrix), intent(in) :: x, y
        real(dp), intent(in) :: yscale
        type(simple_triplet_matrix) :: z
        integer(int64), allocatable :: key(:)
        integer, allocatable :: order(:), ii(:), jj(:), oi(:), oj(:)
        real(dp), allocatable :: vv(:), ov(:)
        integer :: n, k, m, p
        integer(int64) :: current
        real(dp) :: value

        call require_same_shape(x, y, "addition/subtraction")
        n = x%nnz() + y%nnz()
        allocate(key(n), ii(n), jj(n), vv(n), oi(n), oj(n), ov(n))
        if (x%nnz() > 0) then
            ii(1:x%nnz()) = x%i
            jj(1:x%nnz()) = x%j
            vv(1:x%nnz()) = x%v
        end if
        if (y%nnz() > 0) then
            p = x%nnz()
            ii(p+1:n) = y%i
            jj(p+1:n) = y%j
            vv(p+1:n) = yscale * y%v
        end if
        if (n == 0) then
            z = simple_triplet_zero_matrix(x%nrow, x%ncol)
            return
        end if
        key = int(ii, int64) + int(jj - 1, int64) * int(x%nrow, int64)
        call argsort_int64(key, order)
        m = 0
        k = 1
        do while (k <= n)
            current = key(order(k))
            value = 0.0_dp
            do while (k <= n)
                if (key(order(k)) /= current) exit
                value = value + vv(order(k))
                k = k + 1
            end do
            if (ieee_is_nan(value) .or. value /= 0.0_dp) then
                m = m + 1
                oi(m) = int(mod(current - 1_int64, int(x%nrow, int64))) + 1
                oj(m) = int((current - 1_int64) / int(x%nrow, int64)) + 1
                ov(m) = value
            end if
        end do
        z = make_stm(oi(:m), oj(:m), ov(:m), x%nrow, x%ncol)
    end function combine_add_sub

    function stm_hadamard(x, y) result(z)
        type(simple_triplet_matrix), intent(in) :: x, y
        type(simple_triplet_matrix) :: z
        integer, allocatable :: ii(:), jj(:)
        real(dp), allocatable :: vv(:)
        integer :: a, b, n
        real(dp) :: value

        call require_same_shape(x, y, "elementwise multiplication")
        allocate(ii(min(x%nnz(), y%nnz())), jj(min(x%nnz(), y%nnz())), &
                 vv(min(x%nnz(), y%nnz())))
        n = 0
        do a = 1, x%nnz()
            do b = 1, y%nnz()
                if (x%i(a) == y%i(b) .and. x%j(a) == y%j(b)) then
                    value = x%v(a) * y%v(b)
                    if (ieee_is_nan(value) .or. value /= 0.0_dp) then
                        n = n + 1
                        ii(n) = x%i(a)
                        jj(n) = x%j(a)
                        vv(n) = value
                    end if
                    exit
                end if
            end do
        end do
        ! If either stored operand is non-finite, dense zero*Inf/NaN semantics
        ! can create values outside the sparse intersection.
        if (any_nonfinite(x%v) .or. any_nonfinite(y%v)) then
            z = dense_to_stm(x%to_dense() * y%to_dense())
        else
            z = make_stm(ii(:n), jj(:n), vv(:n), x%nrow, x%ncol)
        end if
    end function stm_hadamard

    function stm_times_scalar(x, a) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        real(dp), intent(in) :: a
        type(simple_triplet_matrix) :: y
        if (.not. ieee_is_finite(a)) then
            y = dense_to_stm(x%to_dense() * a)
        else if (a == 0.0_dp) then
            y = simple_triplet_zero_matrix(x%nrow, x%ncol)
        else
            y = x
            y%v = y%v * a
        end if
    end function stm_times_scalar

    function scalar_times_stm(a, x) result(y)
        real(dp), intent(in) :: a
        type(simple_triplet_matrix), intent(in) :: x
        type(simple_triplet_matrix) :: y
        y = stm_times_scalar(x, a)
    end function scalar_times_stm

    function stm_div_scalar(x, a) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        real(dp), intent(in) :: a
        type(simple_triplet_matrix) :: y
        if (a == 0.0_dp .or. .not. ieee_is_finite(a)) then
            y = dense_to_stm(x%to_dense() / a)
        else
            y = x
            y%v = y%v / a
        end if
    end function stm_div_scalar

    function stm_power(x, a) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        real(dp), intent(in) :: a
        type(simple_triplet_matrix) :: y
        if (.not. ieee_is_finite(a) .or. a <= 0.0_dp) &
            error stop "stm_power: exponent must be finite and positive"
        y = x
        y%v = y%v ** a
    end function stm_power

    function stm_row_scale(x, scale) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        real(dp), intent(in) :: scale(:)
        type(simple_triplet_matrix) :: y
        real(dp), allocatable :: a(:,:)
        integer :: k

        if (size(scale) /= x%nrow) error stop "stm_row_scale: scale length mismatch"
        if (any(.not. ieee_is_finite(scale))) then
            a = x%to_dense()
            do k = 1, x%nrow
                a(k,:) = a(k,:) * scale(k)
            end do
            y = dense_to_stm(a)
        else
            y = x
            do k = 1, y%nnz()
                y%v(k) = y%v(k) * scale(y%i(k))
            end do
            call remove_zeros(y)
        end if
    end function stm_row_scale

    function stm_col_scale(x, scale) result(y)
        type(simple_triplet_matrix), intent(in) :: x
        real(dp), intent(in) :: scale(:)
        type(simple_triplet_matrix) :: y
        real(dp), allocatable :: a(:,:)
        integer :: k

        if (size(scale) /= x%ncol) error stop "stm_col_scale: scale length mismatch"
        if (any(.not. ieee_is_finite(scale))) then
            a = x%to_dense()
            do k = 1, x%ncol
                a(:,k) = a(:,k) * scale(k)
            end do
            y = dense_to_stm(a)
        else
            y = x
            do k = 1, y%nnz()
                y%v(k) = y%v(k) * scale(y%j(k))
            end do
            call remove_zeros(y)
        end if
    end function stm_col_scale

    subroutine remove_zeros(x)
        type(simple_triplet_matrix), intent(inout) :: x
        integer, allocatable :: ii(:), jj(:)
        real(dp), allocatable :: vv(:)
        integer :: k, n

        n = count([(ieee_is_nan(x%v(k)) .or. x%v(k) /= 0.0_dp, k=1,x%nnz())])
        if (n == x%nnz()) return
        allocate(ii(n), jj(n), vv(n))
        n = 0
        do k = 1, x%nnz()
            if (ieee_is_nan(x%v(k)) .or. x%v(k) /= 0.0_dp) then
                n = n + 1
                ii(n) = x%i(k); jj(n) = x%j(k); vv(n) = x%v(k)
            end if
        end do
        call move_alloc(ii, x%i)
        call move_alloc(jj, x%j)
        call move_alloc(vv, x%v)
    end subroutine remove_zeros

    function stm_rbind(x, y) result(z)
        type(simple_triplet_matrix), intent(in) :: x, y
        type(simple_triplet_matrix) :: z
        integer, allocatable :: ii(:), jj(:)
        real(dp), allocatable :: vv(:)
        integer :: nx, ny

        if (x%ncol /= y%ncol) error stop "stm_rbind: numbers of columns must match"
        nx = x%nnz(); ny = y%nnz()
        allocate(ii(nx+ny), jj(nx+ny), vv(nx+ny))
        if (nx > 0) then
            ii(:nx)=x%i; jj(:nx)=x%j; vv(:nx)=x%v
        end if
        if (ny > 0) then
            ii(nx+1:)=y%i + x%nrow; jj(nx+1:)=y%j; vv(nx+1:)=y%v
        end if
        z = make_stm(ii, jj, vv, x%nrow+y%nrow, x%ncol)
    end function stm_rbind

    function stm_cbind(x, y) result(z)
        type(simple_triplet_matrix), intent(in) :: x, y
        type(simple_triplet_matrix) :: z
        integer, allocatable :: ii(:), jj(:)
        real(dp), allocatable :: vv(:)
        integer :: nx, ny

        if (x%nrow /= y%nrow) error stop "stm_cbind: numbers of rows must match"
        nx = x%nnz(); ny = y%nnz()
        allocate(ii(nx+ny), jj(nx+ny), vv(nx+ny))
        if (nx > 0) then
            ii(:nx)=x%i; jj(:nx)=x%j; vv(:nx)=x%v
        end if
        if (ny > 0) then
            ii(nx+1:)=y%i; jj(nx+1:)=y%j + x%ncol; vv(nx+1:)=y%v
        end if
        z = make_stm(ii, jj, vv, x%nrow, x%ncol+y%ncol)
    end function stm_cbind

    function stm_extract(self, rows, cols) result(y)
        class(simple_triplet_matrix), intent(in) :: self
        integer, intent(in), optional :: rows(:), cols(:)
        type(simple_triplet_matrix) :: y
        integer, allocatable :: rmap(:), cmap(:), ii(:), jj(:)
        real(dp), allocatable :: vv(:)
        integer :: nr, nc, k, n

        allocate(rmap(self%nrow), cmap(self%ncol))
        rmap = 0; cmap = 0
        if (present(rows)) then
            if (size(rows) > 0) then
                if (any(rows < 1) .or. any(rows > self%nrow)) error stop "stm_extract: row out of bounds"
                if (has_duplicates_int(rows)) error stop "stm_extract: repeated rows not supported"
            end if
            nr = size(rows)
            do k = 1, nr
                rmap(rows(k)) = k
            end do
        else
            nr = self%nrow
            do k = 1, nr
                rmap(k) = k
            end do
        end if
        if (present(cols)) then
            if (size(cols) > 0) then
                if (any(cols < 1) .or. any(cols > self%ncol)) error stop "stm_extract: column out of bounds"
                if (has_duplicates_int(cols)) error stop "stm_extract: repeated columns not supported"
            end if
            nc = size(cols)
            do k = 1, nc
                cmap(cols(k)) = k
            end do
        else
            nc = self%ncol
            do k = 1, nc
                cmap(k) = k
            end do
        end if

        n = 0
        do k = 1, self%nnz()
            if (rmap(self%i(k)) > 0 .and. cmap(self%j(k)) > 0) n = n + 1
        end do
        allocate(ii(n), jj(n), vv(n))
        n = 0
        do k = 1, self%nnz()
            if (rmap(self%i(k)) > 0 .and. cmap(self%j(k)) > 0) then
                n = n + 1
                ii(n) = rmap(self%i(k)); jj(n) = cmap(self%j(k)); vv(n)=self%v(k)
            end if
        end do
        y = make_stm(ii, jj, vv, nr, nc)
    end function stm_extract

    subroutine stm_set(x, rows, cols, values)
        type(simple_triplet_matrix), intent(inout) :: x
        integer, intent(in) :: rows(:), cols(:)
        real(dp), intent(in) :: values(:)
        integer, allocatable :: ii(:), jj(:)
        real(dp), allocatable :: vv(:)
        integer :: k, p, n, cap
        logical :: found
        real(dp) :: value

        if (size(rows) /= size(cols)) error stop "stm_set: coordinate lengths mismatch"
        if (size(values) == 0 .and. size(rows) > 0) error stop "stm_set: empty replacement"
        if (size(rows) > 0) then
            if (any(rows < 1) .or. any(rows > x%nrow) .or. any(cols < 1) .or. any(cols > x%ncol)) &
                error stop "stm_set: subscript out of bounds"
        end if

        cap = x%nnz() + size(rows)
        allocate(ii(cap), jj(cap), vv(cap))
        n = x%nnz()
        if (n > 0) then
            ii(:n)=x%i; jj(:n)=x%j; vv(:n)=x%v
        end if
        do k = 1, size(rows)
            value = values(mod(k - 1, size(values)) + 1)
            found = .false.
            do p = 1, n
                if (ii(p)==rows(k) .and. jj(p)==cols(k)) then
                    vv(p)=value
                    found=.true.
                    exit
                end if
            end do
            if (.not. found) then
                n=n+1; ii(n)=rows(k); jj(n)=cols(k); vv(n)=value
            end if
        end do
        x = make_stm(ii(:n), jj(:n), vv(:n), x%nrow, x%ncol)
        call remove_zeros(x)
    end subroutine stm_set

    function stm_tcrossprod(x, y) result(r)
        type(simple_triplet_matrix), intent(in) :: x
        type(simple_triplet_matrix), intent(in), optional :: y
        real(dp), allocatable :: r(:,:)
        integer, allocatable :: ox(:), oy(:), sx(:), sy(:)
        integer :: col, a, b, fx, lx, fy, ly
        type(simple_triplet_matrix) :: yy

        if (present(y)) then
            if (x%ncol /= y%ncol) error stop "stm_tcrossprod: numbers of columns do not conform"
            if (any_nonfinite(x%v) .or. any_nonfinite(y%v)) then
                r = matmul(x%to_dense(), transpose(y%to_dense()))
                return
            end if
            yy = y
        else
            if (any_nonfinite(x%v)) then
                r = matmul(x%to_dense(), transpose(x%to_dense()))
                return
            end if
            yy = x
        end if

        allocate(r(x%nrow, yy%nrow))
        r = 0.0_dp
        call column_groups(x, ox, sx)
        call column_groups(yy, oy, sy)
        do col = 1, x%ncol
            fx = sx(col-1) + 1; lx = sx(col)
            fy = sy(col-1) + 1; ly = sy(col)
            do a = fx, lx
                do b = fy, ly
                    r(x%i(ox(a)), yy%i(oy(b))) = r(x%i(ox(a)), yy%i(oy(b))) + &
                        x%v(ox(a)) * yy%v(oy(b))
                end do
            end do
        end do
    end function stm_tcrossprod

    function stm_crossprod(x, y) result(r)
        type(simple_triplet_matrix), intent(in) :: x
        type(simple_triplet_matrix), intent(in), optional :: y
        real(dp), allocatable :: r(:,:)
        type(simple_triplet_matrix) :: tx, ty

        tx = x%transpose()
        if (present(y)) then
            if (x%nrow /= y%nrow) error stop "stm_crossprod: numbers of rows do not conform"
            ty = y%transpose()
            r = stm_tcrossprod(tx, ty)
        else
            r = stm_tcrossprod(tx)
        end if
    end function stm_crossprod

    function stm_matprod(x, y) result(r)
        type(simple_triplet_matrix), intent(in) :: x, y
        real(dp), allocatable :: r(:,:)
        integer :: a, b

        if (x%ncol /= y%nrow) error stop "stm_matprod: nonconformable matrices"
        if (any_nonfinite(x%v) .or. any_nonfinite(y%v)) then
            r = matmul(x%to_dense(), y%to_dense())
            return
        end if
        allocate(r(x%nrow, y%ncol))
        r = 0.0_dp
        do a = 1, x%nnz()
            do b = 1, y%nnz()
                if (x%j(a) == y%i(b)) then
                    r(x%i(a), y%j(b)) = r(x%i(a), y%j(b)) + x%v(a) * y%v(b)
                end if
            end do
        end do
    end function stm_matprod

    function stm_tcrossprod_dense(x, y) result(r)
        type(simple_triplet_matrix), intent(in) :: x
        real(dp), intent(in) :: y(:,:)
        real(dp), allocatable :: r(:,:)
        integer :: k

        if (x%ncol /= size(y,2)) error stop "stm_tcrossprod_dense: nonconformable matrices"
        if (any_nonfinite(x%v) .or. any(.not. ieee_is_finite(y))) then
            r = matmul(x%to_dense(), transpose(y))
            return
        end if
        allocate(r(x%nrow, size(y,1)))
        r = 0.0_dp
        do k = 1, x%nnz()
            r(x%i(k),:) = r(x%i(k),:) + x%v(k) * y(:,x%j(k))
        end do
    end function stm_tcrossprod_dense

    function stm_matprod_dense(x, y) result(r)
        type(simple_triplet_matrix), intent(in) :: x
        real(dp), intent(in) :: y(:,:)
        real(dp), allocatable :: r(:,:)
        integer :: k

        if (x%ncol /= size(y,1)) error stop "stm_matprod_dense: nonconformable matrices"
        if (any_nonfinite(x%v) .or. any(.not. ieee_is_finite(y))) then
            r = matmul(x%to_dense(), y)
            return
        end if
        allocate(r(x%nrow, size(y,2)))
        r = 0.0_dp
        do k = 1, x%nnz()
            r(x%i(k),:) = r(x%i(k),:) + x%v(k) * y(x%j(k),:)
        end do
    end function stm_matprod_dense

    function dense_matprod_stm(x, y) result(r)
        real(dp), intent(in) :: x(:,:)
        type(simple_triplet_matrix), intent(in) :: y
        real(dp), allocatable :: r(:,:)
        integer :: k

        if (size(x,2) /= y%nrow) error stop "dense_matprod_stm: nonconformable matrices"
        if (any(.not. ieee_is_finite(x)) .or. any_nonfinite(y%v)) then
            r = matmul(x, y%to_dense())
            return
        end if
        allocate(r(size(x,1), y%ncol))
        r = 0.0_dp
        do k = 1, y%nnz()
            r(:,y%j(k)) = r(:,y%j(k)) + x(:,y%i(k)) * y%v(k)
        end do
    end function dense_matprod_stm

    subroutine column_groups(x, order, cumulative)
        type(simple_triplet_matrix), intent(in) :: x
        integer, allocatable, intent(out) :: order(:)
        integer, allocatable, intent(out) :: cumulative(:)
        integer(int64), allocatable :: key(:)
        integer :: k

        allocate(key(x%nnz()), cumulative(0:x%ncol))
        key = int(x%j, int64)
        call argsort_int64(key, order)
        cumulative = 0
        do k = 1, x%nnz()
            cumulative(x%j(k)) = cumulative(x%j(k)) + 1
        end do
        do k = 1, x%ncol
            cumulative(k) = cumulative(k) + cumulative(k-1)
        end do
    end subroutine column_groups

    pure function any_nonfinite(v) result(ans)
        real(dp), intent(in) :: v(:)
        logical :: ans
        ans = any(.not. ieee_is_finite(v))
    end function any_nonfinite

    subroutine require_same_shape(x, y, where)
        type(simple_triplet_matrix), intent(in) :: x, y
        character(len=*), intent(in) :: where
        if (x%nrow /= y%nrow .or. x%ncol /= y%ncol) &
            error stop "slam_stm: nonconformable matrices in " // where
    end subroutine require_same_shape

    pure function has_duplicates_int(a) result(dup)
        integer, intent(in) :: a(:)
        logical :: dup
        integer :: i, j
        dup = .false.
        do i = 1, size(a)-1
            do j = i+1, size(a)
                if (a(i) == a(j)) then
                    dup = .true.
                    return
                end if
            end do
        end do
    end function has_duplicates_int

end module slam_stm
