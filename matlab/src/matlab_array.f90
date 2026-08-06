! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module matlab_array
    use matlab_kinds, only : dp
    use matlab_types, only : meshgrid2d_result, meshgrid3d_result
    implicit none
    private

    public :: zeros
    public :: ones
    public :: eye
    public :: find
    public :: fliplr
    public :: flipud
    public :: rot90
    public :: repmat
    public :: repmat_vector
    public :: reshape2d
    public :: meshgrid
    public :: meshgrid3
    public :: padarray
    public :: shape_of
    public :: size_dim
    public :: numel
    public :: ndims
    public :: isempty

    interface find
        module procedure find_real_vector
        module procedure find_real_matrix
        module procedure find_logical_vector
        module procedure find_logical_matrix
        module procedure find_integer_vector
        module procedure find_integer_matrix
    end interface find

    interface fliplr
        module procedure fliplr_vector
        module procedure fliplr_matrix
    end interface fliplr

    interface flipud
        module procedure flipud_vector
        module procedure flipud_matrix
    end interface flipud

    interface repmat
        module procedure repmat_matrix
        module procedure repmat_vector
    end interface repmat

    interface reshape2d
        module procedure reshape2d_matrix
        module procedure reshape2d_vector
    end interface reshape2d

contains

    function zeros(m, n) result(a)
        integer, intent(in) :: m
        integer, intent(in), optional :: n
        real(dp), allocatable :: a(:, :)
        integer :: nn

        nn = m
        if (present(n)) nn = n
        allocate(a(max(m, 0), max(nn, 0)))
        a = 0.0_dp
    end function zeros

    function ones(m, n) result(a)
        integer, intent(in) :: m
        integer, intent(in), optional :: n
        real(dp), allocatable :: a(:, :)
        integer :: nn

        nn = m
        if (present(n)) nn = n
        allocate(a(max(m, 0), max(nn, 0)))
        a = 1.0_dp
    end function ones

    function eye(m, n) result(a)
        integer, intent(in) :: m
        integer, intent(in), optional :: n
        real(dp), allocatable :: a(:, :)
        integer :: i, nn

        nn = m
        if (present(n)) nn = n
        allocate(a(max(m, 0), max(nn, 0)))
        a = 0.0_dp
        do i = 1, min(m, nn)
            a(i, i) = 1.0_dp
        end do
    end function eye

    function find_real_vector(x) result(indices)
        real(dp), intent(in) :: x(:)
        integer, allocatable :: indices(:)
        integer :: i, k

        allocate(indices(count(x /= 0.0_dp)))
        k = 0
        do i = 1, size(x)
            if (x(i) /= 0.0_dp) then
                k = k + 1
                indices(k) = i
            end if
        end do
    end function find_real_vector

    function find_logical_vector(x) result(indices)
        logical, intent(in) :: x(:)
        integer, allocatable :: indices(:)
        integer :: i, k

        allocate(indices(count(x)))
        k = 0
        do i = 1, size(x)
            if (x(i)) then
                k = k + 1
                indices(k) = i
            end if
        end do
    end function find_logical_vector


    function find_real_matrix(x) result(indices)
        real(dp), intent(in) :: x(:, :)
        integer, allocatable :: indices(:)
        integer :: i, j, k, linear

        allocate(indices(count(x /= 0.0_dp)))
        k = 0
        linear = 0
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                linear = linear + 1
                if (x(i, j) /= 0.0_dp) then
                    k = k + 1
                    indices(k) = linear
                end if
            end do
        end do
    end function find_real_matrix

    function find_logical_matrix(x) result(indices)
        logical, intent(in) :: x(:, :)
        integer, allocatable :: indices(:)
        integer :: i, j, k, linear

        allocate(indices(count(x)))
        k = 0
        linear = 0
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                linear = linear + 1
                if (x(i, j)) then
                    k = k + 1
                    indices(k) = linear
                end if
            end do
        end do
    end function find_logical_matrix

    function find_integer_vector(x) result(indices)
        integer, intent(in) :: x(:)
        integer, allocatable :: indices(:)
        integer :: i, k

        allocate(indices(count(x /= 0)))
        k = 0
        do i = 1, size(x)
            if (x(i) /= 0) then
                k = k + 1
                indices(k) = i
            end if
        end do
    end function find_integer_vector

    function find_integer_matrix(x) result(indices)
        integer, intent(in) :: x(:, :)
        integer, allocatable :: indices(:)
        integer :: i, j, k, linear

        allocate(indices(count(x /= 0)))
        k = 0
        linear = 0
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                linear = linear + 1
                if (x(i, j) /= 0) then
                    k = k + 1
                    indices(k) = linear
                end if
            end do
        end do
    end function find_integer_matrix

    function fliplr_matrix(a) result(b)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable :: b(:, :)
        integer :: j

        allocate(b(size(a, 1), size(a, 2)))
        do j = 1, size(a, 2)
            b(:, j) = a(:, size(a, 2) - j + 1)
        end do
    end function fliplr_matrix

    function flipud_matrix(a) result(b)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable :: b(:, :)
        integer :: i

        allocate(b(size(a, 1), size(a, 2)))
        do i = 1, size(a, 1)
            b(i, :) = a(size(a, 1) - i + 1, :)
        end do
    end function flipud_matrix


    function fliplr_vector(a) result(b)
        real(dp), intent(in) :: a(:)
        real(dp), allocatable :: b(:)
        integer :: i

        allocate(b(size(a)))
        do i = 1, size(a)
            b(i) = a(size(a) - i + 1)
        end do
    end function fliplr_vector

    function flipud_vector(a) result(b)
        real(dp), intent(in) :: a(:)
        real(dp), allocatable :: b(:)

        b = fliplr_vector(a)
    end function flipud_vector

    function rot90(a, k) result(b)
        real(dp), intent(in) :: a(:, :)
        integer, intent(in), optional :: k
        real(dp), allocatable :: b(:, :)
        integer :: i, j, kk

        kk = 1
        if (present(k)) kk = k
        kk = modulo(kk, 4)

        select case (kk)
        case (0)
            allocate(b(size(a, 1), size(a, 2)))
            b = a
        case (1)
            allocate(b(size(a, 2), size(a, 1)))
            do i = 1, size(a, 1)
                do j = 1, size(a, 2)
                    b(size(a, 2) - j + 1, i) = a(i, j)
                end do
            end do
        case (2)
            allocate(b(size(a, 1), size(a, 2)))
            do i = 1, size(a, 1)
                do j = 1, size(a, 2)
                    b(size(a, 1) - i + 1, size(a, 2) - j + 1) = a(i, j)
                end do
            end do
        case (3)
            allocate(b(size(a, 2), size(a, 1)))
            do i = 1, size(a, 1)
                do j = 1, size(a, 2)
                    b(j, size(a, 1) - i + 1) = a(i, j)
                end do
            end do
        end select
    end function rot90

    function repmat_matrix(a, m, n, source_compatible) result(b)
        real(dp), intent(in) :: a(:, :)
        integer, intent(in) :: m
        integer, intent(in), optional :: n
        logical, intent(in), optional :: source_compatible
        real(dp), allocatable :: b(:, :)
        real(dp), allocatable :: standard(:, :)
        integer :: bi, bj, nn
        logical :: source_mode

        nn = m
        if (present(n)) nn = n
        source_mode = .true.
        if (present(source_compatible)) source_mode = source_compatible

        if (source_mode .and. nn == 1 .and. .not. (m == 1 .and. nn == 1)) then
            allocate(standard(size(a, 1), size(a, 2) * m))
            do bj = 0, m - 1
                standard(:, bj * size(a, 2) + 1:(bj + 1) * size(a, 2)) = a
            end do
            allocate(b(size(standard, 2), size(standard, 1)))
            b = transpose(standard)
            return
        end if

        allocate(b(size(a, 1) * m, size(a, 2) * nn))
        do bi = 0, m - 1
            do bj = 0, nn - 1
                b(bi * size(a, 1) + 1:(bi + 1) * size(a, 1), &
                  bj * size(a, 2) + 1:(bj + 1) * size(a, 2)) = a
            end do
        end do
    end function repmat_matrix

    function repmat_vector(a, m, n, source_compatible) result(b)
        real(dp), intent(in) :: a(:)
        integer, intent(in) :: m
        integer, intent(in), optional :: n
        logical, intent(in), optional :: source_compatible
        real(dp), allocatable :: b(:, :)
        real(dp), allocatable :: column(:, :)
        integer :: nn

        nn = m
        if (present(n)) nn = n
        allocate(column(size(a), 1))
        column(:, 1) = a
        if (present(source_compatible)) then
            b = repmat_matrix(column, m, nn, source_compatible)
        else
            b = repmat_matrix(column, m, nn)
        end if
    end function repmat_vector

    function reshape2d_matrix(a, m, n) result(b)
        real(dp), intent(in) :: a(:, :)
        integer, intent(in) :: m
        integer, intent(in) :: n
        real(dp), allocatable :: b(:, :)

        if (m * n /= size(a)) then
            allocate(b(0, 0))
            return
        end if
        allocate(b(m, n))
        b = reshape(a, [m, n])
    end function reshape2d_matrix


    function reshape2d_vector(a, m, n) result(b)
        real(dp), intent(in) :: a(:)
        integer, intent(in) :: m
        integer, intent(in) :: n
        real(dp), allocatable :: b(:, :)

        if (m * n /= size(a)) then
            allocate(b(0, 0))
            return
        end if
        allocate(b(m, n))
        b = reshape(a, [m, n])
    end function reshape2d_vector

    function meshgrid(x, y) result(grid)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: y(:)
        type(meshgrid2d_result) :: grid
        real(dp), allocatable :: yy(:)
        integer :: i, j

        if (present(y)) then
            allocate(yy(size(y)))
            yy = y
        else
            allocate(yy(size(x)))
            yy = x
        end if
        allocate(grid%x(size(yy), size(x)), grid%y(size(yy), size(x)))
        do i = 1, size(yy)
            do j = 1, size(x)
                grid%x(i, j) = x(j)
                grid%y(i, j) = yy(i)
            end do
        end do
    end function meshgrid

    function meshgrid3(x, y, z) result(grid)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: y(:)
        real(dp), intent(in), optional :: z(:)
        type(meshgrid3d_result) :: grid
        real(dp), allocatable :: yy(:), zz(:)
        integer :: i, j, k

        if (present(y)) then
            allocate(yy(size(y)))
            yy = y
        else
            allocate(yy(size(x)))
            yy = x
        end if
        if (present(z)) then
            allocate(zz(size(z)))
            zz = z
        else
            allocate(zz(size(x)))
            zz = x
        end if

        allocate(grid%x(size(yy), size(x), size(zz)))
        allocate(grid%y(size(yy), size(x), size(zz)))
        allocate(grid%z(size(yy), size(x), size(zz)))
        do k = 1, size(zz)
            do j = 1, size(x)
                do i = 1, size(yy)
                    grid%x(i, j, k) = x(j)
                    grid%y(i, j, k) = yy(i)
                    grid%z(i, j, k) = zz(k)
                end do
            end do
        end do
    end function meshgrid3

    function padarray(a, pad_rows, pad_cols, method, direction, pad_value) result(b)
        real(dp), intent(in) :: a(:, :)
        integer, intent(in) :: pad_rows
        integer, intent(in), optional :: pad_cols
        character(len=*), intent(in), optional :: method
        character(len=*), intent(in), optional :: direction
        real(dp), intent(in), optional :: pad_value
        real(dp), allocatable :: b(:, :)
        character(len=16) :: meth, dir
        integer :: pc, nr, nc, i, j, raw_i, raw_j, ii, jj, shift_r, shift_c
        real(dp) :: value

        pc = pad_rows
        if (present(pad_cols)) pc = pad_cols
        meth = 'constant'
        if (present(method)) meth = lower_string(method)
        dir = 'both'
        if (present(direction)) dir = lower_string(direction)
        value = 0.0_dp
        if (present(pad_value)) value = pad_value

        select case (trim(dir))
        case ('pre', 'post')
            nr = size(a, 1) + pad_rows
            nc = size(a, 2) + pc
        case default
            nr = size(a, 1) + 2 * pad_rows
            nc = size(a, 2) + 2 * pc
            dir = 'both'
        end select
        allocate(b(nr, nc))

        if (trim(meth) == 'constant') then
            b = value
            shift_r = merge(pad_rows, 0, trim(dir) /= 'post')
            shift_c = merge(pc, 0, trim(dir) /= 'post')
            b(shift_r + 1:shift_r + size(a, 1), &
              shift_c + 1:shift_c + size(a, 2)) = a
            return
        end if

        shift_r = merge(pad_rows, 0, trim(dir) /= 'post')
        shift_c = merge(pc, 0, trim(dir) /= 'post')
        do j = 1, nc
            raw_j = j - shift_c
            jj = map_index(raw_j, size(a, 2), trim(meth))
            do i = 1, nr
                raw_i = i - shift_r
                ii = map_index(raw_i, size(a, 1), trim(meth))
                b(i, j) = a(ii, jj)
            end do
        end do
    end function padarray

    function shape_of(x) result(dims)
        class(*), intent(in) :: x(..)
        integer, allocatable :: dims(:)

        select rank (x)
        rank (0)
            allocate(dims(2))
            dims = [1, 1]
        rank (1)
            allocate(dims(2))
            dims = [1, size(x)]
        rank (2)
            allocate(dims(2))
            dims = shape(x)
        rank (3)
            allocate(dims(3))
            dims = shape(x)
        rank default
            allocate(dims(1))
            dims = [size(x)]
        end select
    end function shape_of

    function size_dim(x, dimen) result(n)
        class(*), intent(in) :: x(..)
        integer, intent(in) :: dimen
        integer :: n
        integer, allocatable :: dims(:)

        dims = shape_of(x)
        if (dimen < 1) then
            n = 0
        else if (dimen <= size(dims)) then
            n = dims(dimen)
        else
            n = 1
        end if
    end function size_dim

    function numel(x) result(n)
        class(*), intent(in) :: x(..)
        integer :: n

        n = size(x)
    end function numel

    function ndims(x) result(n)
        class(*), intent(in) :: x(..)
        integer :: n
        integer, allocatable :: dims(:)

        dims = shape_of(x)
        n = size(dims)
    end function ndims

    function isempty(x) result(empty)
        class(*), intent(in) :: x(..)
        logical :: empty

        empty = size(x) == 0
    end function isempty

    pure function map_index(raw, n, method) result(idx)
        integer, intent(in) :: raw
        integer, intent(in) :: n
        character(len=*), intent(in) :: method
        integer :: idx, t

        select case (trim(method))
        case ('circular')
            idx = modulo(raw - 1, n) + 1
        case ('replicate')
            idx = min(max(raw, 1), n)
        case ('symmetric')
            t = modulo(raw - 1, 2 * n) + 1
            if (t <= n) then
                idx = t
            else
                idx = 2 * n - t + 1
            end if
        case default
            idx = min(max(raw, 1), n)
        end select
    end function map_index

    pure function lower_string(text) result(out)
        character(len=*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i, code

        out = text
        do i = 1, len(text)
            code = iachar(out(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                out(i:i) = achar(code + 32)
            end if
        end do
    end function lower_string
end module matlab_array
