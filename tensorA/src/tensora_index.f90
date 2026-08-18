! SPDX-License-Identifier: GPL-2.0-or-later
module tensora_index
    use tensora_kinds, only : axis_name_len
    use tensora_types, only : tensor_t
    implicit none
    private

    public :: strides_of, linear_index, decode_index
    public :: complement, full_order, positions_by_name, names_to_positions
    public :: same_axis_set, contraname, is_covariate_name, is_contravariate_name
    public :: as_covariate_name, as_contravariate_name
    public :: marked_name, product_shape

contains

    pure function strides_of(shape) result(stride)
        integer, intent(in) :: shape(:)
        integer, allocatable :: stride(:)
        integer :: k

        allocate(stride(size(shape)))
        if (size(shape) == 0) return
        stride(1) = 1
        do k = 2, size(shape)
            stride(k) = stride(k - 1) * shape(k - 1)
        end do
    end function strides_of

    pure integer function linear_index(sub, shape) result(idx)
        integer, intent(in) :: sub(:), shape(:)
        integer :: k, stride

        idx = 1
        stride = 1
        do k = 1, size(shape)
            idx = idx + (sub(k) - 1) * stride
            stride = stride * shape(k)
        end do
    end function linear_index

    pure subroutine decode_index(idx, shape, sub)
        integer, intent(in) :: idx
        integer, intent(in) :: shape(:)
        integer, intent(out) :: sub(:)
        integer :: q, k

        q = idx - 1
        do k = 1, size(shape)
            sub(k) = mod(q, shape(k)) + 1
            q = q / shape(k)
        end do
    end subroutine decode_index

    pure function complement(n, used) result(rest)
        integer, intent(in) :: n
        integer, intent(in) :: used(:)
        integer, allocatable :: rest(:)
        logical, allocatable :: take(:)
        integer :: k, m

        allocate(take(n))
        take = .true.
        do k = 1, size(used)
            if (used(k) >= 1 .and. used(k) <= n) take(used(k)) = .false.
        end do
        m = count(take)
        allocate(rest(m))
        m = 0
        do k = 1, n
            if (take(k)) then
                m = m + 1
                rest(m) = k
            end if
        end do
    end function complement

    pure function full_order(n, first, last) result(order)
        integer, intent(in) :: n
        integer, intent(in), optional :: first(:), last(:)
        integer, allocatable :: order(:)
        logical, allocatable :: used(:)
        integer :: k, m, v

        allocate(order(n), used(n))
        used = .false.
        m = 0
        if (present(first)) then
            do k = 1, size(first)
                v = first(k)
                if (v < 1 .or. v > n) cycle
                if (.not. used(v)) then
                    m = m + 1
                    order(m) = v
                    used(v) = .true.
                end if
            end do
        end if
        do v = 1, n
            if (present(last)) then
                if (any(last == v)) cycle
            end if
            if (.not. used(v)) then
                m = m + 1
                order(m) = v
                used(v) = .true.
            end if
        end do
        if (present(last)) then
            do k = 1, size(last)
                v = last(k)
                if (v < 1 .or. v > n) cycle
                if (.not. used(v)) then
                    m = m + 1
                    order(m) = v
                    used(v) = .true.
                end if
            end do
        end if
        if (m /= n) error stop "full_order: invalid permutation"
    end function full_order

    function positions_by_name(t, names, missing_ok, both) result(pos)
        type(tensor_t), intent(in) :: t
        character(len=*), intent(in) :: names(:)
        logical, intent(in), optional :: missing_ok, both
        integer, allocatable :: pos(:)
        logical :: miss, covboth
        integer :: k, j, found
        character(len=axis_name_len) :: target, candidate

        miss = .false.
        if (present(missing_ok)) miss = missing_ok
        covboth = .false.
        if (present(both)) covboth = both
        allocate(pos(size(names)))
        do k = 1, size(names)
            target = trim(names(k))
            found = 0
            do j = 1, t%rank()
                candidate = trim(t%axis(j))
                if (covboth) then
                    if (trim(as_covariate_name(target)) == trim(as_covariate_name(candidate))) then
                        if (found /= 0) error stop "positions_by_name: ambiguous axis"
                        found = j
                    end if
                else if (trim(target) == trim(candidate)) then
                    if (found /= 0) error stop "positions_by_name: ambiguous axis"
                    found = j
                end if
            end do
            if (found == 0 .and. .not. miss) error stop "positions_by_name: axis not found"
            pos(k) = found
        end do
    end function positions_by_name

    function names_to_positions(axis, names, missing_ok, both) result(pos)
        character(len=*), intent(in) :: axis(:), names(:)
        logical, intent(in), optional :: missing_ok, both
        type(tensor_t) :: tmp
        integer, allocatable :: pos(:)

        allocate(tmp%shape(size(axis)), tmp%axis(size(axis)), tmp%data(1))
        tmp%shape = 1
        tmp%axis = axis
        tmp%data = (0.0, 0.0)
        pos = positions_by_name(tmp, names, missing_ok, both)
    end function names_to_positions

    pure logical function same_axis_set(a, b) result(ok)
        character(len=*), intent(in) :: a(:), b(:)
        integer :: k

        ok = size(a) == size(b)
        if (.not. ok) return
        do k = 1, size(a)
            if (.not. any(trim(a(k)) == b)) then
                ok = .false.
                return
            end if
        end do
    end function same_axis_set

    pure function contraname(name) result(out)
        character(len=*), intent(in) :: name
        character(len=axis_name_len) :: out
        integer :: n

        n = len_trim(name)
        if (n > 0 .and. name(1:1) == '^') then
            out = adjustl(name(2:n))
        else
            out = '^' // trim(name)
        end if
    end function contraname

    pure logical function is_covariate_name(name) result(ans)
        character(len=*), intent(in) :: name

        ans = len_trim(name) == 0 .or. name(1:1) /= '^'
    end function is_covariate_name

    pure logical function is_contravariate_name(name) result(ans)
        character(len=*), intent(in) :: name

        ans = len_trim(name) > 0 .and. name(1:1) == '^'
    end function is_contravariate_name


    pure function as_contravariate_name(name) result(out)
        character(len=*), intent(in) :: name
        character(len=axis_name_len) :: out

        if (is_contravariate_name(name)) then
            out = trim(name)
        else
            out = '^' // trim(name)
        end if
    end function as_contravariate_name

    pure function marked_name(name, mark) result(out)
        character(len=*), intent(in) :: name, mark
        character(len=axis_name_len) :: out

        out = trim(name) // trim(mark)
    end function marked_name

    pure function as_covariate_name(name) result(out)
        character(len=*), intent(in) :: name
        character(len=axis_name_len) :: out
        integer :: n

        n = len_trim(name)
        if (n > 0 .and. name(1:1) == '^') then
            out = adjustl(name(2:n))
        else
            out = trim(name)
        end if
    end function as_covariate_name

    pure integer function product_shape(shape) result(p)
        integer, intent(in) :: shape(:)

        if (size(shape) == 0) then
            p = 1
        else
            p = product(shape)
        end if
    end function product_shape

end module tensora_index
