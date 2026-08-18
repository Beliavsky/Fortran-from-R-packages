! SPDX-License-Identifier: GPL-2.0-or-later
module tensora_core
    use tensora_kinds, only : dp, axis_name_len
    use tensora_types, only : tensor_t, tensor, tensor_zeros, tensor_ones, scalar_tensor
    use tensora_index, only : linear_index, decode_index, complement, full_order, &
        positions_by_name, contraname, marked_name, product_shape
    implicit none
    private

    interface reorder_tensor
        module procedure reorder_tensor_pos
        module procedure reorder_tensor_names
    end interface reorder_tensor

    interface mul_tensor
        module procedure mul_tensor_pos
        module procedure mul_tensor_names
    end interface mul_tensor

    public :: reorder_tensor, reorder_tensor_pos, reorder_tensor_names, pos_tensor
    public :: rename_axis, rename_first_axes, mark_tensor
    public :: mul_tensor, mul_tensor_pos, mul_tensor_names
    public :: trace_tensor, margin_tensor, diagmul_tensor
    public :: delta_tensor, diag_tensor, tripledelta_tensor, one_tensor
    public :: add_tensor, sub_tensor, elem_mul_tensor, elem_div_tensor
    public :: scale_tensor, repeat_tensor, slice_tensor, undrop_tensor, untensor_tensor
    public :: bind_tensor, einstein_pair, riemann_pair
    public :: contraname_tensor, is_covariate_tensor, is_contravariate_tensor

contains

    function pos_tensor(shape) result(pos)
        integer, intent(in) :: shape(:)
        integer, allocatable :: pos(:,:)
        integer :: k

        allocate(pos(product_shape(shape), size(shape)))
        do k = 1, product_shape(shape)
            call decode_index(k, shape, pos(k,:))
        end do
    end function pos_tensor

    function reorder_tensor_pos(x, first) result(y)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: first(:)
        type(tensor_t) :: y
        integer, allocatable :: order(:), osub(:), isub(:)
        integer :: k, idx

        order = full_order(x%rank(), first=first)
        allocate(y%shape(x%rank()), y%axis(x%rank()), y%data(x%nelem()))
        y%shape = x%shape(order)
        y%axis = x%axis(order)
        allocate(osub(y%rank()), isub(x%rank()))
        do k = 1, y%nelem()
            call decode_index(k, y%shape, osub)
            isub(order) = osub
            idx = linear_index(isub, x%shape)
            y%data(k) = x%data(idx)
        end do
    end function reorder_tensor_pos

    function reorder_tensor_names(x, first) result(y)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: first(:)
        type(tensor_t) :: y
        integer, allocatable :: pos(:)

        pos = positions_by_name(x, first)
        y = reorder_tensor_pos(x, pos)
    end function reorder_tensor_names

    function rename_first_axes(x, new_names) result(y)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: new_names(:)
        type(tensor_t) :: y
        integer :: k

        if (size(new_names) > x%rank()) error stop "rename_first_axes: too many names"
        y = x
        do k = 1, size(new_names)
            y%axis(k) = trim(new_names(k))
        end do
    end function rename_first_axes

    function rename_axis(x, old_name, new_name) result(y)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: old_name, new_name
        type(tensor_t) :: y
        integer :: p, k

        y = x
        p = x%axis_pos(old_name)
        if (p <= 0) error stop "rename_axis: axis not found or ambiguous"
        do k = 1, x%rank()
            if (k /= p .and. trim(x%axis(k)) == trim(new_name)) then
                error stop "rename_axis: duplicate target axis"
            end if
        end do
        y%axis(p) = trim(new_name)
    end function rename_axis

    function mark_tensor(x, mark, axes) result(y)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in), optional :: mark
        character(len=*), intent(in), optional :: axes(:)
        type(tensor_t) :: y
        character(len=axis_name_len) :: mk
        integer, allocatable :: pos(:)
        integer :: k

        mk = "'"
        if (present(mark)) mk = trim(mark)
        y = x
        if (present(axes)) then
            pos = positions_by_name(x, axes)
        else
            allocate(pos(x%rank()))
            pos = [(k, k=1, x%rank())]
        end if
        do k = 1, size(pos)
            y%axis(pos(k)) = marked_name(y%axis(pos(k)), mk)
        end do
    end function mark_tensor

    function mul_tensor_pos(x, ix, y, iy, byx, byy) result(z)
        type(tensor_t), intent(in) :: x, y
        integer, intent(in) :: ix(:), iy(:)
        integer, intent(in), optional :: byx(:), byy(:)
        type(tensor_t) :: z
        integer, allocatable :: bx(:), by(:), rx(:), ry(:), ox(:), oy(:)
        type(tensor_t) :: xr, yr
        complex(dp), allocatable :: xm(:,:,:), ym(:,:,:), zm(:,:,:)
        integer :: nxr, nyr, nb, inner, outerx, outery, para, q

        if (size(ix) /= size(iy)) error stop "mul_tensor: contraction rank mismatch"
        if (size(ix) > 0) then
            if (any(x%shape(ix) /= y%shape(iy))) error stop "mul_tensor: incompatible contraction dimensions"
        end if
        if (present(byx) .neqv. present(byy)) error stop "mul_tensor: byx and byy must both be present"
        if (present(byx)) then
            if (size(byx) /= size(byy)) error stop "mul_tensor: by rank mismatch"
            allocate(bx(size(byx)), by(size(byy)))
            bx = byx
            by = byy
            if (size(bx) > 0) then
                if (any(x%shape(bx) /= y%shape(by))) error stop "mul_tensor: incompatible by dimensions"
            end if
        else
            allocate(bx(0), by(0))
        end if
        if (has_overlap(ix, bx) .or. has_overlap(iy, by)) error stop "mul_tensor: contracted and by axes overlap"

        rx = complement(x%rank(), [ix, bx])
        ry = complement(y%rank(), [iy, by])
        ox = [rx, ix, bx]
        oy = [iy, ry, by]
        xr = reorder_tensor_pos(x, ox)
        yr = reorder_tensor_pos(y, oy)

        nxr = size(rx)
        nyr = size(ry)
        nb = size(bx)
        outerx = product_shape(x%shape(rx))
        outery = product_shape(y%shape(ry))
        inner = product_shape(x%shape(ix))
        para = product_shape(x%shape(bx))
        allocate(xm(outerx,inner,para), ym(inner,outery,para), zm(outerx,outery,para))
        xm = reshape(xr%data, shape(xm))
        ym = reshape(yr%data, shape(ym))
        do q = 1, para
            zm(:,:,q) = matmul(xm(:,:,q), ym(:,:,q))
        end do

        allocate(z%shape(nxr + nyr + nb), z%axis(nxr + nyr + nb))
        if (nxr > 0) then
            z%shape(1:nxr) = x%shape(rx)
            z%axis(1:nxr) = x%axis(rx)
        end if
        if (nyr > 0) then
            z%shape(nxr+1:nxr+nyr) = y%shape(ry)
            z%axis(nxr+1:nxr+nyr) = y%axis(ry)
        end if
        if (nb > 0) then
            z%shape(nxr+nyr+1:) = x%shape(bx)
            z%axis(nxr+nyr+1:) = x%axis(bx)
        end if
        allocate(z%data(size(zm)))
        z%data = reshape(zm, [size(zm)])
    end function mul_tensor_pos

    function mul_tensor_names(x, inames, y, jnames, by_names) result(z)
        type(tensor_t), intent(in) :: x, y
        character(len=*), intent(in) :: inames(:), jnames(:)
        character(len=*), intent(in), optional :: by_names(:)
        type(tensor_t) :: z
        integer, allocatable :: ix(:), iy(:), bx(:), by(:)

        ix = positions_by_name(x, inames)
        iy = positions_by_name(y, jnames)
        if (present(by_names)) then
            bx = positions_by_name(x, by_names)
            by = positions_by_name(y, by_names)
            z = mul_tensor_pos(x, ix, y, iy, bx, by)
        else
            z = mul_tensor_pos(x, ix, y, iy)
        end if
    end function mul_tensor_names

    function trace_tensor(x, ia, ja) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: ia(:), ja(:)
        type(tensor_t) :: z
        integer, allocatable :: rest(:), zsub(:), xsub(:), csub(:)
        integer :: zi, ci, p, inner

        if (size(ia) /= size(ja)) error stop "trace_tensor: pair-count mismatch"
        if (size(ia) == 0) then
            z = x
            return
        end if
        if (has_overlap(ia, ja)) error stop "trace_tensor: axes must be distinct"
        if (any(x%shape(ia) /= x%shape(ja))) error stop "trace_tensor: incompatible dimensions"
        rest = complement(x%rank(), [ia, ja])
        allocate(z%shape(size(rest)), z%axis(size(rest)))
        if (size(rest) > 0) then
            z%shape = x%shape(rest)
            z%axis = x%axis(rest)
        end if
        allocate(z%data(product_shape(z%shape)))
        z%data = (0.0_dp, 0.0_dp)
        allocate(zsub(z%rank()), xsub(x%rank()), csub(size(ia)))
        inner = product_shape(x%shape(ia))
        do zi = 1, z%nelem()
            if (z%rank() > 0) call decode_index(zi, z%shape, zsub)
            xsub = 1
            do p = 1, size(rest)
                xsub(rest(p)) = zsub(p)
            end do
            do ci = 1, inner
                call decode_index(ci, x%shape(ia), csub)
                do p = 1, size(ia)
                    xsub(ia(p)) = csub(p)
                    xsub(ja(p)) = csub(p)
                end do
                z%data(zi) = z%data(zi) + x%data(linear_index(xsub, x%shape))
            end do
        end do
    end function trace_tensor

    function margin_tensor(x, remove) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: remove(:)
        type(tensor_t) :: z
        integer, allocatable :: rest(:), zsub(:), xsub(:), rsub(:)
        integer :: zi, ri, p, nsum

        rest = complement(x%rank(), remove)
        allocate(z%shape(size(rest)), z%axis(size(rest)))
        if (size(rest) > 0) then
            z%shape = x%shape(rest)
            z%axis = x%axis(rest)
        end if
        allocate(z%data(product_shape(z%shape)))
        z%data = (0.0_dp, 0.0_dp)
        allocate(zsub(z%rank()), xsub(x%rank()), rsub(size(remove)))
        nsum = product_shape(x%shape(remove))
        do zi = 1, z%nelem()
            if (z%rank() > 0) call decode_index(zi, z%shape, zsub)
            xsub = 1
            do p = 1, size(rest)
                xsub(rest(p)) = zsub(p)
            end do
            do ri = 1, nsum
                if (size(remove) > 0) call decode_index(ri, x%shape(remove), rsub)
                do p = 1, size(remove)
                    xsub(remove(p)) = rsub(p)
                end do
                z%data(zi) = z%data(zi) + x%data(linear_index(xsub, x%shape))
            end do
        end do
    end function margin_tensor

    function delta_tensor(shape, axis, mark, by) result(z)
        integer, intent(in) :: shape(:)
        character(len=*), intent(in) :: axis(:)
        character(len=*), intent(in), optional :: mark
        integer, intent(in), optional :: by(:)
        type(tensor_t) :: z
        integer, allocatable :: b(:), nb(:), sub(:), outsub(:)
        character(len=axis_name_len) :: mk
        integer :: k, p, idx

        if (size(shape) /= size(axis)) error stop "delta_tensor: metadata mismatch"
        mk = "'"
        if (present(mark)) mk = trim(mark)
        if (present(by)) then
            allocate(b(size(by)))
            b = by
        else
            allocate(b(0))
        end if
        nb = complement(size(shape), b)
        allocate(z%shape(2*size(nb) + size(b)), z%axis(2*size(nb) + size(b)))
        if (size(nb) > 0) then
            z%shape(1:size(nb)) = shape(nb)
            z%shape(size(nb)+1:2*size(nb)) = shape(nb)
            z%axis(1:size(nb)) = axis(nb)
            do k = 1, size(nb)
                z%axis(size(nb)+k) = marked_name(axis(nb(k)), mk)
            end do
        end if
        if (size(b) > 0) then
            z%shape(2*size(nb)+1:) = shape(b)
            z%axis(2*size(nb)+1:) = axis(b)
        end if
        allocate(z%data(product_shape(z%shape)))
        z%data = (0.0_dp, 0.0_dp)
        allocate(sub(size(shape)), outsub(z%rank()))
        do idx = 1, product_shape(shape)
            call decode_index(idx, shape, sub)
            p = 0
            do k = 1, size(nb)
                p = p + 1
                outsub(p) = sub(nb(k))
            end do
            do k = 1, size(nb)
                p = p + 1
                outsub(p) = sub(nb(k))
            end do
            do k = 1, size(b)
                p = p + 1
                outsub(p) = sub(b(k))
            end do
            z%data(linear_index(outsub, z%shape)) = (1.0_dp, 0.0_dp)
        end do
    end function delta_tensor

    function diag_tensor(x, mark, by) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in), optional :: mark
        integer, intent(in), optional :: by(:)
        type(tensor_t) :: z
        integer, allocatable :: b(:), nb(:), xsub(:), zsub(:)
        character(len=axis_name_len) :: mk
        integer :: idx, k, p

        mk = "'"
        if (present(mark)) mk = trim(mark)
        if (present(by)) then
            allocate(b(size(by)))
            b = by
        else
            allocate(b(0))
        end if
        nb = complement(x%rank(), b)
        allocate(z%shape(2*size(nb)+size(b)), z%axis(2*size(nb)+size(b)))
        if (size(nb) > 0) then
            z%shape(1:size(nb)) = x%shape(nb)
            z%shape(size(nb)+1:2*size(nb)) = x%shape(nb)
            z%axis(1:size(nb)) = x%axis(nb)
            do k = 1, size(nb)
                z%axis(size(nb)+k) = marked_name(x%axis(nb(k)), mk)
            end do
        end if
        if (size(b) > 0) then
            z%shape(2*size(nb)+1:) = x%shape(b)
            z%axis(2*size(nb)+1:) = x%axis(b)
        end if
        allocate(z%data(product_shape(z%shape)))
        z%data = (0.0_dp, 0.0_dp)
        allocate(xsub(x%rank()), zsub(z%rank()))
        do idx = 1, x%nelem()
            call decode_index(idx, x%shape, xsub)
            p = 0
            do k = 1, size(nb)
                p = p + 1
                zsub(p) = xsub(nb(k))
            end do
            do k = 1, size(nb)
                p = p + 1
                zsub(p) = xsub(nb(k))
            end do
            do k = 1, size(b)
                p = p + 1
                zsub(p) = xsub(b(k))
            end do
            z%data(linear_index(zsub, z%shape)) = x%data(idx)
        end do
    end function diag_tensor

    function tripledelta_tensor(shape, axis, mark1, mark2) result(z)
        integer, intent(in) :: shape(:)
        character(len=*), intent(in) :: axis(:)
        character(len=*), intent(in), optional :: mark1, mark2
        type(tensor_t) :: z
        character(len=axis_name_len) :: m1, m2
        integer, allocatable :: sub(:), zsub(:)
        integer :: idx, k, n

        if (size(shape) /= size(axis)) error stop "tripledelta_tensor: metadata mismatch"
        m1 = "'"
        m2 = "*"
        if (present(mark1)) m1 = trim(mark1)
        if (present(mark2)) m2 = trim(mark2)
        n = size(shape)
        allocate(z%shape(3*n), z%axis(3*n), sub(n), zsub(3*n))
        z%shape(1:n) = shape
        z%shape(n+1:2*n) = shape
        z%shape(2*n+1:3*n) = shape
        z%axis(1:n) = axis
        do k = 1, n
            z%axis(n+k) = marked_name(axis(k), m1)
            z%axis(2*n+k) = marked_name(axis(k), m2)
        end do
        allocate(z%data(product_shape(z%shape)))
        z%data = (0.0_dp, 0.0_dp)
        do idx = 1, product_shape(shape)
            call decode_index(idx, shape, sub)
            zsub(1:n) = sub
            zsub(n+1:2*n) = sub
            zsub(2*n+1:3*n) = sub
            z%data(linear_index(zsub, z%shape)) = (1.0_dp, 0.0_dp)
        end do
    end function tripledelta_tensor

    function one_tensor(shape, axis) result(z)
        integer, intent(in) :: shape(:)
        character(len=*), intent(in) :: axis(:)
        type(tensor_t) :: z

        z = tensor_ones(shape, axis)
    end function one_tensor

    function diagmul_tensor(x, ix, d, jd) result(z)
        type(tensor_t), intent(in) :: x, d
        character(len=*), intent(in) :: ix(:), jd(:)
        type(tensor_t) :: z, dr
        integer :: k

        if (size(ix) /= size(jd)) error stop "diagmul_tensor: axis-count mismatch"
        dr = d
        do k = 1, size(ix)
            dr = rename_axis(dr, jd(k), ix(k))
        end do
        z = elem_mul_tensor(x, dr)
    end function diagmul_tensor

    function add_tensor(x, y) result(z)
        type(tensor_t), intent(in) :: x, y
        type(tensor_t) :: z

        z = elementwise_binary(x, y, 1)
    end function add_tensor

    function sub_tensor(x, y) result(z)
        type(tensor_t), intent(in) :: x, y
        type(tensor_t) :: z

        z = elementwise_binary(x, y, 2)
    end function sub_tensor

    function elem_mul_tensor(x, y) result(z)
        type(tensor_t), intent(in) :: x, y
        type(tensor_t) :: z

        z = elementwise_binary(x, y, 3)
    end function elem_mul_tensor

    function elem_div_tensor(x, y) result(z)
        type(tensor_t), intent(in) :: x, y
        type(tensor_t) :: z

        z = elementwise_binary(x, y, 4)
    end function elem_div_tensor

    function scale_tensor(x, a) result(z)
        type(tensor_t), intent(in) :: x
        complex(dp), intent(in) :: a
        type(tensor_t) :: z

        z = x
        z%data = a * x%data
    end function scale_tensor

    function repeat_tensor(x, times, pos, name) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: times
        integer, intent(in), optional :: pos
        character(len=*), intent(in), optional :: name
        type(tensor_t) :: z, tmp
        integer :: p, k, idx, block
        character(len=axis_name_len) :: nm

        if (times < 1) error stop "repeat_tensor: times must be positive"
        p = 1
        if (present(pos)) p = pos
        if (p < 1 .or. p > x%rank()+1) error stop "repeat_tensor: bad insertion position"
        nm = "i"
        if (present(name)) nm = trim(name)
        allocate(tmp%shape(x%rank()+1), tmp%axis(x%rank()+1), tmp%data(x%nelem()*times))
        if (p > 1) then
            tmp%shape(1:p-1) = x%shape(1:p-1)
            tmp%axis(1:p-1) = x%axis(1:p-1)
        end if
        tmp%shape(p) = times
        tmp%axis(p) = nm
        if (p <= x%rank()) then
            tmp%shape(p+1:) = x%shape(p:)
            tmp%axis(p+1:) = x%axis(p:)
        end if
        block = product_shape(x%shape(1:p-1))
        idx = 0
        do k = 1, x%nelem()/block
            call repeat_block(x%data((k-1)*block+1:k*block), times, tmp%data, idx)
        end do
        z = tmp
    end function repeat_tensor

    function slice_tensor(x, axis_pos, select, drop) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: axis_pos
        integer, intent(in) :: select(:)
        logical, intent(in), optional :: drop
        type(tensor_t) :: z
        logical :: dr
        integer, allocatable :: zsub(:), xsub(:), zshape(:)
        character(len=axis_name_len), allocatable :: zaxis(:)
        integer :: zi, k, p

        if (axis_pos < 1 .or. axis_pos > x%rank()) error stop "slice_tensor: bad axis"
        if (size(select) < 1) error stop "slice_tensor: empty selection"
        if (any(select < 1) .or. any(select > x%shape(axis_pos))) error stop "slice_tensor: selection out of range"
        dr = .false.
        if (present(drop)) dr = drop
        if (dr .and. size(select) == 1) then
            allocate(zshape(x%rank()-1), zaxis(x%rank()-1))
            p = 0
            do k = 1, x%rank()
                if (k /= axis_pos) then
                    p = p + 1
                    zshape(p) = x%shape(k)
                    zaxis(p) = x%axis(k)
                end if
            end do
        else
            allocate(zshape(x%rank()), zaxis(x%rank()))
            zshape = x%shape
            zshape(axis_pos) = size(select)
            zaxis = x%axis
        end if
        allocate(z%shape(size(zshape)), z%axis(size(zaxis)), z%data(product_shape(zshape)))
        z%shape = zshape
        z%axis = zaxis
        allocate(zsub(z%rank()), xsub(x%rank()))
        do zi = 1, z%nelem()
            if (z%rank() > 0) call decode_index(zi, z%shape, zsub)
            if (dr .and. size(select) == 1) then
                p = 0
                do k = 1, x%rank()
                    if (k == axis_pos) then
                        xsub(k) = select(1)
                    else
                        p = p + 1
                        xsub(k) = zsub(p)
                    end if
                end do
            else
                xsub = zsub
                xsub(axis_pos) = select(zsub(axis_pos))
            end if
            z%data(zi) = x%data(linear_index(xsub, x%shape))
        end do
    end function slice_tensor

    function undrop_tensor(x, name, pos) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: name
        integer, intent(in), optional :: pos
        type(tensor_t) :: z
        integer :: p

        p = 1
        if (present(pos)) p = pos
        if (p < 1 .or. p > x%rank()+1) error stop "undrop_tensor: bad position"
        allocate(z%shape(x%rank()+1), z%axis(x%rank()+1), z%data(x%nelem()))
        if (p > 1) then
            z%shape(1:p-1) = x%shape(1:p-1)
            z%axis(1:p-1) = x%axis(1:p-1)
        end if
        z%shape(p) = 1
        z%axis(p) = trim(name)
        if (p <= x%rank()) then
            z%shape(p+1:) = x%shape(p:)
            z%axis(p+1:) = x%axis(p:)
        end if
        z%data = x%data
    end function undrop_tensor

    function untensor_tensor(x, collapse, name, pos) result(z)
        type(tensor_t), intent(in) :: x
        integer, intent(in) :: collapse(:)
        character(len=*), intent(in) :: name
        integer, intent(in), optional :: pos
        type(tensor_t) :: z, tmp
        integer, allocatable :: rest(:), order(:)
        integer :: p, k, outp

        if (size(collapse) == 0) then
            z = x
            return
        end if
        rest = complement(x%rank(), collapse)
        order = [collapse, rest]
        tmp = reorder_tensor_pos(x, order)
        p = 1
        if (present(pos)) p = pos
        if (p < 1 .or. p > size(rest)+1) error stop "untensor_tensor: bad insertion position"
        allocate(z%shape(size(rest)+1), z%axis(size(rest)+1), z%data(x%nelem()))
        outp = 0
        do k = 1, p-1
            outp = outp + 1
            z%shape(outp) = x%shape(rest(k))
            z%axis(outp) = x%axis(rest(k))
        end do
        outp = outp + 1
        z%shape(outp) = product_shape(x%shape(collapse))
        z%axis(outp) = trim(name)
        do k = p, size(rest)
            outp = outp + 1
            z%shape(outp) = x%shape(rest(k))
            z%axis(outp) = x%axis(rest(k))
        end do
        ! tmp has collapsed axes first. Move that compound first axis to requested position.
        z%data = reshape_collapsed_data(tmp, x%shape(collapse), x%shape(rest), p)
    end function untensor_tensor

    function bind_tensor(a, da, b, db) result(z)
        type(tensor_t), intent(in) :: a, b
        integer, intent(in) :: da, db
        type(tensor_t) :: z, br
        integer, allocatable :: bfirst(:), bmap(:), zsub(:), asub(:), bsub(:)
        integer :: k, p, zi, ai, bi

        if (da < 1 .or. da > a%rank() .or. db < 1 .or. db > b%rank()) error stop "bind_tensor: bad axis"
        if (a%rank() /= b%rank()) error stop "bind_tensor: ranks differ"
        allocate(bfirst(b%rank()), bmap(a%rank()))
        bfirst(da) = db
        do k = 1, a%rank()
            if (k == da) cycle
            p = b%axis_pos(a%axis(k))
            if (p <= 0 .or. p == db) error stop "bind_tensor: non-bound axes do not match"
            if (a%shape(k) /= b%shape(p)) error stop "bind_tensor: incompatible non-bound dimensions"
            bfirst(k) = p
        end do
        br = reorder_tensor_pos(b, bfirst)
        allocate(z%shape(a%rank()), z%axis(a%rank()), z%data(a%nelem()+b%nelem()))
        z%shape = a%shape
        z%shape(da) = a%shape(da) + br%shape(da)
        z%axis = a%axis
        allocate(zsub(z%rank()), asub(a%rank()), bsub(br%rank()))
        do zi = 1, z%nelem()
            call decode_index(zi, z%shape, zsub)
            if (zsub(da) <= a%shape(da)) then
                asub = zsub
                ai = linear_index(asub, a%shape)
                z%data(zi) = a%data(ai)
            else
                bsub = zsub
                bsub(da) = bsub(da) - a%shape(da)
                bi = linear_index(bsub, br%shape)
                z%data(zi) = br%data(bi)
            end if
        end do
    end function bind_tensor

    function einstein_pair(x, y, only, by) result(z)
        type(tensor_t), intent(in) :: x, y
        character(len=*), intent(in), optional :: only(:), by(:)
        type(tensor_t) :: z
        character(len=axis_name_len), allocatable :: common(:), bnames(:)
        integer :: i, j, n
        logical :: use_name

        allocate(common(min(x%rank(), y%rank())))
        n = 0
        do i = 1, x%rank()
            do j = 1, y%rank()
                if (trim(x%axis(i)) == trim(y%axis(j))) then
                    use_name = .true.
                    if (present(only)) use_name = any(trim(x%axis(i)) == only)
                    if (present(by)) then
                        if (any(trim(x%axis(i)) == by)) use_name = .false.
                    end if
                    if (use_name) then
                        n = n + 1
                        common(n) = x%axis(i)
                    end if
                end if
            end do
        end do
        common = common(:n)
        if (present(by)) then
            call shared_names(x, y, by, bnames)
            z = mul_tensor_names(x, common, y, common, bnames)
        else
            z = mul_tensor_names(x, common, y, common)
        end if
    end function einstein_pair

    function riemann_pair(x, y, only, by) result(z)
        type(tensor_t), intent(in) :: x, y
        character(len=*), intent(in), optional :: only(:), by(:)
        type(tensor_t) :: z
        character(len=axis_name_len), allocatable :: xn(:), yn(:), bnames(:)
        integer :: i, j, n
        logical :: use_name

        allocate(xn(min(x%rank(), y%rank())), yn(min(x%rank(), y%rank())))
        n = 0
        do i = 1, x%rank()
            do j = 1, y%rank()
                if (trim(contraname(x%axis(i))) == trim(y%axis(j))) then
                    use_name = .true.
                    if (present(only)) then
                        use_name = any(trim(x%axis(i)) == only) .or. &
                            any(trim(contraname(x%axis(i))) == only)
                    end if
                    if (present(by)) then
                        if (any(trim(x%axis(i)) == by) .or. &
                            any(trim(contraname(x%axis(i))) == by)) use_name = .false.
                    end if
                    if (use_name) then
                        n = n + 1
                        xn(n) = x%axis(i)
                        yn(n) = y%axis(j)
                    end if
                end if
            end do
        end do
        xn = xn(:n)
        yn = yn(:n)
        if (present(by)) then
            call shared_names(x, y, by, bnames)
            z = mul_tensor_names(x, xn, y, yn, bnames)
        else
            z = mul_tensor_names(x, xn, y, yn)
        end if
    end function riemann_pair

    function contraname_tensor(x) result(y)
        type(tensor_t), intent(in) :: x
        type(tensor_t) :: y
        integer :: k

        y = x
        do k = 1, y%rank()
            y%axis(k) = contraname(y%axis(k))
        end do
    end function contraname_tensor

    pure logical function is_covariate_tensor(x) result(ans)
        type(tensor_t), intent(in) :: x
        integer :: k

        ans = .true.
        do k = 1, x%rank()
            if (x%axis(k)(1:1) == '^') then
                ans = .false.
                return
            end if
        end do
    end function is_covariate_tensor

    pure logical function is_contravariate_tensor(x) result(ans)
        type(tensor_t), intent(in) :: x
        integer :: k

        ans = .true.
        do k = 1, x%rank()
            if (x%axis(k)(1:1) /= '^') then
                ans = .false.
                return
            end if
        end do
    end function is_contravariate_tensor

    function elementwise_binary(x, y, op) result(z)
        type(tensor_t), intent(in) :: x, y
        integer, intent(in) :: op
        type(tensor_t) :: z
        integer, allocatable :: yonly(:), zsub(:), xsub(:), ysub(:)
        integer :: i, j, zi, xi, yi, n, p

        allocate(yonly(y%rank()))
        n = 0
        do j = 1, y%rank()
            p = x%axis_pos(y%axis(j))
            if (p == 0) then
                n = n + 1
                yonly(n) = j
            else
                if (p < 0) error stop "elementwise_binary: ambiguous x axis"
                if (x%shape(p) /= y%shape(j)) error stop "elementwise_binary: shared dimension mismatch"
            end if
        end do
        yonly = yonly(:n)
        allocate(z%shape(x%rank()+n), z%axis(x%rank()+n))
        if (x%rank() > 0) then
            z%shape(1:x%rank()) = x%shape
            z%axis(1:x%rank()) = x%axis
        end if
        if (n > 0) then
            z%shape(x%rank()+1:) = y%shape(yonly)
            z%axis(x%rank()+1:) = y%axis(yonly)
        end if
        allocate(z%data(product_shape(z%shape)), zsub(z%rank()), xsub(x%rank()), ysub(y%rank()))
        do zi = 1, z%nelem()
            if (z%rank() > 0) call decode_index(zi, z%shape, zsub)
            if (x%rank() > 0) xsub = zsub(1:x%rank())
            do j = 1, y%rank()
                p = x%axis_pos(y%axis(j))
                if (p > 0) then
                    ysub(j) = zsub(p)
                else
                    do i = 1, n
                        if (yonly(i) == j) then
                            ysub(j) = zsub(x%rank()+i)
                            exit
                        end if
                    end do
                end if
            end do
            if (x%rank() == 0) then
                xi = 1
            else
                xi = linear_index(xsub, x%shape)
            end if
            if (y%rank() == 0) then
                yi = 1
            else
                yi = linear_index(ysub, y%shape)
            end if
            select case (op)
            case (1)
                z%data(zi) = x%data(xi) + y%data(yi)
            case (2)
                z%data(zi) = x%data(xi) - y%data(yi)
            case (3)
                z%data(zi) = x%data(xi) * y%data(yi)
            case (4)
                z%data(zi) = x%data(xi) / y%data(yi)
            case default
                error stop "elementwise_binary: unknown operation"
            end select
        end do
    end function elementwise_binary

    pure logical function has_overlap(a, b) result(overlap)
        integer, intent(in) :: a(:), b(:)
        integer :: k

        overlap = .false.
        do k = 1, size(a)
            if (any(b == a(k))) then
                overlap = .true.
                return
            end if
        end do
    end function has_overlap

    subroutine repeat_block(block_in, times, out, idx)
        complex(dp), intent(in) :: block_in(:)
        integer, intent(in) :: times
        complex(dp), intent(inout) :: out(:)
        integer, intent(inout) :: idx
        integer :: k, n

        n = size(block_in)
        do k = 1, times
            out(idx+1:idx+n) = block_in
            idx = idx + n
        end do
    end subroutine repeat_block

    function reshape_collapsed_data(tmp, cshape, rshape, pos) result(out)
        type(tensor_t), intent(in) :: tmp
        integer, intent(in) :: cshape(:), rshape(:), pos
        complex(dp), allocatable :: out(:)
        integer, allocatable :: tshape(:), osub(:), tsub(:)
        integer :: zi, k, p, cprod

        cprod = product_shape(cshape)
        allocate(tshape(size(rshape)+1))
        if (pos > 1) tshape(1:pos-1) = rshape(1:pos-1)
        tshape(pos) = cprod
        if (pos <= size(rshape)) tshape(pos+1:) = rshape(pos:)
        allocate(out(product_shape(tshape)), osub(size(tshape)), tsub(tmp%rank()))
        do zi = 1, size(out)
            call decode_index(zi, tshape, osub)
            ! tmp axes are original collapsed axes followed by rest axes.
            call decode_index(osub(pos), cshape, tsub(1:size(cshape)))
            p = size(cshape)
            do k = 1, size(rshape)
                p = p + 1
                if (k < pos) then
                    tsub(p) = osub(k)
                else
                    tsub(p) = osub(k+1)
                end if
            end do
            out(zi) = tmp%data(linear_index(tsub, tmp%shape))
        end do
    end function reshape_collapsed_data

    subroutine shared_names(x, y, requested, shared)
        type(tensor_t), intent(in) :: x, y
        character(len=*), intent(in) :: requested(:)
        character(len=axis_name_len), allocatable, intent(out) :: shared(:)
        integer :: k, n

        allocate(shared(size(requested)))
        n = 0
        do k = 1, size(requested)
            if (x%axis_pos(requested(k)) > 0 .and. y%axis_pos(requested(k)) > 0) then
                n = n + 1
                shared(n) = trim(requested(k))
            end if
        end do
        shared = shared(:n)
    end subroutine shared_names

end module tensora_core
