! SPDX-License-Identifier: GPL-2.0-or-later
module tensora_stats
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use tensora_kinds, only : dp, axis_name_len
    use tensora_types, only : tensor_t
    use tensora_index, only : complement, positions_by_name, contraname, product_shape
    use tensora_core, only : margin_tensor, sub_tensor, scale_tensor, mark_tensor, mul_tensor_names
    use tensora_linalg, only : inv_tensor
    implicit none
    private

    public :: norm_tensor, norm_tensor_along, mean_tensor, var_tensor, cov_tensor
    public :: drag_tensor

contains

    real(dp) function norm_tensor(x) result(nrm)
        type(tensor_t), intent(in) :: x

        nrm = sqrt(sum(abs(x%data)**2))
    end function norm_tensor

    function norm_tensor_along(x, remove_names) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: remove_names(:)
        type(tensor_t) :: z
        integer, allocatable :: remove(:)

        remove = positions_by_name(x, remove_names)
        z = margin_tensor(abs_square_tensor(x), remove)
        z%data = sqrt(z%data)
    end function norm_tensor_along

    function mean_tensor(x, along_names, na_rm) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: along_names(:)
        logical, intent(in), optional :: na_rm
        type(tensor_t) :: z
        integer, allocatable :: along(:), rest(:), zsub(:), xsub(:), asub(:)
        integer :: zi, ai, p, nsum, count_valid, xi
        logical :: skip_bad
        complex(dp) :: total

        skip_bad = .false.
        if (present(na_rm)) skip_bad = na_rm
        along = positions_by_name(x, along_names)
        rest = complement(x%rank(), along)
        allocate(z%shape(size(rest)), z%axis(size(rest)), z%data(product_shape(x%shape(rest))))
        if (size(rest) > 0) then
            z%shape = x%shape(rest)
            z%axis = x%axis(rest)
        end if
        allocate(zsub(z%rank()), xsub(x%rank()), asub(size(along)))
        nsum = product_shape(x%shape(along))
        do zi = 1, z%nelem()
            if (z%rank() > 0) call decode_local(zi, z%shape, zsub)
            xsub = 1
            do p = 1, size(rest)
                xsub(rest(p)) = zsub(p)
            end do
            total = (0.0_dp, 0.0_dp)
            count_valid = 0
            do ai = 1, nsum
                if (size(along) > 0) call decode_local(ai, x%shape(along), asub)
                do p = 1, size(along)
                    xsub(along(p)) = asub(p)
                end do
                xi = linear_local(xsub, x%shape)
                if (complex_is_finite(x%data(xi))) then
                    total = total + x%data(xi)
                    count_valid = count_valid + 1
                else if (.not. skip_bad) then
                    error stop "mean_tensor: non-finite value; set na_rm=.true. to omit"
                end if
            end do
            if (count_valid == 0) then
                z%data(zi) = cmplx(ieee_nan(), 0.0_dp, dp)
            else
                z%data(zi) = total / real(count_valid,dp)
            end if
        end do
    end function mean_tensor

    function var_tensor(x, along_names, by_names, mark, na_rm) result(z)
        type(tensor_t), intent(in) :: x
        character(len=*), intent(in) :: along_names(:)
        character(len=*), intent(in), optional :: by_names(:)
        character(len=*), intent(in), optional :: mark
        logical, intent(in), optional :: na_rm
        type(tensor_t) :: z, mu, centered, y
        character(len=axis_name_len) :: mk
        character(len=axis_name_len), allocatable :: tomark(:), b(:)
        integer, allocatable :: along(:), bypos(:), keep(:)
        integer :: k, n, df
        logical :: skip_bad

        skip_bad = .false.
        if (present(na_rm)) skip_bad = na_rm
        if (skip_bad) error stop "var_tensor: na_rm=.true. is not yet supported for pairwise tensor covariance"
        mk = "'"
        if (present(mark)) mk = trim(mark)
        along = positions_by_name(x, along_names)
        if (present(by_names)) then
            bypos = positions_by_name(x, by_names, missing_ok=.true.)
            bypos = pack(bypos, bypos > 0)
        else
            allocate(bypos(0))
        end if
        keep = complement(x%rank(), [along, bypos])
        allocate(tomark(size(keep)))
        do k = 1, size(keep)
            tomark(k) = x%axis(keep(k))
        end do
        mu = mean_tensor(x, along_names)
        centered = sub_tensor(x, mu)
        y = mark_tensor(centered, mk, tomark)
        df = product_shape(x%shape(along)) - 1
        if (df <= 0) error stop "var_tensor: need at least two observations"
        if (present(by_names)) then
            allocate(b(size(by_names)))
            n = 0
            do k = 1, size(by_names)
                if (x%axis_pos(by_names(k)) > 0 .and. y%axis_pos(by_names(k)) > 0) then
                    n = n + 1
                    b(n) = trim(by_names(k))
                end if
            end do
            b = b(:n)
            z = mul_tensor_names(centered, along_names, y, along_names, b)
        else
            z = mul_tensor_names(centered, along_names, y, along_names)
        end if
        z = scale_tensor(z, cmplx(1.0_dp/real(df,dp),0.0_dp,dp))
    end function var_tensor

    function cov_tensor(x, y, along_names, by_names, mark) result(z)
        type(tensor_t), intent(in) :: x, y
        character(len=*), intent(in) :: along_names(:)
        character(len=*), intent(in), optional :: by_names(:)
        character(len=*), intent(in), optional :: mark
        type(tensor_t) :: z, mux, muy, xc, yc, ym
        character(len=axis_name_len) :: mk
        character(len=axis_name_len), allocatable :: tomark(:), b(:)
        integer, allocatable :: ax(:), ay(:), byp(:), keepy(:)
        integer :: k, n, df

        mk = "'"
        if (present(mark)) mk = trim(mark)
        ax = positions_by_name(x, along_names)
        ay = positions_by_name(y, along_names)
        allocate(keepy(0))
        if (any(x%shape(ax) /= y%shape(ay))) error stop "cov_tensor: observation dimensions differ"
        if (present(by_names)) then
            byp = positions_by_name(y, by_names, missing_ok=.true.)
            keepy = complement(y%rank(), [ay, pack(byp,byp>0)])
        else
            allocate(byp(0))
            keepy = complement(y%rank(), ay)
        end if
        allocate(tomark(size(keepy)))
        do k = 1, size(keepy)
            tomark(k) = y%axis(keepy(k))
        end do
        mux = mean_tensor(x, along_names)
        muy = mean_tensor(y, along_names)
        xc = sub_tensor(x, mux)
        yc = sub_tensor(y, muy)
        ym = mark_tensor(yc, mk, tomark)
        df = product_shape(x%shape(ax)) - 1
        if (df <= 0) error stop "cov_tensor: need at least two observations"
        if (present(by_names)) then
            allocate(b(size(by_names)))
            n = 0
            do k = 1, size(by_names)
                if (xc%axis_pos(by_names(k)) > 0 .and. ym%axis_pos(by_names(k)) > 0) then
                    n = n + 1
                    b(n) = trim(by_names(k))
                end if
            end do
            b = b(:n)
            z = mul_tensor_names(xc, along_names, ym, along_names, b)
        else
            z = mul_tensor_names(xc, along_names, ym, along_names)
        end if
        z = scale_tensor(z, cmplx(1.0_dp/real(df,dp),0.0_dp,dp))
    end function cov_tensor

    function drag_tensor(x, g, dims) result(z)
        type(tensor_t), intent(in) :: x, g
        character(len=*), intent(in) :: dims(:)
        type(tensor_t) :: z, gcov, gcon, metric
        character(len=axis_name_len) :: na, nb
        integer :: k
        logical :: gcovariant

        if (g%rank() /= 2) error stop "drag_tensor: metric must have rank two"
        gcovariant = g%axis(1)(1:1) /= '^' .and. g%axis(2)(1:1) /= '^'
        if (.not. gcovariant) then
            if (.not. (g%axis(1)(1:1) == '^' .and. g%axis(2)(1:1) == '^')) then
                error stop "drag_tensor: metric axes must be both covariant or both contravariant"
            end if
        end if
        if (gcovariant) then
            gcov = g
            gcon = inv_tensor(g, [trim(g%axis(1))])
        else
            gcon = g
            gcov = inv_tensor(g, [trim(g%axis(1))])
        end if
        z = x
        do k = 1, size(dims)
            na = trim(dims(k))
            if (z%axis_pos(na) <= 0) error stop "drag_tensor: tensor axis not found"
            nb = contraname(na)
            if (na(1:1) == '^') then
                metric = gcov
            else
                metric = gcon
            end if
            metric%axis(1) = na
            metric%axis(2) = nb
            z = mul_tensor_names(z, [na], metric, [na])
        end do
    end function drag_tensor

    function abs_square_tensor(x) result(z)
        type(tensor_t), intent(in) :: x
        type(tensor_t) :: z

        z = x
        z%data = cmplx(abs(x%data)**2, 0.0_dp, dp)
    end function abs_square_tensor

    pure logical function complex_is_finite(z) result(ok)
        complex(dp), intent(in) :: z

        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function complex_is_finite

    pure real(dp) function ieee_nan() result(x)
        use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function ieee_nan

    pure subroutine decode_local(idx, shape, sub)
        integer, intent(in) :: idx, shape(:)
        integer, intent(out) :: sub(:)
        integer :: q, k

        q = idx - 1
        do k = 1, size(shape)
            sub(k) = mod(q,shape(k)) + 1
            q = q / shape(k)
        end do
    end subroutine decode_local

    pure integer function linear_local(sub, shape) result(idx)
        integer, intent(in) :: sub(:), shape(:)
        integer :: k, stride

        idx = 1
        stride = 1
        do k = 1, size(shape)
            idx = idx + (sub(k)-1)*stride
            stride = stride*shape(k)
        end do
    end function linear_local

end module tensora_stats
