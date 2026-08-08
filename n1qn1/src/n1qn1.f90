! Modern Fortran translation of Scilab n1qn1/n1qn1a and majour.
!
! Original algorithm:
!   Copyright (C) INRIA, Claude Lemarechal, 1987.
!   Copyright (C) 2012-2016 Scilab Enterprises.
!
! The original source is available under CeCILL v2.1 and, for the Scilab
! source files, under GNU GPL v2.0 pursuant to article 5.3.4 of CeCILL v2.1.
! This translation is distributed under CeCILL v2.1. See LICENCE and NOTICE.
module n1qn1_module
    use, intrinsic :: iso_fortran_env, only : real64, int64, output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    integer, parameter, public :: dp = real64

    integer, parameter, public :: n1qn1_success = 0
    integer, parameter, public :: n1qn1_max_iterations = 1
    integer, parameter, public :: n1qn1_max_evaluations = 2
    integer, parameter, public :: n1qn1_not_descent = 3
    integer, parameter, public :: n1qn1_rank_loss = 4
    integer, parameter, public :: n1qn1_user_stop = 5
    integer, parameter, public :: n1qn1_invalid_input = -1
    integer, parameter, public :: n1qn1_nonfinite = -2

    type, public :: n1qn1_control_t
        real(dp) :: epsilon = epsilon(1.0_dp)
        integer :: max_iterations = 100
        integer :: max_evaluations = 100
        integer :: verbosity = 0
    end type n1qn1_control_t

    type, public :: n1qn1_result_t
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: gradient(:)
        real(dp), allocatable :: hessian(:, :)
        real(dp), allocatable :: c_hess(:)
        real(dp), allocatable :: factor(:)
        real(dp) :: value = huge(1.0_dp)
        real(dp) :: gradient_norm_squared = huge(1.0_dp)
        integer :: iterations = 0
        integer :: function_evaluations = 0
        integer :: gradient_evaluations = 0
        integer :: status = n1qn1_invalid_input
        character(len=:), allocatable :: message
    end type n1qn1_result_t

    abstract interface
        function n1qn1_objective(x, user_data) result(value)
            import :: dp
            real(dp), intent(in) :: x(:)
            class(*), intent(inout), optional :: user_data
            real(dp) :: value
        end function n1qn1_objective

        subroutine n1qn1_gradient(x, gradient, user_data)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: gradient(:)
            class(*), intent(inout), optional :: user_data
        end subroutine n1qn1_gradient

        subroutine n1qn1_progress(iteration, evaluations, x, value, gradient, stop, user_data)
            import :: dp
            integer, intent(in) :: iteration
            integer, intent(in) :: evaluations
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: value
            real(dp), intent(in) :: gradient(:)
            logical, intent(out) :: stop
            class(*), intent(inout), optional :: user_data
        end subroutine n1qn1_progress
    end interface

    public :: n1qn1_minimize
    public :: n1qn1_objective, n1qn1_gradient, n1qn1_progress
    public :: n1qn1_workspace_size
    public :: n1qn1_packed_size
    public :: pack_lower
    public :: unpack_lower
    public :: factor_to_hessian
    public :: n1qn1_status_message

contains

    pure integer function n1qn1_packed_size(n) result(size_h)
        integer, intent(in) :: n
        integer(int64) :: value
        if (n <= 0) then
            size_h = 0
            return
        end if
        value = int(n, int64) * int(n + 1, int64) / 2_int64
        if (value > int(huge(size_h), int64)) then
            size_h = -1
        else
            size_h = int(value)
        end if
    end function n1qn1_packed_size

    pure integer function n1qn1_workspace_size(n) result(size_zm)
        integer, intent(in) :: n
        integer(int64) :: value
        if (n <= 0) then
            size_zm = 0
            return
        end if
        value = int(n, int64) * int(n + 13, int64) / 2_int64
        if (value > int(huge(size_zm), int64)) then
            size_zm = -1
        else
            size_zm = int(value)
        end if
    end function n1qn1_workspace_size

    subroutine n1qn1_minimize(objective, gradient_function, x0, result, control, scale, &
                              initial_hessian, initial_factor, user_data, progress)
        procedure(n1qn1_objective) :: objective
        procedure(n1qn1_gradient) :: gradient_function
        real(dp), intent(in) :: x0(:)
        type(n1qn1_result_t), intent(out) :: result
        type(n1qn1_control_t), intent(in), optional :: control
        real(dp), intent(in), optional :: scale(:)
        real(dp), intent(in), optional :: initial_hessian(:, :)
        real(dp), intent(in), optional :: initial_factor(:)
        class(*), intent(inout), optional :: user_data
        procedure(n1qn1_progress), optional :: progress

        type(n1qn1_control_t) :: ctrl
        real(dp), allocatable :: x(:), g(:), scale_work(:), h(:)
        real(dp), allocatable :: d(:), w(:), xa(:), ga(:), xb(:), gb(:)
        real(dp) :: f, acc
        integer :: n, nh, niter, nsim, mode, status

        ctrl = n1qn1_control_t()
        if (present(control)) ctrl = control
        n = size(x0)

        if (n <= 0) then
            call set_error(result, n1qn1_invalid_input, "x0 must contain at least one element")
            return
        end if
        if (n1qn1_packed_size(n) < 0 .or. n1qn1_workspace_size(n) < 0) then
            call set_error(result, n1qn1_invalid_input, &
                "problem dimension is too large for default-integer packed indexing")
            return
        end if
        if (ctrl%epsilon < 0.0_dp .or. ctrl%max_iterations <= 0 .or. ctrl%max_evaluations <= 0) then
            call set_error(result, n1qn1_invalid_input, "invalid control parameter")
            return
        end if
        if (present(scale)) then
            if (size(scale) /= n .or. any(scale <= 0.0_dp) .or. any(.not. ieee_is_finite(scale))) then
                call set_error(result, n1qn1_invalid_input, &
                    "scale must be finite, positive, and have size(x0) elements")
                return
            end if
        end if
        if (present(initial_hessian) .and. present(initial_factor)) then
            call set_error(result, n1qn1_invalid_input, &
                "provide initial_hessian or initial_factor, not both")
            return
        end if

        nh = n1qn1_packed_size(n)
        allocate(x(n), g(n), scale_work(n), h(nh), d(n), w(n), xa(n), ga(n), xb(n), gb(n))
        x = x0
        g = 0.0_dp
        if (present(scale)) then
            scale_work = scale
        else
            scale_work = 0.1_dp
        end if
        h = 0.0_dp
        mode = 1

        if (present(initial_hessian)) then
            if (size(initial_hessian, 1) /= n .or. size(initial_hessian, 2) /= n) then
                call set_error(result, n1qn1_invalid_input, &
                    "initial_hessian must be a square size(x0) matrix")
                return
            end if
            if (any(.not. ieee_is_finite(initial_hessian))) then
                call set_error(result, n1qn1_invalid_input, "initial_hessian contains non-finite values")
                return
            end if
            call pack_lower(initial_hessian, h)
            mode = 2
        else if (present(initial_factor)) then
            if (size(initial_factor) /= nh .or. any(.not. ieee_is_finite(initial_factor))) then
                call set_error(result, n1qn1_invalid_input, &
                    "initial_factor must be finite and have n*(n+1)/2 elements")
                return
            end if
            h = initial_factor
            mode = 3
        end if

        f = huge(1.0_dp)
        acc = ctrl%epsilon
        niter = ctrl%max_iterations
        nsim = ctrl%max_evaluations

        call n1qn1_kernel(objective, gradient_function, n, x, f, g, scale_work, acc, mode, &
                          niter, nsim, ctrl%verbosity, h, d, w, xa, ga, xb, gb, status, &
                          user_data, progress)

        allocate(result%x(n), result%gradient(n), result%hessian(n, n), result%factor(nh))
        allocate(result%c_hess(n1qn1_workspace_size(n)))
        result%x = x
        result%gradient = g
        result%value = f
        result%gradient_norm_squared = acc
        result%iterations = niter
        result%function_evaluations = nsim
        result%gradient_evaluations = nsim
        result%status = status
        result%message = n1qn1_status_message(status)
        result%factor = h
        call factor_to_hessian(h, result%hessian)
        result%c_hess = 0.0_dp
        call pack_lower(result%hessian, result%c_hess(1:nh))


    end subroutine n1qn1_minimize

    subroutine n1qn1_kernel(objective, gradient_function, n, x, f, g, scale, acc, mode, &
                            niter, nsim, iprint, h, d, w, xa, ga, xb, gb, status, &
                            user_data, progress)
        procedure(n1qn1_objective) :: objective
        procedure(n1qn1_gradient) :: gradient_function
        integer, intent(in) :: n, mode, iprint
        real(dp), intent(inout) :: x(n), f, g(n), acc, h(:)
        real(dp), intent(in) :: scale(n)
        integer, intent(inout) :: niter, nsim
        real(dp), intent(inout) :: d(n), w(n), xa(n), ga(n), xb(n), gb(n)
        integer, intent(out) :: status
        class(*), intent(inout), optional :: user_data
        procedure(n1qn1_progress), optional :: progress

        integer :: i, j, k, i1, ii, ij, ik, jk, ni, ip, ir, np, nip
        integer :: itr, nfun, isfv, ial, indic, niter_limit, nsim_limit
        real(dp) :: c, v, cc, fa, fb, hh, gl1, gl2, dga, dgb, dff
        real(dp) :: fmin, gmin, step, stmin, stepbd, steplb
        logical :: stop_requested

        niter_limit = niter
        nsim_limit = nsim
        status = n1qn1_success

        indic = 4
        if (present(user_data)) then
            f = objective(x, user_data)
            call gradient_function(x, g, user_data)
        else
            f = objective(x)
            call gradient_function(x, g)
        end if
        if (.not. ieee_is_finite(f) .or. any(.not. ieee_is_finite(g))) indic = -1
        if (indic <= 0) then
            acc = 0.0_dp
            niter = 1
            nsim = 1
            if (indic == 0) then
                status = n1qn1_user_stop
            else
                status = n1qn1_nonfinite
            end if
            return
        end if

        nfun = 1
        itr = 0
        np = n + 1

        if (mode < 2) goto 20
        goto 60

20      c = 0.0_dp
        do i = 1, n
            c = max(c, abs(g(i) * scale(i)))
        end do
        if (c <= 0.0_dp) c = 1.0_dp
        k = n1qn1_packed_size(n)
        h(1:k) = 0.0_dp
        k = 1
        do i = 1, n
            h(k) = 0.01_dp * c / (scale(i) * scale(i))
            k = k + np - i
        end do
        goto 100

60      if (mode >= 3) goto 80
        k = n
        if (n > 1) goto 300
        if (h(1) > 0.0_dp) goto 305
        h(1) = 0.0_dp
        k = 0
        goto 305

300     np = n + 1
        ii = 1
        do i = 2, n
            hh = h(ii)
            ni = ii + np - i
            if (hh <= 0.0_dp) then
                h(ii) = 0.0_dp
                k = k - 1
                ii = ni + 1
                cycle
            end if
            ip = ii + 1
            ii = ni + 1
            jk = ii
            do ij = ip, ni
                v = h(ij) / hh
                do ik = ij, ni
                    h(jk) = h(jk) - h(ik) * v
                    jk = jk + 1
                end do
                h(ij) = v
            end do
        end do
        if (h(ii) <= 0.0_dp) then
            h(ii) = 0.0_dp
            k = k - 1
        end if

305     if (k >= n) goto 100
70      if (iprint /= 0) write(output_unit, '(a)') &
            "n1qn1: replacing a non-positive-definite initial Hessian"
        goto 20

80      k = 1
        do i = 1, n
            if (h(k) <= 0.0_dp) goto 70
            k = k + np - i
        end do

100     dff = 0.0_dp
110     fa = f
        isfv = 1
        xa = x
        ga = g

130     itr = itr + 1
        ial = 0
        if (itr > niter_limit) then
            status = n1qn1_max_iterations
            goto 250
        end if

        if (present(progress)) then
            stop_requested = .false.
            if (present(user_data)) then
                call progress(itr, nfun, x, f, g, stop_requested, user_data)
            else
                call progress(itr, nfun, x, f, g, stop_requested)
            end if
            if (stop_requested) then
                status = n1qn1_user_stop
                goto 250
            end if
        end if

        if (iprint >= 2) write(output_unit, '(a,i0,a,i0,a,es14.6)') &
            "n1qn1 iteration ", itr, ", evaluations ", nfun, ", f = ", fa

        d = -ga
        w(1) = d(1)
        if (n > 1) goto 400
        d(1) = d(1) / h(1)
        goto 412

400     do i = 2, n
            ij = i
            i1 = i - 1
            v = d(i)
            do j = 1, i1
                v = v - h(ij) * d(j)
                ij = ij + n - j
            end do
            w(i) = v
            d(i) = v
        end do
        d(n) = d(n) / h(ij)
        np = n + 1
        do nip = 2, n
            i = np - nip
            ii = ij - nip
            v = d(i) / h(ii)
            ip = i + 1
            ij = ii
            do j = ip, n
                ii = ii + 1
                v = v - h(ii) * d(j)
            end do
            d(i) = v
        end do

412     c = 0.0_dp
        dga = 0.0_dp
        do i = 1, n
            c = max(c, abs(d(i) / scale(i)))
            dga = dga + ga(i) * d(i)
        end do
        if (dga >= 0.0_dp) then
            if (maxval(abs(ga)) <= sqrt(epsilon(1.0_dp))) then
                status = n1qn1_success
            else
                status = n1qn1_not_descent
            end if
            goto 240
        end if

        stmin = 0.0_dp
        stepbd = 0.0_dp
        steplb = acc / c
        fmin = fa
        gmin = dga
        step = 1.0_dp
        if (dff <= 0.0_dp) step = min(step, 1.0_dp / c)
        if (dff > 0.0_dp) step = min(step, (dff + dff) / (-dga))

170     c = stmin + step
        if (nfun >= nsim_limit) then
            status = n1qn1_max_evaluations
            goto 250
        end if
        nfun = nfun + 1
        xb = xa + c * d
        indic = 4
        if (present(user_data)) then
            fb = objective(xb, user_data)
            call gradient_function(xb, gb, user_data)
        else
            fb = objective(xb)
            call gradient_function(xb, gb)
        end if
        if (.not. ieee_is_finite(fb) .or. any(.not. ieee_is_finite(gb))) indic = -1

        if (indic > 0) goto 185
        if (indic < 0) goto 183
        x = xb
        g = gb
        status = n1qn1_user_stop
        goto 250

183     stepbd = step
        ial = 1
        step = step / 10.0_dp
        if (iprint >= 3) write(output_unit, '(a,es14.6)') &
            "n1qn1 rejected non-finite step ", c
        if (stepbd > steplb) goto 170
        status = n1qn1_nonfinite
        goto 240

185     isfv = min(2, isfv)
        if (fb > f) goto 220
        if (fb < f) goto 200
        gl1 = sum((scale * g) ** 2)
        gl2 = sum((scale * gb) ** 2)
        if (gl2 >= gl1) goto 220

200     isfv = 3
        f = fb
        x = xb
        g = gb

220     dgb = dot_product(gb, d)
        if (iprint >= 3) write(output_unit, '(a,es14.6,a,es14.6,a,es14.6)') &
            "n1qn1 step = ", c, ", df = ", fb - fa, ", derivative = ", dgb

        if (fb - fa <= 0.1_dp * c * dga) goto 280
        ial = 0
        if (step > steplb) goto 270

240     if (isfv >= 2) goto 110
        goto 250

250     acc = dot_product(g, g)
        niter = itr
        nsim = nfun
        return

270     stepbd = step
        c = gmin + dgb - 3.0_dp * (fb - fmin) / step
        if (c == 0.0_dp) goto 250
        cc = abs(c) - gmin * (dgb / abs(c))
        cc = sqrt(abs(c)) * sqrt(max(0.0_dp, cc))
        c = (c - gmin + cc) / (dgb - gmin + cc + cc)
        step = step * max(0.1_dp, c)
        goto 170

280     if (ial == 0) goto 285
        if (stepbd > steplb) goto 285
        goto 240

285     stepbd = stepbd - step
        stmin = c
        fmin = fb
        gmin = dgb
        step = 9.0_dp * stmin
        if (stepbd > 0.0_dp) step = 0.5_dp * stepbd
        c = dga + 3.0_dp * dgb - 4.0_dp * (fb - fa) / stmin
        if (c > 0.0_dp) step = min(step, stmin * max(1.0_dp, -dgb / c))
        if (dgb < 0.7_dp * dga) goto 170

        isfv = 4 - isfv
        if (stmin + step <= steplb) goto 240

        ir = -n
        xa = xb
        xb = ga
        d = gb - ga
        ga = gb
        call majour(h, xb, w, n, 1.0_dp / dga, ir, 1, 0.0_dp)
        ir = -ir
        call majour(h, d, d, n, 1.0_dp / (stmin * (dgb - dga)), ir, 1, 0.0_dp)
        if (ir < n) then
            status = n1qn1_rank_loss
            goto 250
        end if
        dff = fa - fb
        fa = fb
        goto 130
    end subroutine n1qn1_kernel

    subroutine majour(hm, hd, dd, n, hno, ir, indic, eps_update)
        integer, intent(in) :: n, indic
        integer, intent(inout) :: ir
        real(dp), intent(inout) :: hm(:), hd(n), dd(n)
        real(dp), intent(in) :: hno, eps_update

        real(dp) :: b, r, y, gm, del, hml, hon, honm
        integer :: i, j, ll, mm, np, iplus

        if (n == 1) goto 100
        np = n + 1
        if (hno > 0.0_dp) goto 99
        if (hno == 0.0_dp .or. ir == 0) goto 999
        hon = 1.0_dp / hno
        ll = 1
        if (indic == 0) goto 1

        do i = 1, n
            if (hm(ll) == 0.0_dp) goto 2
            hon = hon + dd(i) * dd(i) / hm(ll)
            ll = ll + np - i
        end do
2       goto 3

1       dd = hd
        do i = 1, n
            iplus = i + 1
            del = dd(i)
            if (hm(ll) <= 0.0_dp) then
                dd(i) = 0.0_dp
                ll = ll + np - i
                cycle
            end if
            hon = hon + del * del / hm(ll)
            if (i /= n) then
                do j = iplus, n
                    ll = ll + 1
                    dd(j) = dd(j) - del * hm(ll)
                end do
            end if
            ll = ll + 1
        end do

3       if (ir <= 0) goto 9
        if (hon > 0.0_dp) goto 10
        if (indic <= 1) goto 99
        goto 11

9       hon = 0.0_dp
        ir = -ir - 1
        goto 11

10      hon = eps_update / hno
        if (eps_update == 0.0_dp) ir = ir - 1

11      mm = 1
        honm = hon
        do i = 1, n
            j = np - i
            ll = ll - i
            if (hm(ll) /= 0.0_dp) honm = hon - dd(j) * dd(j) / hm(ll)
            dd(j) = hon
            hon = honm
        end do
        goto 13

99      mm = 0
        honm = 1.0_dp / hno

13      ll = 1
        do i = 1, n
            iplus = i + 1
            del = hd(i)
            if (hm(ll) <= 0.0_dp) then
                if (ir <= 0 .and. hno >= 0.0_dp .and. del /= 0.0_dp) then
                    ir = 1 - ir
                    hm(ll) = del * del / honm
                    if (i == n) goto 999
                    do j = iplus, n
                        ll = ll + 1
                        hm(ll) = hd(j) / del
                    end do
                    goto 999
                end if
                hon = honm
                ll = ll + np - i
                cycle
            end if

            hml = del / hm(ll)
            if (mm <= 0) then
                hon = honm + del * hml
            else
                hon = dd(i)
            end if
            r = hon / honm
            hm(ll) = hm(ll) * r
            if (r == 0.0_dp .or. i == n) goto 20
            b = hml / hon
            if (r <= 4.0_dp) then
                do j = iplus, n
                    ll = ll + 1
                    hd(j) = hd(j) - del * hm(ll)
                    hm(ll) = hm(ll) + b * hd(j)
                end do
            else
                gm = honm / hon
                do j = iplus, n
                    ll = ll + 1
                    y = hm(ll)
                    hm(ll) = b * hd(j) + y * gm
                    hd(j) = hd(j) - del * y
                end do
            end if
            honm = hon
            ll = ll + 1
        end do

20      if (ir < 0) ir = -ir
        goto 999

100     hm(1) = hm(1) + hno * hd(1) * hd(1)
        ir = 1
        if (hm(1) <= 0.0_dp) then
            hm(1) = 0.0_dp
            ir = 0
        end if
999     continue
    end subroutine majour

    subroutine pack_lower(matrix, packed)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), intent(out) :: packed(:)
        integer :: n, i, j, k

        n = size(matrix, 1)
        if (size(matrix, 2) /= n .or. size(packed) < n1qn1_packed_size(n)) error stop &
            "pack_lower: inconsistent dimensions"
        k = 0
        do j = 1, n
            do i = j, n
                k = k + 1
                packed(k) = matrix(i, j)
            end do
        end do
    end subroutine pack_lower

    subroutine unpack_lower(packed, matrix)
        real(dp), intent(in) :: packed(:)
        real(dp), intent(out) :: matrix(:, :)
        integer :: n, i, j, k

        n = size(matrix, 1)
        if (size(matrix, 2) /= n .or. size(packed) < n1qn1_packed_size(n)) error stop &
            "unpack_lower: inconsistent dimensions"
        matrix = 0.0_dp
        k = 0
        do j = 1, n
            do i = j, n
                k = k + 1
                matrix(i, j) = packed(k)
                matrix(j, i) = packed(k)
            end do
        end do
    end subroutine unpack_lower

    subroutine factor_to_hessian(factor, hessian)
        real(dp), intent(in) :: factor(:)
        real(dp), intent(out) :: hessian(:, :)
        real(dp), allocatable :: l(:, :), dmat(:, :)
        integer :: n, i, j, k

        n = size(hessian, 1)
        if (size(hessian, 2) /= n .or. size(factor) < n1qn1_packed_size(n)) error stop &
            "factor_to_hessian: inconsistent dimensions"
        allocate(l(n, n), dmat(n, n))
        l = 0.0_dp
        dmat = 0.0_dp
        do i = 1, n
            l(i, i) = 1.0_dp
        end do
        k = 0
        do j = 1, n
            do i = j, n
                k = k + 1
                if (i == j) then
                    dmat(i, i) = factor(k)
                else
                    l(i, j) = factor(k)
                end if
            end do
        end do
        hessian = matmul(matmul(l, dmat), transpose(l))
    end subroutine factor_to_hessian

    function n1qn1_status_message(status) result(message)
        integer, intent(in) :: status
        character(len=:), allocatable :: message

        select case (status)
        case (n1qn1_success)
            message = "converged"
        case (n1qn1_max_iterations)
            message = "maximum iteration count reached"
        case (n1qn1_max_evaluations)
            message = "maximum objective/gradient evaluation count reached"
        case (n1qn1_not_descent)
            message = "the quasi-Newton direction was not a descent direction"
        case (n1qn1_rank_loss)
            message = "the updated Hessian factor lost rank"
        case (n1qn1_user_stop)
            message = "stopped by the user callback"
        case (n1qn1_nonfinite)
            message = "objective or gradient returned a non-finite value"
        case default
            message = "invalid input"
        end select
    end function n1qn1_status_message

    subroutine set_error(result, status, message)
        type(n1qn1_result_t), intent(out) :: result
        integer, intent(in) :: status
        character(len=*), intent(in) :: message
        result%status = status
        result%message = message
    end subroutine set_error

end module n1qn1_module
