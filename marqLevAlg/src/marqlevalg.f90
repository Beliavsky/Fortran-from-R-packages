! SPDX-License-Identifier: GPL-2.0-or-later
module marqlevalg
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use mla_kinds, only : dp
    use mla_interfaces, only : objective_fn, gradient_fn, hessian_fn
    use mla_linalg, only : packed_size, pack_upper, unpack_upper, invert_packed_spd, solve_packed_spd
    use mla_derivatives, only : numerical_derivatives, numerical_information_from_gradient, &
        numerical_derivatives_scaled, numerical_information_from_gradient_scaled
    implicit none
    private

    type, public :: mla_control
        integer :: maxiter = 500
        real(dp) :: epsa = 1.0e-4_dp
        real(dp) :: epsb = 1.0e-4_dp
        real(dp) :: epsd = 1.0e-4_dp
        logical :: minimize = .true.
        logical :: blinding = .true.
        integer :: multiple_try = 25
        integer, allocatable :: partial_h(:)
    end type mla_control

    type, public :: mla_result
        real(dp), allocatable :: par(:)
        real(dp), allocatable :: grad(:)
        real(dp), allocatable :: vcov(:, :)
        real(dp) :: fn_value = 0.0_dp
        real(dp) :: ca = huge(1.0_dp)
        real(dp) :: cb = huge(1.0_dp)
        real(dp) :: rdm = huge(1.0_dp)
        integer :: iterations = 0
        integer :: ier = 0
        integer :: istop = 0
    end type mla_result

    interface marqlev_optimize
        module procedure marqlev_numeric
        module procedure marqlev_with_gradient
        module procedure marqlev_with_hessian
    end interface marqlev_optimize

    public :: dp, marqlev_optimize
    public :: deriva, deriva_grad

contains

    function dispatch_objective(x, fn) result(f)
        real(dp), intent(in) :: x(:)
        procedure(objective_fn) :: fn
        real(dp) :: f
        f = fn(x)
    end function dispatch_objective

    subroutine dispatch_gradient(x, g, gr)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        procedure(gradient_fn) :: gr
        call gr(x, g)
    end subroutine dispatch_gradient

    subroutine dispatch_hessian(x, h, hess)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: h(:, :)
        procedure(hessian_fn) :: hess
        call hess(x, h)
    end subroutine dispatch_hessian

    subroutine deriva(x, fn, f, info, grad)
        real(dp), intent(in) :: x(:)
        procedure(objective_fn) :: fn
        real(dp), intent(out) :: f
        real(dp), intent(out) :: info(:, :), grad(:)
        call numerical_derivatives(x, fn, f, info, grad)
    end subroutine deriva

    subroutine deriva_grad(x, gr, info)
        real(dp), intent(in) :: x(:)
        procedure(gradient_fn) :: gr
        real(dp), intent(out) :: info(:, :)
        call numerical_information_from_gradient(x, gr, info)
    end subroutine deriva_grad

    subroutine marqlev_numeric(x0, fn, result, control)
        real(dp), intent(in) :: x0(:)
        procedure(objective_fn) :: fn
        type(mla_result), intent(out) :: result
        type(mla_control), intent(in), optional :: control
        type(mla_control) :: ctl

        ctl = mla_control()
        if (present(control)) ctl = control
        call optimize_core(x0, fn, result, ctl, deriv_mode=0)
    end subroutine marqlev_numeric

    subroutine marqlev_with_gradient(x0, fn, gr, result, control)
        real(dp), intent(in) :: x0(:)
        procedure(objective_fn) :: fn
        procedure(gradient_fn) :: gr
        type(mla_result), intent(out) :: result
        type(mla_control), intent(in), optional :: control
        type(mla_control) :: ctl

        ctl = mla_control()
        if (present(control)) ctl = control
        call optimize_core(x0, fn, result, ctl, deriv_mode=1, gr=gr)
    end subroutine marqlev_with_gradient

    subroutine marqlev_with_hessian(x0, fn, gr, hess, result, control)
        real(dp), intent(in) :: x0(:)
        procedure(objective_fn) :: fn
        procedure(gradient_fn) :: gr
        procedure(hessian_fn) :: hess
        type(mla_result), intent(out) :: result
        type(mla_control), intent(in), optional :: control
        type(mla_control) :: ctl

        ctl = mla_control()
        if (present(control)) ctl = control
        call optimize_core(x0, fn, result, ctl, deriv_mode=2, gr=gr, hess=hess)
    end subroutine marqlev_with_hessian

    subroutine optimize_core(x0, fn, result, ctl, deriv_mode, gr, hess)
        real(dp), intent(in) :: x0(:)
        procedure(objective_fn) :: fn
        type(mla_result), intent(out) :: result
        type(mla_control), intent(in) :: ctl
        integer, intent(in) :: deriv_mode
        procedure(gradient_fn), optional :: gr
        procedure(hessian_fn), optional :: hess

        real(dp), allocatable :: x(:), oldx(:), g(:), info(:, :), info_work(:, :)
        real(dp), allocatable :: pinfo(:), pinv(:), delta(:), trial(:), invmat(:, :)
        real(dp) :: f, internal_f, internal_trial, ca, cb, dd, olddd
        real(dp) :: da, dm, ga, tr, th, eps, ep, maxd, vw, step
        integer :: n, np, i, ier, idpos, ncount, tries, ni
        logical :: converged, line_ok, finite_trial

        n = size(x0)
        np = packed_size(n)
        allocate(x(n), oldx(n), g(n), info(n, n), info_work(n, n))
        allocate(pinfo(np), pinv(np), delta(n), trial(n), invmat(n, n))
        x = x0
        oldx = x
        ca = ctl%epsa + 1.0_dp
        cb = ctl%epsb + 1.0_dp
        dd = ctl%epsd + 1.0_dp
        olddd = 1.0_dp
        da = 1.0e-2_dp
        dm = 5.0_dp
        th = 1.0e-5_dp
        eps = 1.0e-7_dp
        ep = 1.0e-20_dp
        ni = 0
        ier = 0
        result%istop = 0
        converged = .false.

        do
            if (.not. all(ieee_is_finite(x))) then
                result%istop = 4
                exit
            end if

            call evaluate_derivatives(x, fn, ctl%minimize, deriv_mode, internal_f, info, g, gr, hess)
            if (deriv_mode == 0 .and. ni == 0 .and. ctl%multiple_try > 1) then
                tries = 0
                do while (.not. ieee_is_finite(internal_f) .and. tries < ctl%multiple_try)
                    tries = tries + 1
                    x = x / 2.0_dp
                    call evaluate_derivatives(x, fn, ctl%minimize, deriv_mode, internal_f, info, g, gr, hess)
                end do
            end if

            if (.not. ieee_is_finite(internal_f) .or. .not. all(ieee_is_finite(g)) .or. &
                .not. all(ieee_is_finite(info))) then
                result%istop = 4
                exit
            end if

            f = merge(-internal_f, internal_f, ctl%minimize)
            call pack_upper(info, pinfo)
            call invert_packed_spd(pinfo, pinv, ier)
            if (ier == -1) then
                dd = ctl%epsd + 1.0_dp
            else
                call unpack_upper(pinv, invmat)
                dd = dot_product(g, matmul(invmat, g)) / real(n, dp)
                if (.not. ieee_is_finite(dd)) dd = ctl%epsd + 1.0_dp
            end if

            oldx = x
            if (dd <= olddd) olddd = dd
            if (ca < ctl%epsa .and. cb < ctl%epsb .and. dd < ctl%epsd) then
                converged = .true.
                result%istop = 1
                exit
            end if

            if (ca < ctl%epsa .and. cb < ctl%epsb .and. ier == -1 .and. allocated(ctl%partial_h)) then
                if (size(ctl%partial_h) > 0) then
                    call partial_rdm(info, g, ctl%partial_h, dd, ier, invmat)
                    if (ier == 0 .and. dd < ctl%epsd) then
                        converged = .true.
                        result%istop = 3
                        exit
                    end if
                end if
            end if

            tr = 0.0_dp
            do i = 1, n
                tr = tr + abs(info(i, i))
            end do
            tr = tr / real(n, dp)
            ga = 0.01_dp
            ncount = 0

            do
                info_work = info
                do i = 1, n
                    if (abs(info(i, i)) > 0.0_dp) then
                        info_work(i, i) = info(i, i) + da * &
                            ((1.0_dp - ga) * abs(info(i, i)) + ga * tr)
                    else
                        info_work(i, i) = da * ga * tr
                    end if
                end do
                call pack_upper(info_work, pinfo)
                call solve_packed_spd(pinfo, g, delta, idpos)
                if (idpos == 0) exit

                ncount = ncount + 1
                if (ncount <= 3 .or. ga >= 1.0_dp) then
                    da = da * dm
                else
                    ga = min(1.0_dp, ga * dm)
                end if
                if (ncount >= 10) exit
            end do

            if (idpos /= 0) then
                result%istop = 4
                exit
            end if

            trial = x + delta
            internal_trial = internal_objective(trial, fn, ctl%minimize)
            finite_trial = ieee_is_finite(internal_trial)
            if (.not. finite_trial) then
                if (ctl%blinding) then
                    internal_trial = -500000.0_dp
                else
                    result%istop = 4
                    exit
                end if
            end if

            if (internal_f < internal_trial) then
                if (da < eps) then
                    da = eps
                else
                    da = da / (dm + 2.0_dp)
                end if
            else
                maxd = maxval(abs(delta))
                if (maxd <= tiny(1.0_dp)) then
                    vw = th
                else
                    vw = th / maxd
                end if
                step = log(1.5_dp)
                call search_step(vw, step, x, delta, fn, ctl%minimize, internal_trial, line_ok)
                if (.not. line_ok) then
                    result%istop = 4
                    exit
                end if
                delta = vw * delta
                da = (dm - 3.0_dp) * da
            end if

            trial = x + delta
            internal_trial = internal_objective(trial, fn, ctl%minimize)
            if (.not. ieee_is_finite(internal_trial)) then
                result%istop = 4
                exit
            end if
            ca = dot_product(delta, delta)
            cb = abs(internal_f - internal_trial)
            x = trial
            ni = ni + 1
            if (ni >= ctl%maxiter) then
                result%istop = 2
                exit
            end if
        end do

        call evaluate_derivatives(x, fn, ctl%minimize, deriv_mode, internal_f, info, g, gr, hess)
        f = merge(-internal_f, internal_f, ctl%minimize)
        call pack_upper(info, pinfo)
        call invert_packed_spd(pinfo, pinv, ier)
        if (ier == 0) then
            call unpack_upper(pinv, invmat)
        else
            invmat = info
        end if

        allocate(result%par(n), result%grad(n), result%vcov(n, n))
        result%par = x
        result%fn_value = f
        result%grad = merge(-g, g, ctl%minimize)
        result%vcov = invmat
        result%iterations = ni
        result%ier = ier
        result%ca = ca
        result%cb = cb
        result%rdm = dd
        if (result%istop == 0) result%istop = merge(1, 2, converged)
    end subroutine optimize_core

    subroutine evaluate_derivatives(x, fn, minimize, mode, f, info, g, gr, hess)
        real(dp), intent(in) :: x(:)
        procedure(objective_fn) :: fn
        logical, intent(in) :: minimize
        integer, intent(in) :: mode
        real(dp), intent(out) :: f, info(:, :), g(:)
        procedure(gradient_fn), optional :: gr
        procedure(hessian_fn), optional :: hess
        real(dp), allocatable :: userg(:), userh(:, :)
        real(dp) :: scale

        allocate(userg(size(x)), userh(size(x), size(x)))
        scale = merge(-1.0_dp, 1.0_dp, minimize)
        select case (mode)
        case (0)
            call numerical_derivatives_scaled(x, fn, scale, f, info, g)
        case (1)
            f = scale * dispatch_objective(x, fn)
            call dispatch_gradient(x, userg, gr)
            g = scale * userg
            call numerical_information_from_gradient_scaled(x, gr, scale, info)
        case (2)
            f = scale * dispatch_objective(x, fn)
            call dispatch_gradient(x, userg, gr)
            g = scale * userg
            call dispatch_hessian(x, userh, hess)
            info = userh
        end select
    end subroutine evaluate_derivatives

    function internal_objective(x, fn, minimize) result(f)
        real(dp), intent(in) :: x(:)
        procedure(objective_fn) :: fn
        logical, intent(in) :: minimize
        real(dp) :: f

        f = dispatch_objective(x, fn)
        if (minimize) f = -f
    end function internal_objective

    subroutine partial_rdm(info, g, drop, rdm, ier, invfull)
        real(dp), intent(in) :: info(:, :), g(:)
        integer, intent(in) :: drop(:)
        real(dp), intent(out) :: rdm
        integer, intent(out) :: ier
        real(dp), intent(out) :: invfull(:, :)
        integer, allocatable :: keep(:)
        real(dp), allocatable :: sub(:, :), invsub(:, :), p(:), pi(:), gs(:)
        integer :: i, j, n, nr
        logical :: is_drop

        n = size(g)
        allocate(keep(n))
        nr = 0
        do i = 1, n
            is_drop = any(drop == i)
            if (.not. is_drop) then
                nr = nr + 1
                keep(nr) = i
            end if
        end do
        if (nr <= 0) then
            ier = -1
            rdm = huge(1.0_dp)
            invfull = 0.0_dp
            return
        end if
        allocate(sub(nr, nr), invsub(nr, nr), p(nr * (nr + 1) / 2), &
                 pi(nr * (nr + 1) / 2), gs(nr))
        do i = 1, nr
            gs(i) = g(keep(i))
            do j = 1, nr
                sub(i, j) = info(keep(i), keep(j))
            end do
        end do
        call pack_upper(sub, p)
        call invert_packed_spd(p, pi, ier)
        invfull = 0.0_dp
        if (ier == 0) then
            call unpack_upper(pi, invsub)
            rdm = dot_product(gs, matmul(invsub, gs)) / real(nr, dp)
            do i = 1, nr
                do j = 1, nr
                    invfull(keep(i), keep(j)) = invsub(i, j)
                end do
            end do
        else
            rdm = huge(1.0_dp)
        end if
    end subroutine partial_rdm

    function line_value_dispatch(t, x, delta, fn, minimize) result(v)
        real(dp), intent(in) :: t, x(:), delta(:)
        procedure(objective_fn) :: fn
        logical, intent(in) :: minimize
        real(dp) :: v
        real(dp) :: z(size(x))

        z = x + exp(t) * delta
        v = -internal_objective(z, fn, minimize)
    end function line_value_dispatch

    subroutine search_step(vw, step_in, x, delta, fn, minimize, internal_best, ok)
        real(dp), intent(inout) :: vw
        real(dp), intent(in) :: step_in
        real(dp), intent(in) :: x(:), delta(:)
        procedure(objective_fn) :: fn
        logical, intent(in) :: minimize
        real(dp), intent(out) :: internal_best
        logical, intent(out) :: ok
        real(dp) :: step, vlw1, vlw2, vlw3, fi1, fi2, fi3, vm, fim, denom
        integer :: iter

        ok = .false.
        if (.not. ieee_is_finite(vw) .or. vw <= 0.0_dp) return
        step = step_in
        vlw1 = log(vw)
        vlw2 = vlw1 + step
        fi1 = line_value_dispatch(vlw1, x, delta, fn, minimize)
        fi2 = line_value_dispatch(vlw2, x, delta, fn, minimize)
        if (.not. ieee_is_finite(fi1) .or. .not. ieee_is_finite(fi2)) return

        if (fi2 >= fi1) then
            vlw3 = vlw2
            vlw2 = vlw1
            fi3 = fi2
            fi2 = fi1
            step = -step
            vlw1 = vlw2 + step
            fi1 = line_value_dispatch(vlw1, x, delta, fn, minimize)
            denom = 2.0_dp * (fi1 - 2.0_dp * fi2 + fi3)
            if (abs(denom) > tiny(1.0_dp)) then
                vm = vlw2 - step * (fi1 - fi3) / denom
                fim = line_value_dispatch(vm, x, delta, fn, minimize)
            else
                vm = vlw2
                fim = fi2
            end if
            if (.not. ieee_is_finite(fim)) fim = huge(1.0_dp)
            if (fim <= fi2) then
                vw = exp(vm)
            else
                fim = fi2
                vw = exp(vlw2)
            end if
        else
            call swap(vlw1, vlw2)
            call swap(fi1, fi2)
            fim = fi2
            vm = vlw2
            do iter = 1, 40
                vlw3 = vlw2
                vlw2 = vlw1
                fi3 = fi2
                fi2 = fi1
                vlw1 = vlw2 + step
                fi1 = line_value_dispatch(vlw1, x, delta, fn, minimize)
                if (.not. ieee_is_finite(fi1)) exit
                if (fi1 > fi2) then
                    denom = 2.0_dp * (fi1 - 2.0_dp * fi2 + fi3)
                    if (abs(denom) > tiny(1.0_dp)) then
                        vm = vlw2 - step * (fi1 - fi3) / denom
                        fim = line_value_dispatch(vm, x, delta, fn, minimize)
                    else
                        vm = vlw2
                        fim = fi2
                    end if
                    if (.not. ieee_is_finite(fim) .or. fim > fi2) then
                        vm = vlw2
                        fim = fi2
                    end if
                    exit
                end if
                if (abs(fi1 - fi2) <= epsilon(1.0_dp) * max(1.0_dp, abs(fi2))) then
                    vm = vlw2
                    fim = fi2
                    exit
                end if
                vm = vlw1
                fim = fi1
            end do
            vw = exp(vm)
        end if

        internal_best = -fim
        ok = ieee_is_finite(internal_best) .and. ieee_is_finite(vw)
    contains
        subroutine swap(a, b)
            real(dp), intent(inout) :: a, b
            real(dp) :: t
            t = a
            a = b
            b = t
        end subroutine swap
    end subroutine search_step

end module marqlevalg
