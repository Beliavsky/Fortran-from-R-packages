! SPDX-License-Identifier: GPL-3.0-only
!
! Modern Fortran translation of the computational code in the R package
! subplex 1.9. The original Subplex algorithm was coded by Tom Rowan.
! See original/subplex-master and LICENSE for provenance and licensing.
module subplex
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

    abstract interface
        function subplex_objective(x) result(f)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: f
        end function subplex_objective
    end interface

    type, public :: subplex_control
        real(dp) :: reltol = epsilon(1.0_dp)
        integer :: maxeval = 10000
        real(dp), allocatable :: parscale(:)
    end type subplex_control

    type, public :: subplex_result
        real(dp), allocatable :: x(:)
        real(dp) :: value = huge(1.0_dp)
        integer :: counts = 0
        integer :: convergence = -2
        character(len=:), allocatable :: message
        real(dp), allocatable :: hessian(:,:)
    end type subplex_result

    type :: subplex_state
        real(dp) :: alpha = 1.0_dp
        real(dp) :: beta = 0.5_dp
        real(dp) :: gamma = 2.0_dp
        real(dp) :: delta = 0.5_dp
        real(dp) :: psi = 0.25_dp
        real(dp) :: omega = 0.1_dp
        integer :: nsmin = 1
        integer :: nsmax = 1
    end type subplex_state

    public :: subplex_minimize, numerical_hessian

contains

    subroutine subplex_minimize(fn, x0, result, control, compute_hessian)
        procedure(subplex_objective) :: fn
        real(dp), intent(in) :: x0(:)
        type(subplex_result), intent(out) :: result
        type(subplex_control), intent(in), optional :: control
        logical, intent(in), optional :: compute_hessian

        type(subplex_control) :: ctl
        real(dp), allocatable :: scale(:), hstep(:)
        logical :: want_hessian
        integer :: n

        n = size(x0)
        allocate(result%x(n))
        result%x = x0
        result%value = huge(1.0_dp)
        result%counts = 0
        result%convergence = -2
        result%message = "invalid input"

        if (n < 1) return

        ctl%reltol = epsilon(1.0_dp)
        ctl%maxeval = 10000
        if (present(control)) then
            ctl%reltol = control%reltol
            ctl%maxeval = control%maxeval
            if (allocated(control%parscale)) then
                allocate(ctl%parscale(size(control%parscale)))
                ctl%parscale = control%parscale
            end if
        end if

        if (ctl%reltol < 0.0_dp .or. ctl%maxeval <= 0) return

        if (allocated(ctl%parscale)) then
            if (size(ctl%parscale) /= 1 .and. size(ctl%parscale) /= n) return
            if (size(ctl%parscale) == 1) then
                allocate(scale(n))
                scale = abs(ctl%parscale(1))
            else
                allocate(scale(n))
                scale = abs(ctl%parscale)
            end if
        else
            allocate(scale(n))
            scale = 1.0_dp
        end if

        call subplx_core(fn, ctl%reltol, ctl%maxeval, scale, result%x, &
            result%value, result%counts, result%convergence)

        select case (result%convergence)
        case (-2)
            result%message = "parscale is too small relative to par"
        case (-1)
            result%message = "number of function evaluations exceeds maxeval"
        case (0)
            result%message = "success! tolerance satisfied"
        case (1)
            result%message = "limit of machine precision reached"
        case default
            result%message = "unexpected Subplex convergence code"
        end select

        want_hessian = .false.
        if (present(compute_hessian)) want_hessian = compute_hessian
        if (want_hessian .and. result%convergence /= -2) then
            allocate(hstep(n), result%hessian(n,n))
            hstep = abs(scale) * epsilon(1.0_dp)**(1.0_dp/3.0_dp)
            call numerical_hessian(fn, result%x, hstep, result%hessian)
        end if
    end subroutine subplex_minimize

    subroutine subplx_core(fn, tol, maxnfe, scale, x, fx, nfe, iflag)
        procedure(subplex_objective) :: fn
        real(dp), intent(in) :: tol
        integer, intent(in) :: maxnfe
        real(dp), intent(in) :: scale(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(out) :: fx
        integer, intent(out) :: nfe, iflag

        type(subplex_state) :: state
        integer :: n, nsubs, ins, ipptr, ns, i
        integer, allocatable :: ip(:), nsvals(:)
        real(dp), allocatable :: step(:), oldx(:), absdx(:)
        real(dp) :: sfx

        n = size(x)
        if (size(scale) /= n) then
            iflag = -2
            fx = huge(1.0_dp)
            nfe = 0
            return
        end if

        do i = 1, n
            if (same_real(x(i) + scale(i), x(i))) then
                iflag = -2
                fx = huge(1.0_dp)
                nfe = 0
                return
            end if
        end do

        call set_default_state(n, state)
        allocate(step(n), oldx(n), absdx(n), ip(n), nsvals(n))
        step = scale
        ip = [(i, i=1,n)]
        nfe = 0
        call eval_objective(fn, x, sfx, nfe)

        do
            absdx = abs(step)
            call sort_descending_indices(absdx, ip)
            call partition_subspaces(state, ip, absdx, nsubs, nsvals)
            oldx = x
            ipptr = 1

            do ins = 1, nsubs
                ns = nsvals(ins)
                call simplex_minimize(fn, state, step, ip(ipptr:ipptr+ns-1), &
                    maxnfe, x, sfx, nfe, iflag)
                if (iflag /= 0) then
                    fx = sfx
                    return
                end if
                ipptr = ipptr + ns
            end do

            absdx = x - oldx
            do i = 1, n
                if (max(abs(absdx(i)), abs(step(i))*state%psi) / &
                    max(abs(x(i)), 1.0_dp) > tol) exit
            end do
            if (i > n) then
                iflag = 0
                fx = sfx
                return
            end if
            call set_step(state, nsubs, absdx, step)
        end do
    end subroutine subplx_core

    subroutine set_default_state(n, state)
        integer, intent(in) :: n
        type(subplex_state), intent(out) :: state

        state%alpha = 1.0_dp
        state%beta = 0.5_dp
        state%gamma = 2.0_dp
        state%delta = 0.5_dp
        state%psi = 0.25_dp
        state%omega = 0.1_dp
        state%nsmin = min(2, n)
        state%nsmax = min(5, n)
    end subroutine set_default_state

    subroutine eval_objective(fn, x, fx, nfe)
        procedure(subplex_objective) :: fn
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fx
        integer, intent(inout) :: nfe

        fx = fn(x)
        nfe = nfe + 1
    end subroutine eval_objective

    subroutine eval_subspace(fn, ips, xs, x, fx, nfe)
        procedure(subplex_objective) :: fn
        integer, intent(in) :: ips(:)
        real(dp), intent(in) :: xs(:)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(out) :: fx
        integer, intent(inout) :: nfe

        x(ips) = xs
        call eval_objective(fn, x, fx, nfe)
    end subroutine eval_subspace

    subroutine simplex_minimize(fn, state, step, ips, maxnfe, x, fx, nfe, iflag)
        procedure(subplex_objective) :: fn
        type(subplex_state), intent(in) :: state
        real(dp), intent(in) :: step(:)
        integer, intent(in) :: ips(:)
        integer, intent(in) :: maxnfe
        real(dp), intent(inout) :: x(:)
        real(dp), intent(inout) :: fx
        integer, intent(inout) :: nfe
        integer, intent(out) :: iflag

        integer :: ns, npts, icent, itemp, ih, il, inew, is, j
        real(dp), allocatable :: s(:,:), fs(:)
        real(dp) :: fc, fe, fr, stol
        logical :: small, updatc

        ns = size(ips)
        npts = ns + 1
        icent = ns + 2
        itemp = ns + 3
        allocate(s(ns,ns+3), fs(ns+1))

        updatc = .false.
        call start_simplex(x, step, ips, s, small)
        if (small) then
            iflag = 1
            return
        end if

        fs(1) = fx
        do j = 2, npts
            call eval_subspace(fn, ips, s(:,j), x, fs(j), nfe)
        end do
        il = 1
        call order_simplex(fs, il, is, ih)
        stol = state%psi * stable_distance(s(:,ih), s(:,il))
        inew = ih

        do
            call simplex_centroid(s, ih, inew, updatc, s(:,icent))
            updatc = .true.
            inew = ih

            call new_point(state%alpha, s(:,icent), s(:,ih), s(:,itemp), small)
            if (.not. small) then
                call eval_subspace(fn, ips, s(:,itemp), x, fr, nfe)
                if (fr < fs(il)) then
                    call new_point(-state%gamma, s(:,icent), s(:,itemp), s(:,ih), small)
                    if (.not. small) then
                        call eval_subspace(fn, ips, s(:,ih), x, fe, nfe)
                        if (fe < fr) then
                            fs(ih) = fe
                        else
                            s(:,ih) = s(:,itemp)
                            fs(ih) = fr
                        end if
                    end if
                else if (fr < fs(is)) then
                    s(:,ih) = s(:,itemp)
                    fs(ih) = fr
                else
                    if (fr > fs(ih)) then
                        call new_point(-state%beta, s(:,icent), s(:,ih), s(:,itemp), small)
                    else
                        call new_point_inplace(-state%beta, s(:,icent), s(:,itemp), small)
                    end if
                    if (.not. small) then
                        call eval_subspace(fn, ips, s(:,itemp), x, fc, nfe)
                        if (fc < min(fr, fs(ih))) then
                            s(:,ih) = s(:,itemp)
                            fs(ih) = fc
                        else
                            do j = 1, npts
                                if (j /= il) then
                                    call new_point_inplace(-state%delta, s(:,il), s(:,j), small)
                                    if (small) exit
                                    call eval_subspace(fn, ips, s(:,j), x, fs(j), nfe)
                                end if
                            end do
                        end if
                        updatc = .false.
                    end if
                end if
                if (.not. small) call order_simplex(fs, il, is, ih)
            end if

            fx = fs(il)
            if (nfe >= maxnfe) then
                iflag = -1
                exit
            else if (stable_distance(s(:,ih), s(:,il)) <= stol .or. small) then
                iflag = 0
                exit
            end if
        end do

        x(ips) = s(:,il)
    end subroutine simplex_minimize

    subroutine start_simplex(x, step, ips, s, small)
        real(dp), intent(in) :: x(:), step(:)
        integer, intent(in) :: ips(:)
        real(dp), intent(out) :: s(:,:)
        logical, intent(out) :: small
        integer :: ns, j

        ns = size(ips)
        s(:,1) = x(ips)
        do j = 2, ns + 1
            s(:,j) = s(:,1)
            s(j-1,j) = s(j-1,1) + step(ips(j-1))
        end do
        small = .false.
        do j = 2, ns + 1
            if (same_real(s(j-1,j), s(j-1,1))) then
                small = .true.
                return
            end if
        end do
    end subroutine start_simplex

    subroutine simplex_centroid(s, ih, inew, update, c)
        real(dp), intent(in) :: s(:,:)
        integer, intent(in) :: ih, inew
        logical, intent(in) :: update
        real(dp), intent(inout) :: c(:)
        integer :: j, ns

        ns = size(c)
        if (update) then
            if (ih /= inew) c = c + (s(:,inew) - s(:,ih)) / real(ns, dp)
        else
            c = 0.0_dp
            do j = 1, ns + 1
                if (j /= ih) c = c + s(:,j)
            end do
            c = c / real(ns, dp)
        end if
    end subroutine simplex_centroid

    subroutine new_point(coef, xbase, xold, xnew, small)
        real(dp), intent(in) :: coef
        real(dp), intent(in) :: xbase(:), xold(:)
        real(dp), intent(out) :: xnew(:)
        logical, intent(out) :: small

        xnew = xbase + coef * (xbase - xold)
        small = arrays_same_exact(xnew, xbase) .or. &
            arrays_same_exact(xnew, xold)
    end subroutine new_point

    subroutine new_point_inplace(coef, xbase, xold, small)
        real(dp), intent(in) :: coef
        real(dp), intent(in) :: xbase(:)
        real(dp), intent(inout) :: xold(:)
        logical, intent(out) :: small
        real(dp), allocatable :: prior(:)

        allocate(prior(size(xold)))
        prior = xold
        xold = xbase + coef * (xbase - prior)
        small = arrays_same_exact(xold, xbase) .or. &
            arrays_same_exact(xold, prior)
    end subroutine new_point_inplace

    subroutine order_simplex(fs, il, is, ih)
        real(dp), intent(in) :: fs(:)
        integer, intent(inout) :: il
        integer, intent(out) :: is, ih
        integer :: i, il0, j, npts

        npts = size(fs)
        il0 = il
        j = modulo(il0, npts) + 1
        if (fs(j) >= fs(il)) then
            ih = j
            is = il0
        else
            ih = il0
            is = j
            il = j
        end if
        do i = il0 + 1, il0 + npts - 2
            j = modulo(i, npts) + 1
            if (fs(j) >= fs(ih)) then
                is = ih
                ih = j
            else if (fs(j) > fs(is)) then
                is = j
            else if (fs(j) < fs(il)) then
                il = j
            end if
        end do
    end subroutine order_simplex

    pure real(dp) function stable_distance(x, y) result(d)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: absxmy, scale, ssq
        integer :: i

        absxmy = abs(x(1) - y(1))
        if (absxmy <= 1.0_dp) then
            ssq = absxmy * absxmy
            scale = 1.0_dp
        else
            ssq = 1.0_dp
            scale = absxmy
        end if
        do i = 2, size(x)
            absxmy = abs(x(i) - y(i))
            if (absxmy <= scale) then
                ssq = ssq + (absxmy / scale)**2
            else
                ssq = 1.0_dp + ssq * (scale / absxmy)**2
                scale = absxmy
            end if
        end do
        d = scale * sqrt(ssq)
    end function stable_distance

    subroutine sort_descending_indices(key, ix)
        real(dp), intent(in) :: key(:)
        integer, intent(inout) :: ix(:)
        integer :: i, ifirst, ilast, iswap, ixi, ixip1

        ifirst = 1
        iswap = 1
        ilast = size(ix) - 1
        do while (ifirst <= ilast)
            do i = ifirst, ilast
                ixi = ix(i)
                ixip1 = ix(i+1)
                if (key(ixi) < key(ixip1)) then
                    ix(i) = ixip1
                    ix(i+1) = ixi
                    iswap = i
                end if
            end do
            ilast = iswap - 1
            do i = ilast, ifirst, -1
                ixi = ix(i)
                ixip1 = ix(i+1)
                if (key(ixi) < key(ixip1)) then
                    ix(i) = ixip1
                    ix(i+1) = ixi
                    iswap = i
                end if
            end do
            ifirst = iswap + 1
        end do
    end subroutine sort_descending_indices

    subroutine partition_subspaces(state, ip, absdx, nsubs, nsvals)
        type(subplex_state), intent(in) :: state
        integer, intent(in) :: ip(:)
        real(dp), intent(in) :: absdx(:)
        integer, intent(out) :: nsubs
        integer, intent(out) :: nsvals(:)

        integer :: i, n, nleft, ns1, ns2, nused
        real(dp) :: asleft, as1, as1max, as2, gap, gapmax

        n = size(absdx)
        nsubs = 0
        nused = 0
        nleft = n
        asleft = sum(absdx)
        do while (nused < n)
            nsubs = nsubs + 1
            as1 = 0.0_dp
            do i = 1, state%nsmin - 1
                as1 = as1 + absdx(ip(nused+i))
            end do
            gapmax = -1.0_dp
            as1max = as1
            do ns1 = state%nsmin, min(state%nsmax, nleft)
                as1 = as1 + absdx(ip(nused+ns1))
                ns2 = nleft - ns1
                if (ns2 > 0) then
                    if (ns2 >= ((ns2-1)/state%nsmax + 1)*state%nsmin) then
                        as2 = asleft - as1
                        gap = as1/real(ns1,dp) - as2/real(ns2,dp)
                        if (gap > gapmax) then
                            gapmax = gap
                            nsvals(nsubs) = ns1
                            as1max = as1
                        end if
                    end if
                else
                    if (as1/real(ns1,dp) > gapmax) then
                        nsvals(nsubs) = ns1
                        exit
                    end if
                end if
            end do
            nused = nused + nsvals(nsubs)
            nleft = n - nused
            if (nused < n) asleft = asleft - as1max
        end do
    end subroutine partition_subspaces

    subroutine set_step(state, nsubs, deltax, step)
        type(subplex_state), intent(in) :: state
        integer, intent(in) :: nsubs
        real(dp), intent(in) :: deltax(:)
        real(dp), intent(inout) :: step(:)
        real(dp) :: denom, stpfac
        integer :: i

        if (nsubs > 1) then
            denom = sum(abs(step))
            if (denom > 0.0_dp) then
                stpfac = min(max(sum(abs(deltax))/denom, state%omega), &
                    1.0_dp/state%omega)
            else
                stpfac = state%omega
            end if
        else
            stpfac = state%psi
        end if
        step = stpfac * step
        do i = 1, size(step)
            if (deltax(i) < 0.0_dp .or. deltax(i) > 0.0_dp) then
                step(i) = sign(step(i), deltax(i))
            else
                step(i) = -step(i)
            end if
        end do
    end subroutine set_step

    subroutine numerical_hessian(fn, x, h, hess)
        procedure(subplex_objective) :: fn
        real(dp), intent(in) :: x(:), h(:)
        real(dp), intent(out) :: hess(:,:)
        real(dp), allocatable :: xx(:), df1(:), df2(:)
        integer :: i

        allocate(xx(size(x)), df1(size(x)), df2(size(x)))
        xx = x
        do i = 1, size(x)
            xx(i) = x(i) + h(i)
            call numerical_derivative(fn, xx, h, df1)
            xx(i) = x(i) - h(i)
            call numerical_derivative(fn, xx, h, df2)
            hess(:,i) = (df1 - df2) / (2.0_dp*h(i))
            xx(i) = x(i)
        end do
    end subroutine numerical_hessian

    subroutine numerical_derivative(fn, x, h, grad)
        procedure(subplex_objective) :: fn
        real(dp), intent(in) :: x(:), h(:)
        real(dp), intent(out) :: grad(:)
        real(dp), allocatable :: xx(:)
        real(dp) :: f1, f2
        integer :: i

        allocate(xx(size(x)))
        xx = x
        do i = 1, size(x)
            xx(i) = x(i) + h(i)
            f1 = fn(xx)
            xx(i) = x(i) - h(i)
            f2 = fn(xx)
            grad(i) = (f1 - f2) / (2.0_dp*h(i))
            xx(i) = x(i)
        end do
    end subroutine numerical_derivative

    pure logical function same_real(a, b) result(same)
        real(dp), intent(in) :: a, b
        same = (a <= b) .and. (a >= b)
    end function same_real

    pure logical function arrays_same_exact(a, b) result(same)
        real(dp), intent(in) :: a(:), b(:)
        integer :: i

        same = .false.
        if (size(a) /= size(b)) return
        same = .true.
        do i = 1, size(a)
            if (.not. same_real(a(i), b(i))) then
                same = .false.
                return
            end if
        end do
    end function arrays_same_exact

end module subplex
