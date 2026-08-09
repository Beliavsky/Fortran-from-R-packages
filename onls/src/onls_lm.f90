! SPDX-License-Identifier: GPL-2.0-or-later
module onls_lm
    use onls_kinds, only : dp
    use onls_linalg, only : solve_spd
    implicit none
    private
    public :: lm_control, lm_result, residual_fn, lm_solve

    type :: lm_control
        integer :: maxiter = 1024
        integer :: maxfev = 100000
        real(dp) :: ftol = sqrt(epsilon(1.0_dp))
        real(dp) :: xtol = sqrt(epsilon(1.0_dp))
        real(dp) :: gtol = 1.0e-8_dp
        real(dp) :: fd_step = sqrt(epsilon(1.0_dp))
        real(dp) :: lambda0 = 1.0e-3_dp
    end type lm_control

    type :: lm_result
        real(dp), allocatable :: par(:)
        real(dp), allocatable :: fvec(:)
        real(dp) :: rss = huge(1.0_dp)
        integer :: niter = 0
        integer :: nfev = 0
        integer :: info = 0
        character(len=128) :: message = ''
    end type lm_result

    abstract interface
        subroutine residual_fn(par, r, ctx, ierr)
            import dp
            real(dp), intent(in) :: par(:)
            real(dp), intent(out) :: r(:)
            class(*), intent(inout) :: ctx
            integer, intent(out) :: ierr
        end subroutine residual_fn
    end interface
contains
    subroutine lm_solve(fn, ctx, start, m, result, control, lower, upper)
        procedure(residual_fn) :: fn
        class(*), intent(inout) :: ctx
        real(dp), intent(in) :: start(:)
        integer, intent(in) :: m
        type(lm_result), intent(out) :: result
        type(lm_control), intent(in), optional :: control
        real(dp), intent(in), optional :: lower(:), upper(:)
        type(lm_control) :: ctl
        real(dp), allocatable :: p(:), pnew(:), r(:), rnew(:), jmat(:,:), h(:,:), g(:), step(:), d(:)
        real(dp), allocatable :: lo(:), hi(:), pp(:), pm(:), rp(:), rm(:)
        real(dp) :: rss, rssnew, lambda, hstep, pred, ared, rho, scale
        integer :: n, i, k, ierr
        logical :: ok

        ctl = lm_control()
        if (present(control)) ctl = control
        n = size(start)
        allocate(p(n), pnew(n), r(m), rnew(m), jmat(m,n), h(n,n), g(n), step(n), d(n))
        allocate(lo(n), hi(n), pp(n), pm(n), rp(m), rm(m))
        lo = -huge(1.0_dp)
        hi = huge(1.0_dp)
        if (present(lower)) then
            if (size(lower) == 1) then
                lo = lower(1)
            else if (size(lower) == n) then
                lo = lower
            end if
        end if
        if (present(upper)) then
            if (size(upper) == 1) then
                hi = upper(1)
            else if (size(upper) == n) then
                hi = upper
            end if
        end if
        p = min(hi, max(lo, start))
        call fn(p, r, ctx, ierr)
        result%nfev = 1
        if (ierr /= 0 .or. any(.not.(abs(r) <= huge(1.0_dp)))) then
            result%info = -1
            result%message = 'Initial residual evaluation failed'
            result%par = p
            result%fvec = r
            return
        end if
        rss = dot_product(r,r)
        lambda = ctl%lambda0
        do k = 1, ctl%maxiter
            result%niter = k
            if (result%nfev >= ctl%maxfev) then
                result%info = 5
                result%message = 'Maximum function evaluations reached'
                exit
            end if
            do i = 1, n
                hstep = ctl%fd_step * max(1.0_dp, abs(p(i)))
                if (p(i) + hstep <= hi(i) .and. p(i) - hstep >= lo(i)) then
                    pp = p
                    pm = p
                    pp(i) = pp(i) + hstep
                    pm(i) = pm(i) - hstep
                    call fn(pp, rp, ctx, ierr)
                    if (ierr /= 0) exit
                    call fn(pm, rm, ctx, ierr)
                    if (ierr /= 0) exit
                    result%nfev = result%nfev + 2
                    jmat(:,i) = (rp-rm)/(2.0_dp*hstep)
                else if (p(i) + hstep <= hi(i)) then
                    pp = p
                    pp(i) = pp(i) + hstep
                    call fn(pp, rp, ctx, ierr)
                    if (ierr /= 0) exit
                    result%nfev = result%nfev + 1
                    jmat(:,i) = (rp-r)/hstep
                else if (p(i) - hstep >= lo(i)) then
                    pm = p
                    pm(i) = pm(i) - hstep
                    call fn(pm, rm, ctx, ierr)
                    if (ierr /= 0) exit
                    result%nfev = result%nfev + 1
                    jmat(:,i) = (r-rm)/hstep
                else
                    jmat(:,i) = 0.0_dp
                end if
            end do
            if (ierr /= 0) then
                result%info = -2
                result%message = 'Residual evaluation failed during differentiation'
                exit
            end if
            g = matmul(transpose(jmat), r)
            if (maxval(abs(g)) <= ctl%gtol * max(1.0_dp, sqrt(rss))) then
                result%info = 4
                result%message = 'Gradient orthogonality tolerance reached'
                exit
            end if
            h = matmul(transpose(jmat), jmat)
            d = max(1.0_dp, [(h(i,i), i=1,n)])
            h = h + lambda * diag_matrix(d)
            call solve_spd(h, -g, step, ok)
            if (.not. ok) then
                lambda = 10.0_dp * max(lambda, 1.0e-12_dp)
                cycle
            end if
            pnew = min(hi, max(lo, p + step))
            if (maxval(abs(pnew-p)) <= ctl%xtol * (ctl%xtol + maxval(abs(p)))) then
                p = pnew
                result%info = 2
                result%message = 'Parameter tolerance reached'
                exit
            end if
            call fn(pnew, rnew, ctx, ierr)
            result%nfev = result%nfev + 1
            if (ierr /= 0) then
                lambda = 10.0_dp * max(lambda, 1.0e-12_dp)
                cycle
            end if
            rssnew = dot_product(rnew,rnew)
            ared = 0.5_dp * (rss-rssnew)
            pred = -dot_product(g, pnew-p) - 0.5_dp * dot_product(pnew-p, &
                   matmul(matmul(transpose(jmat),jmat), pnew-p))
            if (pred <= 0.0_dp) then
                rho = -1.0_dp
            else
                rho = ared/pred
            end if
            if (rssnew < rss) then
                scale = max(1.0_dp, rss)
                p = pnew
                r = rnew
                if (abs(rss-rssnew) <= ctl%ftol * scale) then
                    rss = rssnew
                    result%info = 1
                    result%message = 'Function tolerance reached'
                    exit
                end if
                rss = rssnew
                if (rho > 0.75_dp) then
                    lambda = max(1.0e-15_dp, 0.3333333333333333_dp * lambda)
                else if (rho < 0.25_dp) then
                    lambda = min(1.0e15_dp, 2.0_dp * lambda)
                end if
            else
                lambda = min(1.0e15_dp, 10.0_dp * max(lambda, 1.0e-12_dp))
            end if
        end do
        if (result%info == 0) then
            result%info = 5
            result%message = 'Maximum iterations reached'
        end if
        result%par = p
        result%fvec = r
        result%rss = dot_product(r,r)
    contains
        pure function diag_matrix(v) result(a)
            real(dp), intent(in) :: v(:)
            real(dp) :: a(size(v),size(v))
            integer :: ii
            a = 0.0_dp
            do ii = 1, size(v)
                a(ii,ii) = v(ii)
            end do
        end function diag_matrix
    end subroutine lm_solve
end module onls_lm
