! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_optimize
    use hyper2_kinds, only : dp
    use hyper2_types, only : hyper2_model, hyper3_model
    use hyper2_likelihood, only : loglik_h2, loglik_h3, gradient_h2, gradient_h3, &
        hessian_independent_h2, hessian_independent_h3, fillup, equalp
    implicit none
    private

    type, public :: fit_result
        real(dp), allocatable :: p(:)
        real(dp) :: log_likelihood = -huge(1.0_dp)
        integer :: iterations = 0
        logical :: converged = .false.
        integer :: status = 0
    end type fit_result

    interface maxp
        module procedure maxp_h2
        module procedure maxp_h3
    end interface
    interface maxp_simplex
        module procedure maxp_simplex_h2
        module procedure maxp_simplex_h3
    end interface

    public :: maxp, maxp_h2, maxp_h3, maxp_simplex, feasible_simplex

contains

    logical function feasible_simplex(x, small, fcm, fcv)
        real(dp), intent(in) :: x(:), small
        real(dp), intent(in), optional :: fcm(:,:), fcv(:)
        feasible_simplex = all(x >= small) .and. (1.0_dp - sum(x) >= small)
        if (.not. feasible_simplex) return
        if (present(fcm) .neqv. present(fcv)) then
            feasible_simplex = .false.
            return
        end if
        if (present(fcm)) then
            if (size(fcm,2) /= size(x) .or. size(fcm,1) /= size(fcv)) then
                feasible_simplex = .false.
                return
            end if
            feasible_simplex = all(matmul(fcm, x) >= fcv)
        end if
    end function feasible_simplex

    subroutine solve_linear(a_in, b_in, x, ok)
        real(dp), intent(in) :: a_in(:,:), b_in(:)
        real(dp), allocatable, intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: a(:,:), b(:), row(:)
        real(dp) :: piv, fac, tb
        integer :: n, i, j, k, imax
        n = size(b_in)
        allocate(x(n))
        x = 0.0_dp
        ok = size(a_in,1) == n .and. size(a_in,2) == n
        if (.not. ok) return
        a = a_in
        b = b_in
        allocate(row(n))
        do k = 1, n
            imax = k
            do i = k + 1, n
                if (abs(a(i,k)) > abs(a(imax,k))) imax = i
            end do
            if (abs(a(imax,k)) <= 1.0e-14_dp) then
                ok = .false.
                return
            end if
            if (imax /= k) then
                row = a(k,:)
                a(k,:) = a(imax,:)
                a(imax,:) = row
                tb = b(k)
                b(k) = b(imax)
                b(imax) = tb
            end if
            piv = a(k,k)
            do i = k + 1, n
                fac = a(i,k) / piv
                a(i,k:n) = a(i,k:n) - fac * a(k,k:n)
                b(i) = b(i) - fac * b(k)
            end do
        end do
        do i = n, 1, -1
            x(i) = b(i)
            do j = i + 1, n
                x(i) = x(i) - a(i,j) * x(j)
            end do
            x(i) = x(i) / a(i,i)
        end do
    end subroutine solve_linear

    subroutine initial_x(n, startp, small, x)
        integer, intent(in) :: n
        real(dp), intent(in), optional :: startp(:)
        real(dp), intent(in) :: small
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), allocatable :: p(:)
        integer :: m
        m = n - 1
        allocate(x(max(0,m)), p(max(0,n)))
        if (m == 0) return
        if (present(startp)) then
            if (size(startp) == n) then
                x = startp(1:m)
            else if (size(startp) == m) then
                x = startp
            else
                p = 1.0_dp / real(n, dp)
                x = p(1:m)
            end if
        else
            p = 1.0_dp / real(n, dp)
            x = p(1:m)
        end if
        if (.not. feasible_simplex(x, small)) then
            p = 1.0_dp / real(n, dp)
            x = p(1:m)
        end if
    end subroutine initial_x

    function maxp_h2(h, startp, fcm, fcv, small, tol, max_iter) result(res)
        type(hyper2_model), intent(in) :: h
        real(dp), intent(in), optional :: startp(:), fcm(:,:), fcv(:), small, tol
        integer, intent(in), optional :: max_iter
        type(fit_result) :: res
        real(dp), allocatable :: x(:), p(:), g(:), hs(:,:), a(:,:), d(:), xn(:), pn(:)
        real(dp) :: sm, tt, ll, lln, alpha, slope, ridge
        integer :: it, mi, n, j
        logical :: ok, accepted
        sm = 1.0e-8_dp
        if (present(small)) sm = small
        tt = 1.0e-9_dp
        if (present(tol)) tt = tol
        mi = 500
        if (present(max_iter)) mi = max_iter
        n = h%size()
        if (n < 1) then
            allocate(res%p(0))
            res%status = 2
            return
        end if
        if (n == 1) then
            allocate(res%p(1))
            res%p = 1.0_dp
            res%log_likelihood = loglik_h2(res%p,h)
            res%converged = .true.
            return
        end if
        call initial_x(n, startp, sm, x)
        if (.not. feasible_simplex(x, sm, fcm, fcv)) then
            res%status = 3
            res%p = fillup(x)
            return
        end if
        p = fillup(x)
        ll = loglik_h2(p, h)
        do it = 1, mi
            g = gradient_h2(h, p)
            if (maxval(abs(g)) <= tt) then
                res%converged = .true.
                exit
            end if
            hs = hessian_independent_h2(h, p)
            allocate(a(size(hs,1),size(hs,2)))
            a = -hs
            ridge = 1.0e-10_dp
            do j = 1, size(a,1)
                a(j,j) = a(j,j) + ridge
            end do
            call solve_linear(a, g, d, ok)
            deallocate(a)
            if (.not. ok .or. dot_product(g,d) <= 0.0_dp) then
                d = g / max(1.0_dp, maxval(abs(g)))
            end if
            slope = dot_product(g, d)
            alpha = 1.0_dp
            accepted = .false.
            do while (alpha >= 1.0e-12_dp)
                xn = x + alpha*d
                if (feasible_simplex(xn, sm, fcm, fcv)) then
                    pn = fillup(xn)
                    lln = loglik_h2(pn, h)
                    if (lln >= ll + 1.0e-4_dp*alpha*slope) then
                        accepted = .true.
                        exit
                    end if
                end if
                alpha = alpha * 0.5_dp
            end do
            if (.not. accepted) then
                d = g / max(1.0_dp, maxval(abs(g)))
                alpha = 1.0e-2_dp
                do while (alpha >= 1.0e-14_dp)
                    xn = x + alpha*d
                    if (feasible_simplex(xn, sm, fcm, fcv)) then
                        pn = fillup(xn)
                        lln = loglik_h2(pn, h)
                        if (lln > ll) then
                            accepted = .true.
                            exit
                        end if
                    end if
                    alpha = alpha * 0.5_dp
                end do
            end if
            if (.not. accepted) exit
            if (maxval(abs(xn-x)) <= tt .and. abs(lln-ll) <= tt) then
                x = xn
                p = pn
                ll = lln
                res%converged = .true.
                exit
            end if
            x = xn
            p = pn
            ll = lln
        end do
        res%p = p
        res%log_likelihood = ll
        res%iterations = min(it, mi)
        if (.not. res%converged) then
            g = gradient_h2(h, p)
            if (maxval(abs(g)) < sqrt(tt)) res%converged = .true.
        end if
        if (.not. res%converged) res%status = 1
    end function maxp_h2

    function maxp_h3(h, startp, fcm, fcv, small, tol, max_iter) result(res)
        type(hyper3_model), intent(in) :: h
        real(dp), intent(in), optional :: startp(:), fcm(:,:), fcv(:), small, tol
        integer, intent(in), optional :: max_iter
        type(fit_result) :: res
        real(dp), allocatable :: x(:), p(:), g(:), hs(:,:), a(:,:), d(:), xn(:), pn(:)
        real(dp) :: sm, tt, ll, lln, alpha, slope
        integer :: it, mi, n, j
        logical :: ok, accepted
        sm = 1.0e-8_dp
        if (present(small)) sm = small
        tt = 1.0e-9_dp
        if (present(tol)) tt = tol
        mi = 500
        if (present(max_iter)) mi = max_iter
        n = h%size()
        if (n < 1) then
            allocate(res%p(0))
            res%status = 2
            return
        end if
        if (n == 1) then
            allocate(res%p(1))
            res%p = 1.0_dp
            res%log_likelihood = loglik_h3(res%p,h)
            res%converged = .true.
            return
        end if
        call initial_x(n, startp, sm, x)
        if (.not. feasible_simplex(x, sm, fcm, fcv)) then
            res%status = 3
            res%p = fillup(x)
            return
        end if
        p = fillup(x)
        ll = loglik_h3(p, h)
        do it = 1, mi
            g = gradient_h3(h, p)
            if (maxval(abs(g)) <= tt) then
                res%converged = .true.
                exit
            end if
            hs = hessian_independent_h3(h, p)
            allocate(a(size(hs,1),size(hs,2)))
            a = -hs
            do j = 1, size(a,1)
                a(j,j) = a(j,j) + 1.0e-10_dp
            end do
            call solve_linear(a, g, d, ok)
            deallocate(a)
            if (.not. ok .or. dot_product(g,d) <= 0.0_dp) d = g / max(1.0_dp,maxval(abs(g)))
            slope = dot_product(g,d)
            alpha = 1.0_dp
            accepted = .false.
            do while (alpha >= 1.0e-12_dp)
                xn = x + alpha*d
                if (feasible_simplex(xn, sm, fcm, fcv)) then
                    pn = fillup(xn)
                    lln = loglik_h3(pn,h)
                    if (lln >= ll + 1.0e-4_dp*alpha*slope) then
                        accepted = .true.
                        exit
                    end if
                end if
                alpha = alpha*0.5_dp
            end do
            if (.not. accepted) exit
            if (maxval(abs(xn-x)) <= tt .and. abs(lln-ll) <= tt) then
                x=xn
                p=pn
                ll=lln
                res%converged=.true.
                exit
            end if
            x=xn
            p=pn
            ll=lln
        end do
        res%p=p
        res%log_likelihood=ll
        res%iterations=min(it,mi)
        if (.not. res%converged) then
            g = gradient_h3(h, p)
            if (maxval(abs(g)) < sqrt(tt)) res%converged = .true.
        end if
        if (.not. res%converged) res%status=1
    end function maxp_h3

    function random_simplex(n) result(p)
        integer, intent(in) :: n
        real(dp), allocatable :: p(:)
        real(dp) :: u
        integer :: i
        allocate(p(n))
        do i=1,n
            call random_number(u)
            u=max(u,tiny(1.0_dp))
            p(i)=-log(u)
        end do
        p=p/sum(p)
    end function random_simplex

    function maxp_simplex_h2(h, nstart, small, tol, max_iter) result(best)
        type(hyper2_model), intent(in) :: h
        integer, intent(in), optional :: nstart, max_iter
        real(dp), intent(in), optional :: small, tol
        type(fit_result) :: best, cur
        real(dp), allocatable :: p(:)
        integer :: i, ns
        ns=100
        if(present(nstart)) ns=nstart
        best%log_likelihood=-huge(1.0_dp)
        do i=1,ns
            p=random_simplex(h%size())
            cur=maxp_h2(h,startp=p,small=small,tol=tol,max_iter=max_iter)
            if(cur%log_likelihood>best%log_likelihood) best=cur
        end do
    end function maxp_simplex_h2

    function maxp_simplex_h3(h, nstart, small, tol, max_iter) result(best)
        type(hyper3_model), intent(in) :: h
        integer, intent(in), optional :: nstart, max_iter
        real(dp), intent(in), optional :: small, tol
        type(fit_result) :: best, cur
        real(dp), allocatable :: p(:)
        integer :: i, ns
        ns=100
        if(present(nstart)) ns=nstart
        best%log_likelihood=-huge(1.0_dp)
        do i=1,ns
            p=random_simplex(h%size())
            cur=maxp_h3(h,startp=p,small=small,tol=tol,max_iter=max_iter)
            if(cur%log_likelihood>best%log_likelihood) best=cur
        end do
    end function maxp_simplex_h3

end module hyper2_optimize
