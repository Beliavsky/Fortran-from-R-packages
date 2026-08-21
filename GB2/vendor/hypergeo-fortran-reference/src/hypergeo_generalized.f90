! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_generalized
    use hypergeo_kinds, only : dp
    use hypergeo_special, only : nan_complex, finite_complex
    use hypergeo_types, only : hypergeo_info
    use contfrac, only : gcf, contfrac_info
    implicit none
    private

    public :: genhypergeo, genhypergeo_series, genhypergeo_contfrac_single
    public :: genhypergeo_contfrac, genhypergeo_shanks, shanks_transform

contains

    function genhypergeo(u, l, z, tol, maxiter, series, info) result(value)
        complex(dp), intent(in) :: u(:), l(:), z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: series
        type(hypergeo_info), intent(out), optional :: info
        complex(dp) :: value
        logical :: use_series

        use_series = .true.
        if (present(series)) use_series = series
        if (use_series) then
            value = genhypergeo_series(u, l, z, tol=tol, maxiter=maxiter, info=info)
        else
            value = genhypergeo_contfrac_single(u, l, z, tol=tol, maxiter=maxiter, info=info)
        end if
    end function genhypergeo

    function genhypergeo_series(u, l, z, tol, maxiter, check_mod, polynomial, info) result(value)
        complex(dp), intent(in) :: u(:), l(:), z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: check_mod, polynomial
        type(hypergeo_info), intent(out), optional :: info
        complex(dp) :: value

        complex(dp), allocatable :: uu(:), ll(:)
        complex(dp) :: term, old, num, den
        real(dp) :: t, residual
        integer :: n, nmax, j, p, q
        logical :: check, poly, converged

        t = 0.0_dp
        if (present(tol)) t = max(0.0_dp, tol)
        nmax = 2000
        if (present(maxiter)) nmax = max(0, maxiter)
        check = .true.
        if (present(check_mod)) check = check_mod
        poly = .false.
        if (present(polynomial)) poly = polynomial

        p = size(u)
        q = size(l)
        if (check) then
            if (p > q + 1) then
                if (abs(z) > 0.0_dp) then
                    value = nan_complex()
                    call fill_info(info, .false., 0, huge(1.0_dp), 'series-domain', 1)
                    return
                end if
            else if (p == q + 1) then
                if (abs(z) > 1.0_dp) then
                    value = nan_complex()
                    call fill_info(info, .false., 0, huge(1.0_dp), 'series-domain', 1)
                    return
                end if
            end if
        end if

        allocate(uu(p), ll(q))
        uu = u
        ll = l
        term = (1.0_dp, 0.0_dp)
        value = term
        residual = huge(1.0_dp)
        converged = nmax == 0

        do n = 1, nmax
            num = (1.0_dp, 0.0_dp)
            den = (1.0_dp, 0.0_dp)
            do j = 1, p
                num = num * uu(j)
            end do
            do j = 1, q
                den = den * ll(j)
            end do
            if (q > 0 .and. abs(den) <= tiny(1.0_dp)) then
                value = nan_complex()
                call fill_info(info, .false., n, huge(1.0_dp), 'series-pole', 1)
                return
            end if
            old = value
            term = term * num / den * z / real(n, dp)
            value = old + term
            residual = abs(value - old)
            if (.not. finite_complex(value)) exit
            if (residual <= t) then
                converged = .true.
                exit
            end if
            if (.not. (abs(term) > tiny(1.0_dp))) then
                converged = .true.
                exit
            end if
            uu = uu + 1.0_dp
            ll = ll + 1.0_dp
        end do

        if (.not. converged .and. .not. poly) value = nan_complex()
        call fill_info(info, converged .or. poly, min(n, nmax), residual, 'series', 1)
    end function genhypergeo_series

    function genhypergeo_contfrac_single(u, l, z, tol, maxiter, info) result(value)
        complex(dp), intent(in) :: u(:), l(:), z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        type(hypergeo_info), intent(out), optional :: info
        complex(dp) :: value

        complex(dp), allocatable :: alpha(:), aa(:), bb(:)
        complex(dp) :: prod_u, prod_l, den, cf
        integer :: k, j, nmax
        real(dp) :: t
        type(contfrac_info) :: cfi

        nmax = 2000
        if (present(maxiter)) nmax = max(1, maxiter)
        t = epsilon(1.0_dp)
        if (present(tol)) then
            if (tol > 0.0_dp) t = tol
        end if
        allocate(alpha(nmax), aa(nmax), bb(nmax))
        do k = 1, nmax
            prod_u = (1.0_dp, 0.0_dp)
            prod_l = cmplx(real(k + 1, dp), 0.0_dp, kind=dp)
            do j = 1, size(u)
                prod_u = prod_u * (u(j) + real(k, dp))
            end do
            do j = 1, size(l)
                prod_l = prod_l * (l(j) + real(k, dp))
            end do
            alpha(k) = z * prod_u / prod_l
        end do
        aa = -alpha
        bb = 1.0_dp + alpha
        cf = gcf(aa, bb, b0=cmplx(0.0_dp, 0.0_dp, kind=dp), tol=t, info=cfi)

        prod_u = (1.0_dp, 0.0_dp)
        prod_l = (1.0_dp, 0.0_dp)
        do j = 1, size(u)
            prod_u = prod_u * u(j)
        end do
        do j = 1, size(l)
            prod_l = prod_l * l(j)
        end do
        den = prod_l * (1.0_dp + cf)
        if (abs(den) <= tiny(1.0_dp) .or. .not. finite_complex(cf)) then
            value = nan_complex()
        else
            value = 1.0_dp + z * prod_u / den
        end if
        call fill_info(info, finite_complex(value), cfi%iterations, cfi%residual, 'continued-fraction', 2)
    end function genhypergeo_contfrac_single

    function genhypergeo_contfrac(u, l, z, tol, maxiter) result(values)
        complex(dp), intent(in) :: u(:), l(:), z(:)
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: values(size(z))
        integer :: i
        do i = 1, size(z)
            values(i) = genhypergeo_contfrac_single(u, l, z(i), tol=tol, maxiter=maxiter)
        end do
    end function genhypergeo_contfrac

    pure function shanks_transform(last, this, next) result(value)
        complex(dp), intent(in) :: last, this, next
        complex(dp) :: value, den
        if (.not. (abs(next - this) > tiny(1.0_dp))) then
            value = next
            return
        end if
        den = next - 2.0_dp * this + last
        if (abs(den) <= tiny(1.0_dp)) then
            value = next
        else
            value = (next * last - this * this) / den
        end if
    end function shanks_transform

    function genhypergeo_shanks(u, l, z, maxiter) result(value)
        complex(dp), intent(in) :: u(:), l(:), z
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        complex(dp), allocatable :: uu(:), ll(:)
        complex(dp) :: term, old, num, den, last, this, next, accelerated
        integer :: n, j, nmax

        nmax = 20
        if (present(maxiter)) nmax = max(0, maxiter)
        allocate(uu(size(u)), ll(size(l)))
        uu = u
        ll = l
        term = 1.0_dp
        value = 1.0_dp
        last = 0.0_dp
        this = 1.0_dp
        next = 2.0_dp
        accelerated = shanks_transform(last, this, next)
        if (nmax == 0) return
        do n = 1, nmax
            num = 1.0_dp
            den = 1.0_dp
            do j = 1, size(uu)
                num = num * uu(j)
            end do
            do j = 1, size(ll)
                den = den * ll(j)
            end do
            old = value
            term = term * num / den * z / real(n, dp)
            value = old + term
            last = this
            this = next
            next = value
            accelerated = shanks_transform(last, this, next)
            uu = uu + 1.0_dp
            ll = ll + 1.0_dp
        end do
        ! The upstream genhypergeo_shanks() computes the accelerated value
        ! diagnostically but returns the ordinary partial sum. Preserve that API.
        if (.not. finite_complex(accelerated)) continue
    end function genhypergeo_shanks

    subroutine fill_info(info, converged, iterations, residual, method, code)
        type(hypergeo_info), intent(out), optional :: info
        logical, intent(in) :: converged
        integer, intent(in) :: iterations, code
        real(dp), intent(in) :: residual
        character(len=*), intent(in) :: method
        if (.not. present(info)) return
        info%converged = converged
        info%iterations = iterations
        info%residual = residual
        info%method = method
        info%method_code = code
    end subroutine fill_info

end module hypergeo_generalized
