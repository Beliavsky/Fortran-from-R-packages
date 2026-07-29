! SPDX-License-Identifier: GPL-2.0-or-later
module evir_distributions
    use evir_kinds, only : dp
    use evir_types, only : evir_rng
    use evir_math, only : random_uniform, safe_nan
    implicit none
    private

    public :: dgev, pgev, qgev, rgev
    public :: dgpd, pgpd, qgpd, rgpd

contains

    elemental real(dp) function dgev(x, xi, mu, sigma) result(d)
        real(dp), intent(in) :: x
        real(dp), intent(in), optional :: xi, mu, sigma
        real(dp) :: sh, loc, sc, z, t
        sh = 1.0_dp
        loc = 0.0_dp
        sc = 1.0_dp
        if (present(xi)) sh = xi
        if (present(mu)) loc = mu
        if (present(sigma)) sc = sigma
        if (sc <= 0.0_dp) then
            d = safe_nan()
            return
        end if
        z = (x-loc)/sc
        if (abs(sh) <= 1.0e-10_dp) then
            d = exp(-z-exp(-z))/sc
        else
            t = 1.0_dp + sh*z
            if (t <= 0.0_dp) then
                d = 0.0_dp
            else
                d = t**(-1.0_dp/sh-1.0_dp)*exp(-t**(-1.0_dp/sh))/sc
            end if
        end if
    end function dgev

    elemental real(dp) function pgev(q, xi, mu, sigma) result(p)
        real(dp), intent(in) :: q
        real(dp), intent(in), optional :: xi, mu, sigma
        real(dp) :: sh, loc, sc, z, t
        sh = 1.0_dp
        loc = 0.0_dp
        sc = 1.0_dp
        if (present(xi)) sh = xi
        if (present(mu)) loc = mu
        if (present(sigma)) sc = sigma
        if (sc <= 0.0_dp) then
            p = safe_nan()
            return
        end if
        z = (q-loc)/sc
        if (abs(sh) <= 1.0e-10_dp) then
            p = exp(-exp(-z))
        else
            t = 1.0_dp + sh*z
            if (t <= 0.0_dp) then
                if (sh > 0.0_dp) then
                    p = 0.0_dp
                else
                    p = 1.0_dp
                end if
            else
                p = exp(-t**(-1.0_dp/sh))
            end if
        end if
    end function pgev

    elemental real(dp) function qgev(p, xi, mu, sigma) result(q)
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: xi, mu, sigma
        real(dp) :: sh, loc, sc
        sh = 1.0_dp
        loc = 0.0_dp
        sc = 1.0_dp
        if (present(xi)) sh = xi
        if (present(mu)) loc = mu
        if (present(sigma)) sc = sigma
        if (sc <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            q = safe_nan()
        else if (p <= 0.0_dp) then
            if (sh > 0.0_dp) then
                q = loc-sc/sh
            else
                q = -huge(1.0_dp)
            end if
        else if (p >= 1.0_dp) then
            if (sh < 0.0_dp) then
                q = loc-sc/sh
            else
                q = huge(1.0_dp)
            end if
        else if (abs(sh) <= 1.0e-10_dp) then
            q = loc-sc*log(-log(p))
        else
            q = loc+(sc/sh)*((-log(p))**(-sh)-1.0_dp)
        end if
    end function qgev

    subroutine rgev(n, values, rng, xi, mu, sigma)
        integer, intent(in) :: n
        real(dp), intent(out) :: values(n)
        type(evir_rng), intent(inout) :: rng
        real(dp), intent(in), optional :: xi, mu, sigma
        integer :: i
        do i = 1, n
            values(i) = qgev(random_uniform(rng), xi, mu, sigma)
        end do
    end subroutine rgev

    elemental real(dp) function dgpd(x, xi, mu, beta) result(d)
        real(dp), intent(in) :: x, xi
        real(dp), intent(in), optional :: mu, beta
        real(dp) :: loc, sc, z, t
        loc = 0.0_dp
        sc = 1.0_dp
        if (present(mu)) loc = mu
        if (present(beta)) sc = beta
        if (sc <= 0.0_dp) then
            d = safe_nan()
            return
        end if
        z = (x-loc)/sc
        if (z < 0.0_dp) then
            d = 0.0_dp
        else if (abs(xi) <= 1.0e-10_dp) then
            d = exp(-z)/sc
        else
            t = 1.0_dp+xi*z
            if (t <= 0.0_dp) then
                d = 0.0_dp
            else
                d = t**(-1.0_dp/xi-1.0_dp)/sc
            end if
        end if
    end function dgpd

    elemental real(dp) function pgpd(q, xi, mu, beta) result(p)
        real(dp), intent(in) :: q, xi
        real(dp), intent(in), optional :: mu, beta
        real(dp) :: loc, sc, z, t
        loc = 0.0_dp
        sc = 1.0_dp
        if (present(mu)) loc = mu
        if (present(beta)) sc = beta
        if (sc <= 0.0_dp) then
            p = safe_nan()
            return
        end if
        z = (q-loc)/sc
        if (z < 0.0_dp) then
            p = 0.0_dp
        else if (abs(xi) <= 1.0e-10_dp) then
            p = 1.0_dp-exp(-z)
        else
            t = 1.0_dp+xi*z
            if (t <= 0.0_dp) then
                p = 1.0_dp
            else
                p = 1.0_dp-t**(-1.0_dp/xi)
            end if
        end if
    end function pgpd

    elemental real(dp) function qgpd(p, xi, mu, beta) result(q)
        real(dp), intent(in) :: p, xi
        real(dp), intent(in), optional :: mu, beta
        real(dp) :: loc, sc
        loc = 0.0_dp
        sc = 1.0_dp
        if (present(mu)) loc = mu
        if (present(beta)) sc = beta
        if (sc <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            q = safe_nan()
        else if (p <= 0.0_dp) then
            q = loc
        else if (p >= 1.0_dp) then
            if (xi < 0.0_dp) then
                q = loc-sc/xi
            else
                q = huge(1.0_dp)
            end if
        else if (abs(xi) <= 1.0e-10_dp) then
            q = loc-sc*log(1.0_dp-p)
        else
            q = loc+(sc/xi)*((1.0_dp-p)**(-xi)-1.0_dp)
        end if
    end function qgpd

    subroutine rgpd(n, values, rng, xi, mu, beta)
        integer, intent(in) :: n
        real(dp), intent(out) :: values(n)
        type(evir_rng), intent(inout) :: rng
        real(dp), intent(in) :: xi
        real(dp), intent(in), optional :: mu, beta
        integer :: i
        do i = 1, n
            values(i) = qgpd(random_uniform(rng), xi, mu, beta)
        end do
    end subroutine rgpd

end module evir_distributions
