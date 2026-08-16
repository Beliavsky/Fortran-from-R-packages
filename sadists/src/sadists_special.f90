! SPDX-License-Identifier: LGPL-3.0-or-later
! Numerical support specific to the sadists translation.  Common normal
! probability utilities and moment/cumulant conversion come from the
! PDQutils-fortran dependency; the remaining routines here are sadists-specific.
module sadists_special
    use sadists_kinds, only : dp, pi, log2
    use pdqutils, only : moment2cumulant, cumulant2moment
    use pdqutils_special, only : nan_dp, pos_inf_dp, neg_inf_dp, factorial_dp, binomial_dp, &
        normal_pdf, normal_cdf, normal_quantile
    implicit none
    private

    public :: nan_dp, pos_inf_dp, neg_inf_dp, factorial_dp, binomial_dp
    public :: normal_pdf, normal_cdf, normal_quantile
    public :: log_gamma_ratio, digamma_dp, polygamma_dp
    public :: moments_to_cumulants, cumulants_to_moments
    public :: normal_moments, chisq_moments, chisq_log_moment
    public :: logchisq_cumulants, lognoncentral_chisq_cumulants
    public :: random_normal, random_gamma, random_poisson, random_chisq
    public :: fill_normal, fill_chisq

contains

    pure real(dp) function log_gamma_ratio(x, s) result(v)
        ! log(Gamma(x+s)/Gamma(x)); asymptotic difference avoids catastrophic
        ! cancellation for very large x, relevant to sadists >= 0.2.5.
        real(dp), intent(in) :: x, s
        real(dp) :: ix
        if (x <= 0.0_dp .or. x+s <= 0.0_dp) then
            v = pos_inf_dp()
            return
        end if
        if (x < 1.0e5_dp) then
            v = log_gamma(x+s) - log_gamma(x)
        else
            ix = 1.0_dp/x
            v = s*log(x) &
                + (0.5_dp*s*s-0.5_dp*s)*ix &
                + (-s**3/6.0_dp+s*s/4.0_dp-s/12.0_dp)*ix**2 &
                + (s**4/12.0_dp-s**3/6.0_dp+s*s/12.0_dp)*ix**3 &
                + (-s**5/20.0_dp+s**4/8.0_dp-s**3/12.0_dp+s/120.0_dp)*ix**4 &
                + (s**6/30.0_dp-s**5/10.0_dp+s**4/12.0_dp-s*s/60.0_dp)*ix**5
        end if
    end function log_gamma_ratio

    pure real(dp) function digamma_dp(xin) result(v)
        real(dp), intent(in) :: xin
        real(dp) :: x, inv, inv2
        if (xin <= 0.0_dp) then
            v = nan_dp(); return
        end if
        x = xin
        v = 0.0_dp
        do while (x < 12.0_dp)
            v = v - 1.0_dp/x
            x = x + 1.0_dp
        end do
        inv = 1.0_dp/x
        inv2 = inv*inv
        v = v + log(x) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp &
            - inv2*(1.0_dp/120.0_dp - inv2*(1.0_dp/252.0_dp &
            - inv2*(1.0_dp/240.0_dp - inv2*(5.0_dp/660.0_dp)))))
    end function digamma_dp

    pure real(dp) function polygamma_dp(m, xin) result(v)
        integer, intent(in) :: m
        real(dp), intent(in) :: xin
        real(dp), parameter :: bern(6) = [1.0_dp/6.0_dp, -1.0_dp/30.0_dp, &
            1.0_dp/42.0_dp, -1.0_dp/30.0_dp, 5.0_dp/66.0_dp, -691.0_dp/2730.0_dp]
        real(dp) :: x, term, rf
        integer :: k, j
        if (m == 0) then
            v = digamma_dp(xin); return
        end if
        if (m < 0 .or. xin <= 0.0_dp) then
            v = nan_dp(); return
        end if
        x = xin
        v = 0.0_dp
        do while (x < 12.0_dp)
            v = v + (-1.0_dp)**(m+1) * factorial_dp(m) / x**(m+1)
            x = x + 1.0_dp
        end do
        ! m-th derivative of log(x) - 1/(2x) - sum B_2k/(2k x^(2k)).
        v = v + (-1.0_dp)**(m-1) * factorial_dp(m-1) / x**m
        v = v + (-1.0_dp)**(m+1) * factorial_dp(m) / (2.0_dp*x**(m+1))
        do k = 1, size(bern)
            rf = 1.0_dp
            do j = 0, m-1
                rf = rf * real(2*k+j,dp)
            end do
            term = -bern(k)/real(2*k,dp) * (-1.0_dp)**m * rf / x**(2*k+m)
            v = v + term
        end do
    end function polygamma_dp

    pure subroutine moments_to_cumulants(moms, kappa)
        real(dp), intent(in) :: moms(:)
        real(dp), intent(out) :: kappa(size(moms))
        kappa = moment2cumulant(moms)
    end subroutine moments_to_cumulants

    pure subroutine cumulants_to_moments(kappa, moms)
        real(dp), intent(in) :: kappa(:)
        real(dp), intent(out) :: moms(size(kappa))
        moms = cumulant2moment(kappa)
    end subroutine cumulants_to_moments

    pure subroutine normal_moments(mu, sigma, moms)
        real(dp), intent(in) :: mu, sigma
        real(dp), intent(out) :: moms(:)
        real(dp) :: m0, m1, mn
        integer :: n
        if (size(moms) == 0) return
        m0 = 1.0_dp
        m1 = mu
        moms(1) = m1
        do n = 2, size(moms)
            mn = mu*m1 + real(n-1,dp)*sigma*sigma*m0
            moms(n) = mn
            m0 = m1
            m1 = mn
        end do
    end subroutine normal_moments

    pure real(dp) function chisq_log_moment(df, ncp, ord) result(v)
        real(dp), intent(in) :: df, ncp, ord
        real(dp) :: a, lambda, width, lt, maxlt, sumexp
        integer :: j, jlo, jhi
        if (df < 0.0_dp .or. ncp < 0.0_dp) then
            v = nan_dp(); return
        end if
        a = 0.5_dp*df
        if (ncp == 0.0_dp) then
            if (df == 0.0_dp) then
                if (ord > 0.0_dp) then
                    v = neg_inf_dp()
                else if (ord < 0.0_dp) then
                    v = pos_inf_dp()
                else
                    v = 0.0_dp
                end if
                return
            end if
            if (a+ord <= 0.0_dp) then
                v = pos_inf_dp(); return
            end if
            v = ord*log2 + log_gamma_ratio(a,ord)
            return
        end if
        if (a+ord <= 0.0_dp) then
            v = pos_inf_dp(); return
        end if
        lambda = 0.5_dp*ncp
        width = 12.0_dp*sqrt(lambda+1.0_dp) + 30.0_dp
        jlo = max(0, int(floor(lambda-width)))
        jhi = max(jlo+1, int(ceiling(lambda+width)))
        maxlt = neg_inf_dp()
        do j = jlo, jhi
            lt = -lambda + real(j,dp)*log(lambda) - log_gamma(real(j+1,dp)) &
                + log_gamma_ratio(a+real(j,dp),ord)
            maxlt = max(maxlt, lt)
        end do
        sumexp = 0.0_dp
        do j = jlo, jhi
            lt = -lambda + real(j,dp)*log(lambda) - log_gamma(real(j+1,dp)) &
                + log_gamma_ratio(a+real(j,dp),ord)
            sumexp = sumexp + exp(lt-maxlt)
        end do
        v = ord*log2 + maxlt + log(sumexp)
    end function chisq_log_moment

    pure subroutine chisq_moments(df, ncp, orders, moms, log_values)
        real(dp), intent(in) :: df, ncp, orders(:)
        real(dp), intent(out) :: moms(size(orders))
        logical, intent(in), optional :: log_values
        logical :: lv
        integer :: i
        lv = .false.; if (present(log_values)) lv = log_values
        do i = 1, size(orders)
            moms(i) = chisq_log_moment(df,ncp,orders(i))
            if (.not. lv) moms(i) = exp(moms(i))
        end do
    end subroutine chisq_moments

    pure subroutine logchisq_cumulants(df, kappa)
        real(dp), intent(in) :: df
        real(dp), intent(out) :: kappa(:)
        integer :: r
        if (df <= 0.0_dp) then
            kappa = nan_dp(); return
        end if
        kappa(1) = digamma_dp(0.5_dp*df) + log2
        do r = 2, size(kappa)
            kappa(r) = polygamma_dp(r-1, 0.5_dp*df)
        end do
    end subroutine logchisq_cumulants

    pure subroutine lognoncentral_chisq_cumulants(df, ncp, kappa)
        real(dp), intent(in) :: df, ncp
        real(dp), intent(out) :: kappa(:)
        real(dp), allocatable :: ck(:), mk(:), moms(:)
        real(dp) :: lambda, width, logw, maxlogw, w, sumw
        integer :: j, jlo, jhi
        if (ncp < 0.0_dp .or. df < 0.0_dp .or. (df == 0.0_dp .and. ncp == 0.0_dp)) then
            kappa = nan_dp(); return
        end if
        if (ncp == 0.0_dp) then
            call logchisq_cumulants(df,kappa)
            return
        end if
        allocate(ck(size(kappa)),mk(size(kappa)),moms(size(kappa)))
        moms = 0.0_dp
        lambda = 0.5_dp*ncp
        width = 12.0_dp*sqrt(lambda+1.0_dp) + 30.0_dp
        jlo = max(0, int(floor(lambda-width)))
        jhi = max(jlo+1, int(ceiling(lambda+width)))
        maxlogw = neg_inf_dp()
        do j = jlo, jhi
            logw = -lambda + real(j,dp)*log(lambda) - log_gamma(real(j+1,dp))
            maxlogw = max(maxlogw,logw)
        end do
        sumw = 0.0_dp
        do j = jlo, jhi
            logw = -lambda + real(j,dp)*log(lambda) - log_gamma(real(j+1,dp))
            w = exp(logw-maxlogw)
            if (df+2.0_dp*real(j,dp) == 0.0_dp) cycle
            call logchisq_cumulants(df+2.0_dp*real(j,dp),ck)
            call cumulants_to_moments(ck,mk)
            moms = moms + w*mk
            sumw = sumw + w
        end do
        moms = moms/sumw
        call moments_to_cumulants(moms,kappa)
    end subroutine lognoncentral_chisq_cumulants

    real(dp) function random_normal(mu, sigma) result(x)
        real(dp), intent(in), optional :: mu, sigma
        real(dp) :: u1, u2, m, s
        m = 0.0_dp; s = 1.0_dp
        if (present(mu)) m=mu
        if (present(sigma)) s=sigma
        call random_number(u1); call random_number(u2)
        u1 = max(u1,tiny(1.0_dp))
        x = m + s*sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
    end function random_normal

    recursive real(dp) function random_gamma(shape, scale) result(x)
        real(dp), intent(in) :: shape
        real(dp), intent(in), optional :: scale
        real(dp) :: sc, d, c, z, u, v
        sc = 1.0_dp; if (present(scale)) sc=scale
        if (shape < 0.0_dp .or. sc < 0.0_dp) then
            x = nan_dp(); return
        else if (shape == 0.0_dp .or. sc == 0.0_dp) then
            x = 0.0_dp; return
        else if (shape < 1.0_dp) then
            call random_number(u)
            x = random_gamma(shape+1.0_dp,sc) * u**(1.0_dp/shape)
            return
        end if
        d = shape - 1.0_dp/3.0_dp
        c = 1.0_dp/sqrt(9.0_dp*d)
        do
            do
                z = random_normal()
                v = 1.0_dp + c*z
                if (v > 0.0_dp) exit
            end do
            v = v*v*v
            call random_number(u)
            if (u < 1.0_dp-0.0331_dp*z**4) exit
            if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
        end do
        x = sc*d*v
    end function random_gamma

    integer function random_poisson(lambda) result(k)
        real(dp), intent(in) :: lambda
        real(dp) :: l, p, u, v, slam, loglam, b, a, inv_alpha, vr, us
        if (lambda < 0.0_dp) then
            k = -1; return
        else if (lambda == 0.0_dp) then
            k = 0; return
        else if (lambda < 30.0_dp) then
            l = exp(-lambda)
            k = 0
            p = 1.0_dp
            do
                call random_number(u)
                p = p*u
                if (p <= l) exit
                k = k+1
            end do
            return
        end if
        slam = sqrt(lambda)
        loglam = log(lambda)
        b = 0.931_dp + 2.53_dp*slam
        a = -0.059_dp + 0.02483_dp*b
        inv_alpha = 1.1239_dp + 1.1328_dp/(b-3.4_dp)
        vr = 0.9277_dp - 3.6224_dp/(b-2.0_dp)
        do
            call random_number(u); call random_number(v)
            u = u - 0.5_dp
            us = 0.5_dp - abs(u)
            if (us <= 0.0_dp) cycle
            k = int(floor((2.0_dp*a/us+b)*u + lambda + 0.43_dp))
            if (us >= 0.07_dp .and. v <= vr) return
            if (k < 0) cycle
            if (us < 0.013_dp .and. v > us) cycle
            if (log(v*inv_alpha/(a/(us*us)+b)) <= &
                -lambda + real(k,dp)*loglam - log_gamma(real(k+1,dp))) return
        end do
    end function random_poisson

    real(dp) function random_chisq(df, ncp) result(x)
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: ncp
        real(dp) :: nc, shape
        integer :: j
        nc=0.0_dp; if (present(ncp)) nc=ncp
        if (df < 0.0_dp .or. nc < 0.0_dp) then
            x=nan_dp(); return
        end if
        if (nc > 0.0_dp) then
            j = random_poisson(0.5_dp*nc)
        else
            j = 0
        end if
        shape = 0.5_dp*df + real(j,dp)
        if (shape == 0.0_dp) then
            x=0.0_dp
        else
            x=random_gamma(shape,2.0_dp)
        end if
    end function random_chisq

    subroutine fill_normal(x, mu, sigma)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in), optional :: mu, sigma
        integer :: i
        do i=1,size(x)
            x(i)=random_normal(mu,sigma)
        end do
    end subroutine fill_normal

    subroutine fill_chisq(x, df, ncp)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: ncp
        real(dp) :: nc
        integer :: i
        nc=0.0_dp; if(present(ncp)) nc=ncp
        do i=1,size(x)
            x(i)=random_chisq(df,nc)
        end do
    end subroutine fill_chisq

end module sadists_special
