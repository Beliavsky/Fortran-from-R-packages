! SPDX-License-Identifier: GPL-2.0-or-later
module msm_stats
    use msm_kinds, only : dp, msm_pi
    implicit none
    private
    public :: normal_pdf, normal_cdf, lognormal_pdf, exponential_pdf, gamma_pdf
    public :: weibull_pdf, poisson_pmf, binomial_pmf, beta_binomial_pmf
    public :: negbinomial_pmf, beta_pdf, student_t_pdf, digamma_approx
    public :: truncated_normal_pdf, me_truncated_normal_pdf, me_uniform_pdf
    public :: set_random_seed, rand_uniform, rand_normal, rand_exponential, rand_gamma
    public :: logit, expit
contains
    pure function logit(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = log(x/(1.0_dp-x))
    end function logit

    pure function expit(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        if (x >= 0.0_dp) then
            y = 1.0_dp/(1.0_dp + exp(-x))
        else
            y = exp(x)/(1.0_dp + exp(x))
        end if
    end function expit

    pure function normal_pdf(x, mean, sd) result(f)
        real(dp), intent(in) :: x, mean, sd
        real(dp) :: f, z
        if (sd <= 0.0_dp) then
            f = 0.0_dp
            return
        end if
        z = (x-mean)/sd
        f = exp(-0.5_dp*z*z)/(sd*sqrt(2.0_dp*msm_pi))
    end function normal_pdf

    pure function normal_cdf(x, mean, sd) result(p)
        real(dp), intent(in) :: x, mean, sd
        real(dp) :: p
        if (sd <= 0.0_dp) then
            p = merge(1.0_dp, 0.0_dp, x >= mean)
        else
            p = 0.5_dp*erfc(-(x-mean)/(sd*sqrt(2.0_dp)))
        end if
    end function normal_cdf

    pure function lognormal_pdf(x, meanlog, sdlog) result(f)
        real(dp), intent(in) :: x, meanlog, sdlog
        real(dp) :: f
        if (x <= 0.0_dp .or. sdlog <= 0.0_dp) then
            f = 0.0_dp
        else
            f = normal_pdf(log(x), meanlog, sdlog)/x
        end if
    end function lognormal_pdf

    pure function exponential_pdf(x, rate) result(f)
        real(dp), intent(in) :: x, rate
        real(dp) :: f
        if (x < 0.0_dp .or. rate < 0.0_dp) then
            f = 0.0_dp
        else
            f = rate*exp(-rate*x)
        end if
    end function exponential_pdf

    pure function gamma_pdf(x, shape, rate) result(f)
        real(dp), intent(in) :: x, shape, rate
        real(dp) :: f
        if (x <= 0.0_dp .or. shape <= 0.0_dp .or. rate <= 0.0_dp) then
            f = 0.0_dp
        else
            f = exp(shape*log(rate) + (shape-1.0_dp)*log(x) - rate*x - log_gamma(shape))
        end if
    end function gamma_pdf

    pure function weibull_pdf(x, shape, scale) result(f)
        real(dp), intent(in) :: x, shape, scale
        real(dp) :: f, z
        if (x < 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
            f = 0.0_dp
        else if (abs(x) <= tiny(1.0_dp) .and. shape < 1.0_dp) then
            f = huge(1.0_dp)
        else
            z = x/scale
            f = shape/scale*z**(shape-1.0_dp)*exp(-(z**shape))
        end if
    end function weibull_pdf

    pure function poisson_pmf(x, lambda) result(f)
        real(dp), intent(in) :: x, lambda
        real(dp) :: f
        integer :: k
        k = nint(x)
        if (lambda < 0.0_dp .or. k < 0 .or. abs(x-real(k,dp)) > 1.0e-10_dp) then
            f = 0.0_dp
        else if (lambda <= 0.0_dp) then
            f = merge(1.0_dp, 0.0_dp, k == 0)
        else
            f = exp(real(k,dp)*log(lambda)-lambda-log_gamma(real(k+1,dp)))
        end if
    end function poisson_pmf

    pure function binomial_pmf(x, size, prob) result(f)
        real(dp), intent(in) :: x, size, prob
        real(dp) :: f, n
        integer :: k, ni
        k = nint(x); ni = nint(size); n = real(ni,dp)
        if (ni < 0 .or. k < 0 .or. k > ni .or. prob < 0.0_dp .or. prob > 1.0_dp .or. &
            abs(x-real(k,dp)) > 1.0e-10_dp) then
            f = 0.0_dp
        else if (prob <= 0.0_dp) then
            f = merge(1.0_dp,0.0_dp,k==0)
        else if (prob >= 1.0_dp) then
            f = merge(1.0_dp,0.0_dp,k==ni)
        else
            f = exp(log_gamma(n+1.0_dp)-log_gamma(real(k+1,dp))-log_gamma(n-real(k,dp)+1.0_dp) + &
                real(k,dp)*log(prob)+(n-real(k,dp))*log(1.0_dp-prob))
        end if
    end function binomial_pmf

    pure function beta_binomial_pmf(x, size, meanp, sdp) result(f)
        real(dp), intent(in) :: x, size, meanp, sdp
        real(dp) :: f, a, b, n, k
        integer :: ki, ni
        ki = nint(x); ni = nint(size); n = real(ni,dp); k = real(ki,dp)
        if (sdp <= 0.0_dp .or. meanp <= 0.0_dp .or. meanp >= 1.0_dp .or. ki < 0 .or. ki > ni) then
            f = 0.0_dp
            return
        end if
        a = meanp/sdp
        b = (1.0_dp-meanp)/sdp
        f = exp(log_gamma(n+1.0_dp)-log_gamma(k+1.0_dp)-log_gamma(n-k+1.0_dp) + &
            log_gamma(k+a)+log_gamma(n-k+b)-log_gamma(n+a+b) - &
            (log_gamma(a)+log_gamma(b)-log_gamma(a+b)))
    end function beta_binomial_pmf

    pure function negbinomial_pmf(x, size, prob) result(f)
        real(dp), intent(in) :: x, size, prob
        real(dp) :: f, k
        integer :: ki
        ki = nint(x); k = real(ki,dp)
        if (ki < 0 .or. size <= 0.0_dp .or. prob <= 0.0_dp .or. prob > 1.0_dp) then
            f = 0.0_dp
        else
            f = exp(log_gamma(k+size)-log_gamma(size)-log_gamma(k+1.0_dp) + &
                size*log(prob)+k*log(1.0_dp-prob))
        end if
    end function negbinomial_pmf

    pure function beta_pdf(x, shape1, shape2) result(f)
        real(dp), intent(in) :: x, shape1, shape2
        real(dp) :: f
        if (x <= 0.0_dp .or. x >= 1.0_dp .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp) then
            f = 0.0_dp
        else
            f = exp((shape1-1.0_dp)*log(x)+(shape2-1.0_dp)*log(1.0_dp-x) - &
                log_gamma(shape1)-log_gamma(shape2)+log_gamma(shape1+shape2))
        end if
    end function beta_pdf

    pure function student_t_pdf(x, mean, scale, df) result(f)
        real(dp), intent(in) :: x, mean, scale, df
        real(dp) :: f, z
        if (scale <= 0.0_dp .or. df <= 0.0_dp) then
            f = 0.0_dp
            return
        end if
        z = (x-mean)/scale
        f = exp(log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df)) / &
            (scale*sqrt(df*msm_pi)) * (1.0_dp+z*z/df)**(-0.5_dp*(df+1.0_dp))
    end function student_t_pdf

    pure function digamma_approx(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y, z, inv, inv2
        if (x <= 0.0_dp) then
            y = huge(1.0_dp)
            return
        end if
        z = x
        y = 0.0_dp
        do while (z < 8.0_dp)
            y = y - 1.0_dp/z
            z = z + 1.0_dp
        end do
        inv = 1.0_dp/z
        inv2 = inv*inv
        y = y + log(z) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - &
            inv2*(1.0_dp/252.0_dp - inv2*(1.0_dp/240.0_dp))))
    end function digamma_approx

    pure function truncated_normal_pdf(x, mean, sd, lower, upper) result(f)
        real(dp), intent(in) :: x, mean, sd, lower, upper
        real(dp) :: f, den
        if (x < lower .or. x > upper) then
            f = 0.0_dp
            return
        end if
        den = normal_cdf(upper,mean,sd)-normal_cdf(lower,mean,sd)
        if (den <= 0.0_dp) then
            f = 0.0_dp
        else
            f = normal_pdf(x,mean,sd)/den
        end if
    end function truncated_normal_pdf

    pure function me_truncated_normal_pdf(x, mean, sd, lower, upper, sderr, meanerr) result(f)
        real(dp), intent(in) :: x, mean, sd, lower, upper, sderr, meanerr
        real(dp) :: f, sumsq, sigtmp, mutmp, nc, nctmp
        sumsq = sd*sd+sderr*sderr
        if (sumsq <= 0.0_dp .or. sd <= 0.0_dp .or. sderr <= 0.0_dp) then
            f = 0.0_dp
            return
        end if
        sigtmp = sd*sderr/sqrt(sumsq)
        mutmp = ((x-meanerr)*sd*sd+mean*sderr*sderr)/sumsq
        nc = normal_cdf(upper,mean,sd)-normal_cdf(lower,mean,sd)
        nctmp = normal_cdf(upper,mutmp,sigtmp)-normal_cdf(lower,mutmp,sigtmp)
        f = nctmp*normal_pdf(x,meanerr+mean,sqrt(sumsq))/nc
    end function me_truncated_normal_pdf

    pure function me_uniform_pdf(x, lower, upper, sderr, meanerr) result(f)
        real(dp), intent(in) :: x, lower, upper, sderr, meanerr
        real(dp) :: f
        if (upper <= lower .or. sderr <= 0.0_dp) then
            f = 0.0_dp
        else
            f = (normal_cdf(x,meanerr+lower,sderr)-normal_cdf(x,meanerr+upper,sderr))/(upper-lower)
        end if
    end function me_uniform_pdf

    subroutine set_random_seed(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)
        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729*i + 7919*i*i, huge(1)-1)
            if (put(i) == 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine set_random_seed

    function rand_uniform() result(u)
        real(dp) :: u
        call random_number(u)
        if (u <= 0.0_dp) u = tiny(1.0_dp)
        if (u >= 1.0_dp) u = 1.0_dp-epsilon(1.0_dp)
    end function rand_uniform

    function rand_normal() result(z)
        real(dp) :: z, u1, u2
        u1 = rand_uniform(); u2 = rand_uniform()
        z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*msm_pi*u2)
    end function rand_normal

    recursive function rand_gamma(shape, rate) result(x)
        real(dp), intent(in) :: shape, rate
        real(dp) :: x, d, c, z, u
        if (shape <= 0.0_dp .or. rate <= 0.0_dp) then
            x = 0.0_dp; return
        end if
        if (shape < 1.0_dp) then
            x = rand_gamma(shape+1.0_dp,rate)*rand_uniform()**(1.0_dp/shape)
            return
        end if
        d = shape-1.0_dp/3.0_dp; c = 1.0_dp/sqrt(9.0_dp*d)
        do
            z=rand_normal(); if (1.0_dp+c*z <= 0.0_dp) cycle
            u=rand_uniform()
            if (log(u) < 0.5_dp*z*z + d*(1.0_dp-(1.0_dp+c*z)**3+3.0_dp*log(1.0_dp+c*z))) exit
        end do
        x=d*(1.0_dp+c*z)**3/rate
    end function rand_gamma

    function rand_exponential(rate) result(x)
        real(dp), intent(in) :: rate
        real(dp) :: x
        if (rate <= 0.0_dp) then
            x = huge(1.0_dp)
        else
            x = -log(rand_uniform())/rate
        end if
    end function rand_exponential
end module msm_stats
