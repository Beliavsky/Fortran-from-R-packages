! SPDX-License-Identifier: GPL-3.0-or-later
module nbconv_rng
    use nbconv_kinds, only : dp
    implicit none
    private

    integer(kind=8), parameter :: pm_mod = 2147483647_8
    integer(kind=8), parameter :: pm_mult = 16807_8
    integer(kind=8), parameter :: pm_q = 127773_8
    integer(kind=8), parameter :: pm_r = 2836_8
    integer(kind=8), save :: rng_state = 123456789_8
    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: nbconv_seed
    public :: uniform_rng
    public :: normal_rng
    public :: gamma_rng
    public :: poisson_rng
    public :: negbin_rng_mu

contains

    subroutine nbconv_seed(seed)
        integer(kind=8), intent(in) :: seed
        rng_state = modulo(abs(seed), pm_mod - 1_8) + 1_8
    end subroutine nbconv_seed

    function uniform_rng() result(u)
        real(dp) :: u
        integer(kind=8) :: hi, lo, test

        hi = rng_state / pm_q
        lo = modulo(rng_state, pm_q)
        test = pm_mult * lo - pm_r * hi
        if (test > 0_8) then
            rng_state = test
        else
            rng_state = test + pm_mod
        end if
        u = real(rng_state, dp) / real(pm_mod, dp)
    end function uniform_rng

    function normal_rng() result(z)
        real(dp) :: z
        real(dp) :: u1, u2

        u1 = max(uniform_rng(), tiny(1.0_dp))
        u2 = uniform_rng()
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
    end function normal_rng

    recursive function gamma_rng(shape, scale) result(x)
        real(dp), intent(in) :: shape
        real(dp), intent(in), optional :: scale
        real(dp) :: x
        real(dp) :: d, c, z, v, u, scl

        if (shape <= 0.0_dp) error stop "gamma_rng: shape must be positive"
        scl = 1.0_dp
        if (present(scale)) then
            if (scale <= 0.0_dp) error stop "gamma_rng: scale must be positive"
            scl = scale
        end if

        if (shape < 1.0_dp) then
            u = max(uniform_rng(), tiny(1.0_dp))
            x = gamma_rng(shape + 1.0_dp) * u**(1.0_dp / shape)
            x = x * scl
            return
        end if

        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            do
                z = normal_rng()
                v = 1.0_dp + c * z
                if (v > 0.0_dp) exit
            end do
            v = v * v * v
            u = uniform_rng()
            if (u < 1.0_dp - 0.0331_dp * z**4) exit
            if (log(max(u, tiny(1.0_dp))) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
        end do
        x = scl * d * v
    end function gamma_rng

    function poisson_rng(lambda) result(k)
        real(dp), intent(in) :: lambda
        integer :: k
        real(dp) :: l, p, u
        real(dp) :: a, b, inv_alpha, v_r, us, v, lhs, rhs
        integer :: kk

        if (lambda < 0.0_dp) error stop "poisson_rng: lambda must be nonnegative"
        if (lambda <= 0.0_dp) then
            k = 0
            return
        end if

        if (lambda < 30.0_dp) then
            l = exp(-lambda)
            p = 1.0_dp
            k = 0
            do
                k = k + 1
                p = p * uniform_rng()
                if (p <= l) exit
            end do
            k = k - 1
            return
        end if

        b = 0.931_dp + 2.53_dp * sqrt(lambda)
        a = -0.059_dp + 0.02483_dp * b
        inv_alpha = 1.1239_dp + 1.1328_dp / (b - 3.4_dp)
        v_r = 0.9277_dp - 3.6224_dp / (b - 2.0_dp)

        do
            u = uniform_rng() - 0.5_dp
            v = uniform_rng()
            us = 0.5_dp - abs(u)
            if (us <= tiny(1.0_dp)) cycle
            kk = floor((2.0_dp * a / us + b) * u + lambda + 0.43_dp)
            if (us >= 0.07_dp .and. v <= v_r .and. kk >= 0) then
                k = kk
                return
            end if
            if (kk < 0) cycle
            if (us < 0.013_dp .and. v > us) cycle
            lhs = log(max(v * inv_alpha / (a / (us * us) + b), tiny(1.0_dp)))
            rhs = -lambda + real(kk, dp) * log(lambda) - log_gamma(real(kk + 1, dp))
            if (lhs <= rhs) then
                k = kk
                return
            end if
        end do
    end function poisson_rng

    function negbin_rng_mu(mu, phi) result(k)
        real(dp), intent(in) :: mu, phi
        integer :: k
        real(dp) :: lambda

        if (mu < 0.0_dp) error stop "negbin_rng_mu: mu must be nonnegative"
        if (phi <= 0.0_dp) error stop "negbin_rng_mu: phi must be positive"
        if (mu <= 0.0_dp) then
            k = 0
            return
        end if
        lambda = gamma_rng(phi, mu / phi)
        k = poisson_rng(lambda)
    end function negbin_rng_mu

end module nbconv_rng
