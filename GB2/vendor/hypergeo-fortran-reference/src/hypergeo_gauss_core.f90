! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_gauss_core
    use hypergeo_kinds, only : dp, pi, ci
    use hypergeo_special, only : complex_gamma, complex_log_gamma, complex_digamma, &
        gamma_ratio3, gamma_ratio4, is_near_integer, is_nonpos_integer, &
        is_zero_parameter, factorial_int, pochhammer_int, nan_complex, finite_complex
    use hypergeo_generalized, only : genhypergeo, genhypergeo_series
    use hypergeo_types, only : hypergeo_info
    implicit none
    private

    complex(dp), save :: integral_a = (0.0_dp, 0.0_dp)
    complex(dp), save :: integral_b = (0.0_dp, 0.0_dp)
    complex(dp), save :: integral_c = (1.0_dp, 0.0_dp)
    complex(dp), save :: integral_z = (0.0_dp, 0.0_dp)

    public :: hypergeo_core, hypergeo_powerseries, hypergeo_general
    public :: hypergeo_taylor, hypergeo_a_nonpos_int, hypergeo_aorb_nonpos_int
    public :: f15_1_1, f15_3_1, f15_3_3, f15_3_4, f15_3_5
    public :: f15_3_6, f15_3_7, f15_3_8, f15_3_9
    public :: i15_3_6, i15_3_7, i15_3_8, i15_3_9
    public :: j15_3_6, j15_3_7, j15_3_8, j15_3_9
    public :: thingfun, crit_points
    public :: hypergeo_cover1, hypergeo_cover2, hypergeo_cover3
    public :: f15_3_10, f15_3_10_a, f15_3_10_b
    public :: f15_3_11, f15_3_11_bit1, f15_3_11_bit2_a, f15_3_11_bit2_b
    public :: f15_3_12, f15_3_12_bit1, f15_3_12_bit2_a, f15_3_12_bit2_b
    public :: f15_3_13, f15_3_13_a, f15_3_13_b
    public :: f15_3_14, f15_3_14_bit1_a, f15_3_14_bit1_b, f15_3_14_bit2
    public :: f15_3_10_11_12, f15_3_13_14
    public :: w07_23_06_0029_01, w07_23_06_0031_01, w07_23_06_0026_01
    public :: w07_23_06_0031_01_bit1, w07_23_06_0031_01_bit2
    public :: w07_23_06_0026_01_bit1, w07_23_06_0026_01_bit2
    public :: w07_23_06_0026_01_bit3_a, w07_23_06_0026_01_bit3_b, w07_23_06_0026_01_bit3_c
    public :: hypergeo_gosper

contains

    pure function crit_points() result(c)
        complex(dp) :: c(2)
        c(1) = cmplx(0.5_dp, sqrt(3.0_dp) / 2.0_dp, kind=dp)
        c(2) = conjg(c(1))
    end function crit_points

    pure function thingfun(z) result(v)
        complex(dp), intent(in) :: z
        complex(dp) :: v(6)
        v(1) = z
        v(2) = z / (z - 1.0_dp)
        v(3) = 1.0_dp - z
        v(4) = 1.0_dp / z
        v(5) = 1.0_dp / (1.0_dp - z)
        v(6) = 1.0_dp - 1.0_dp / z
    end function thingfun

    function f15_1_1(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = genhypergeo([a, b], [c], z, tol=tol, maxiter=maxiter)
    end function f15_1_1

    function f15_3_1(a, b, c, z, h) result(value)
        complex(dp), intent(in) :: a, b, c, z
        complex(dp), intent(in), optional :: h
        complex(dp) :: value, mult, hp

        mult = gamma_ratio3(c, b, c - b)
        integral_a = a
        integral_b = b
        integral_c = c
        integral_z = z
        if (.not. present(h)) then
            value = mult * tanh_sinh_segment((0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp))
        else
            hp = h
            value = mult * (tanh_sinh_segment((0.0_dp, 0.0_dp), hp) &
                + tanh_sinh_segment(hp, (1.0_dp, 0.0_dp)))
        end if
    end function f15_3_1

    function tanh_sinh_segment(z0, z1) result(value)
        complex(dp), intent(in) :: z0, z1
        complex(dp) :: value, previous
        real(dp) :: step
        integer :: level

        previous = (0.0_dp, 0.0_dp)
        step = 0.25_dp
        do level = 1, 5
            value = tanh_sinh_sum(z0, z1, step)
            if (level > 1) then
                if (abs(value - previous) <= 5.0e-13_dp * max(1.0_dp, abs(value))) return
            end if
            previous = value
            step = 0.5_dp * step
        end do
    end function tanh_sinh_segment

    function tanh_sinh_sum(z0, z1, step) result(value)
        complex(dp), intent(in) :: z0, z1
        real(dp), intent(in) :: step
        complex(dp) :: value, point, one_minus_point, one_minus_zpoint, dz, term
        real(dp) :: u, sh, ch, v, x, omx, q, sech2, weight
        integer :: k, kmax

        value = (0.0_dp, 0.0_dp)
        dz = z1 - z0
        kmax = ceiling(4.0_dp / step)
        do k = -kmax, kmax
            u = real(k, dp) * step
            sh = sinh(u)
            ch = cosh(u)
            v = 0.5_dp * pi * sh
            if (v >= 0.0_dp) then
                q = exp(-2.0_dp * v)
                x = 1.0_dp / (1.0_dp + q)
                omx = q / (1.0_dp + q)
            else
                q = exp(2.0_dp * v)
                x = q / (1.0_dp + q)
                omx = 1.0_dp / (1.0_dp + q)
            end if
            sech2 = 4.0_dp * q / ((1.0_dp + q) * (1.0_dp + q))
            weight = 0.25_dp * pi * ch * sech2
            if (.not. (weight > tiny(1.0_dp))) cycle

            ! Form both t and 1-t independently.  This avoids cancellation near
            ! the algebraically singular endpoints even after x rounds to 0 or 1.
            point = z0 * omx + z1 * x
            one_minus_point = (1.0_dp - z0) * omx + (1.0_dp - z1) * x
            one_minus_zpoint = (1.0_dp - integral_z * z0) * omx &
                + (1.0_dp - integral_z * z1) * x
            term = point ** (integral_b - 1.0_dp) &
                * one_minus_point ** (integral_c - integral_b - 1.0_dp) &
                * one_minus_zpoint ** (-integral_a) * dz * weight
            if (finite_complex(term)) value = value + term
        end do
        value = value * step
    end function tanh_sinh_sum

    function f15_3_3(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = (1.0_dp - z) ** (c - a - b) &
            * genhypergeo([c - a, c - b], [c], z, tol=tol, maxiter=maxiter)
    end function f15_3_3

    function f15_3_4(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = (1.0_dp - z) ** (-a) &
            * genhypergeo([a, c - b], [c], z / (z - 1.0_dp), tol=tol, maxiter=maxiter)
    end function f15_3_4

    function f15_3_5(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = (1.0_dp - z) ** (-b) &
            * genhypergeo([b, c - a], [c], z / (z - 1.0_dp), tol=tol, maxiter=maxiter)
    end function f15_3_5

    pure function i15_3_6(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        complex(dp) :: v(2)
        if (is_nonpos_integer(c - a) .or. is_nonpos_integer(c - b)) then
            v(1) = 0.0_dp
        else
            v(1) = gamma_ratio4(c, c - a - b, c - a, c - b)
        end if
        if (is_nonpos_integer(a) .or. is_nonpos_integer(b)) then
            v(2) = 0.0_dp
        else
            v(2) = gamma_ratio4(c, a + b - c, a, b)
        end if
    end function i15_3_6

    pure function j15_3_6(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        logical :: v(4)
        v = [is_nonpos_integer(c), is_nonpos_integer(c - a - b), &
             is_nonpos_integer(c), is_nonpos_integer(a + b - c)]
    end function j15_3_6

    function f15_3_6(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, coef(2)
        coef = i15_3_6(a, b, c)
        value = coef(1) * genhypergeo([a, b], [a + b - c + 1.0_dp], 1.0_dp - z, &
            tol=tol, maxiter=maxiter) &
            + coef(2) * genhypergeo([c - a, c - b], [c - a - b + 1.0_dp], 1.0_dp - z, &
            tol=tol, maxiter=maxiter) * (1.0_dp - z) ** (c - a - b)
    end function f15_3_6

    pure function i15_3_7(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        complex(dp) :: v(2)
        if (is_nonpos_integer(b) .or. is_nonpos_integer(c - a)) then
            v(1) = 0.0_dp
        else
            v(1) = gamma_ratio4(c, b - a, b, c - a)
        end if
        if (is_nonpos_integer(a) .or. is_nonpos_integer(c - b)) then
            v(2) = 0.0_dp
        else
            v(2) = gamma_ratio4(c, a - b, a, c - b)
        end if
    end function i15_3_7

    pure function j15_3_7(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        logical :: v(4)
        v = [is_nonpos_integer(c), is_nonpos_integer(b - a), &
             is_nonpos_integer(c), is_nonpos_integer(a - b)]
    end function j15_3_7

    function f15_3_7(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, coef(2)
        coef = i15_3_7(a, b, c)
        value = coef(1) * (-z) ** (-a) &
            * genhypergeo([a, 1.0_dp - c + a], [1.0_dp - b + a], 1.0_dp / z, &
                tol=tol, maxiter=maxiter) &
            + coef(2) * (-z) ** (-b) &
            * genhypergeo([b, 1.0_dp - c + b], [1.0_dp - a + b], 1.0_dp / z, &
                tol=tol, maxiter=maxiter)
    end function f15_3_7

    pure function i15_3_8(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        complex(dp) :: v(2)
        v = i15_3_7(a, b, c)
    end function i15_3_8

    pure function j15_3_8(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        logical :: v(4)
        v = j15_3_7(a, b, c)
    end function j15_3_8

    function f15_3_8(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, coef(2)
        coef = i15_3_8(a, b, c)
        value = coef(1) * (1.0_dp - z) ** (-a) &
            * genhypergeo([a, c - b], [a - b + 1.0_dp], 1.0_dp / (1.0_dp - z), &
                tol=tol, maxiter=maxiter) &
            + coef(2) * (1.0_dp - z) ** (-b) &
            * genhypergeo([b, c - a], [b - a + 1.0_dp], 1.0_dp / (1.0_dp - z), &
                tol=tol, maxiter=maxiter)
    end function f15_3_8

    pure function i15_3_9(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        complex(dp) :: v(2)
        v = i15_3_6(a, b, c)
    end function i15_3_9

    pure function j15_3_9(a, b, c) result(v)
        complex(dp), intent(in) :: a, b, c
        logical :: v(4)
        v = j15_3_6(a, b, c)
    end function j15_3_9

    function f15_3_9(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, coef(2)
        coef = i15_3_9(a, b, c)
        value = coef(1) * z ** (-a) &
            * genhypergeo([a, a - c + 1.0_dp], [a + b - c + 1.0_dp], 1.0_dp - 1.0_dp / z, &
                tol=tol, maxiter=maxiter) &
            + coef(2) * (1.0_dp - z) ** (c - a - b) * z ** (a - c) &
            * genhypergeo([c - a, 1.0_dp - a], [c - a - b + 1.0_dp], 1.0_dp - 1.0_dp / z, &
                tol=tol, maxiter=maxiter)
    end function f15_3_9

    function f15_3_10(a, b, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, b, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value
        complex(dp) :: ua, ub, term, old, pn, pa, pb, l1mz
        real(dp) :: t
        integer :: n, nmax

        t = effective_tol(tol)
        nmax = effective_maxiter(maxiter)
        ua = a
        ub = b
        term = 1.0_dp
        pn = complex_digamma(cmplx(1.0_dp, 0.0_dp, kind=dp))
        pa = complex_digamma(a)
        pb = complex_digamma(b)
        l1mz = log(1.0_dp - z)
        value = 2.0_dp * pn - pa - pb - l1mz
        do n = 1, nmax
            old = value
            term = term * ua * ub * (1.0_dp - z) / real(n * n, dp)
            pn = pn + 1.0_dp / real(n, dp)
            pa = pa + 1.0_dp / (a + real(n - 1, dp))
            pb = pb + 1.0_dp / (b + real(n - 1, dp))
            value = old + term * (2.0_dp * pn - pa - pb - l1mz)
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
        end do
        value = value * gamma_ratio3(a + b, a, b)
        if (present(method)) continue
    end function f15_3_10

    function f15_3_11(a, b, m_in, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, b, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value
        integer :: m
        m = nint(real(m_in, dp))
        if (m <= 0) then
            value = nan_complex()
            return
        end if
        value = f15_3_11_bit1_impl(a, b, m, z, tol) - &
            f15_3_11_bit2_impl(a, b, m, z, tol, maxiter)
        if (present(method)) continue
    end function f15_3_11

    function f15_3_11_bit1_impl(a, b, m, z, tol) result(value)
        complex(dp), intent(in) :: a, b, z
        integer, intent(in) :: m
        real(dp), intent(in), optional :: tol
        complex(dp) :: value, ua, ub, lower, term, old, mult
        integer :: n
        real(dp) :: t
        t = effective_tol(tol)
        mult = gamma_ratio4(cmplx(real(m, dp), 0.0_dp, dp), a + b + real(m, dp), &
            a + real(m, dp), b + real(m, dp))
        ua = a
        ub = b
        lower = cmplx(1.0_dp - real(m, dp), 0.0_dp, kind=dp)
        term = 1.0_dp
        value = 1.0_dp
        do n = 1, m - 1
            old = value
            term = term * ua * ub / lower * (1.0_dp - z) / real(n, dp)
            value = old + term
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
            lower = lower + 1.0_dp
        end do
        value = value * mult
    end function f15_3_11_bit1_impl

    function f15_3_11_bit2_impl(a, b, m, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, z
        integer, intent(in) :: m
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, ua, ub, term, old, pn, pm, pa, pb, l1mz, mult
        integer :: n, nmax
        real(dp) :: t, parity
        t = effective_tol(tol)
        nmax = effective_maxiter(maxiter)
        parity = parity_sign(m)
        mult = parity * (z - 1.0_dp) ** real(m, dp) * gamma_ratio3(a + b + real(m, dp), a, b)
        ua = a + real(m, dp)
        ub = b + real(m, dp)
        term = 1.0_dp / factorial_int(m)
        pn = complex_digamma(cmplx(1.0_dp, 0.0_dp, dp))
        pm = complex_digamma(cmplx(real(m + 1, dp), 0.0_dp, dp))
        pa = complex_digamma(a + real(m, dp))
        pb = complex_digamma(b + real(m, dp))
        l1mz = log(1.0_dp - z)
        value = term * (l1mz - pn - pm + pa + pb)
        do n = 1, nmax
            old = value
            term = term * ua * ub * (1.0_dp - z) / real(n * (n + m), dp)
            pn = pn + 1.0_dp / real(n, dp)
            pm = pm + 1.0_dp / real(n + m, dp)
            pa = pa + 1.0_dp / (a + real(n + m - 1, dp))
            pb = pb + 1.0_dp / (b + real(n + m - 1, dp))
            value = old + term * (l1mz - pn - pm + pa + pb)
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
        end do
        value = mult / parity * value
        ! Upstream has (z-1)^m rather than (-1)^m(1-z)^m; mult/parity keeps that form.
    end function f15_3_11_bit2_impl

    function f15_3_12(a, b, m_in, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, b, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value
        integer :: m
        m = nint(real(m_in, dp))
        if (m <= 0) then
            value = nan_complex()
            return
        end if
        value = f15_3_12_bit1_impl(a, b, m, z, tol) - &
            f15_3_12_bit2_impl(a, b, m, z, tol, maxiter)
        if (present(method)) continue
    end function f15_3_12

    function f15_3_12_bit1_impl(a, b, m, z, tol) result(value)
        complex(dp), intent(in) :: a, b, z
        integer, intent(in) :: m
        real(dp), intent(in), optional :: tol
        complex(dp) :: value, ua, ub, lower, term, old, mult
        integer :: n
        real(dp) :: t
        t = effective_tol(tol)
        mult = gamma_ratio4(cmplx(real(m, dp), 0.0_dp, dp), a + b - real(m, dp), a, b) &
            / (1.0_dp - z) ** real(m, dp)
        ua = a - real(m, dp)
        ub = b - real(m, dp)
        lower = cmplx(1.0_dp - real(m, dp), 0.0_dp, dp)
        term = 1.0_dp
        value = 1.0_dp
        do n = 1, m - 1
            old = value
            term = term * ua * ub / lower * (1.0_dp - z) / real(n, dp)
            value = old + term
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
            lower = lower + 1.0_dp
        end do
        value = value * mult
    end function f15_3_12_bit1_impl

    function f15_3_12_bit2_impl(a, b, m, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, z
        integer, intent(in) :: m
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, ua, ub, term, old, pn, pm, pa, pb, l1mz, mult
        integer :: n, nmax
        real(dp) :: t
        if (is_nonpos_integer(a - real(m, dp)) .or. is_nonpos_integer(b - real(m, dp))) then
            value = (0.0_dp, 0.0_dp)
            return
        end if
        t = effective_tol(tol)
        nmax = effective_maxiter(maxiter)
        mult = parity_sign(m) * gamma_ratio3(a + b - real(m, dp), a - real(m, dp), b - real(m, dp))
        ua = a
        ub = b
        term = 1.0_dp / factorial_int(m)
        pn = complex_digamma(cmplx(1.0_dp, 0.0_dp, dp))
        pm = complex_digamma(cmplx(real(m + 1, dp), 0.0_dp, dp))
        pa = complex_digamma(a)
        pb = complex_digamma(b)
        l1mz = log(1.0_dp - z)
        value = term * (l1mz - pn - pm + pa + pb)
        do n = 1, nmax
            old = value
            term = term * ua * ub * (1.0_dp - z) / real(n * (n + m), dp)
            pn = pn + 1.0_dp / real(n, dp)
            pm = pm + 1.0_dp / real(n + m, dp)
            pa = pa + 1.0_dp / (a + real(n - 1, dp))
            pb = pb + 1.0_dp / (b + real(n - 1, dp))
            value = old + term * (l1mz - pn - pm + pa + pb)
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
        end do
        value = mult * value
    end function f15_3_12_bit2_impl

    function f15_3_13(a, c, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value, ua, ub, term, old, pn, pa, pc, lmz, mult
        integer :: n, nmax
        real(dp) :: t
        t = effective_tol(tol)
        nmax = effective_maxiter(maxiter)
        ua = a
        ub = 1.0_dp - c + a
        term = 1.0_dp
        pn = complex_digamma(cmplx(1.0_dp, 0.0_dp, dp))
        pa = complex_digamma(a)
        pc = complex_digamma(c - a)
        lmz = log(-z)
        value = lmz + 2.0_dp * pn - pa - pc
        do n = 1, nmax
            old = value
            term = term * ua * ub / (z * real(n * n, dp))
            pn = pn + 1.0_dp / real(n, dp)
            pa = pa + 1.0_dp / (a + real(n - 1, dp))
            pc = pc - 1.0_dp / (c - a - real(n, dp))
            value = old + term * (lmz + 2.0_dp * pn - pa - pc)
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
        end do
        mult = gamma_ratio3(c, a, c - a) * (-z) ** (-a)
        value = mult * value
        if (present(method)) continue
    end function f15_3_13

    function f15_3_14(a, c, m_in, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, c, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value
        integer :: m
        m = nint(real(m_in, dp))
        if (m <= 0) then
            value = nan_complex()
            return
        end if
        value = f15_3_14_bit1_impl(a, c, m, z, tol, maxiter) &
            + f15_3_14_bit2_impl(a, c, m, z, tol)
        if (present(method)) continue
    end function f15_3_14

    function f15_3_14_bit1_impl(a, c, m, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, c, z
        integer, intent(in) :: m
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, ua, ub, term, old, pm, pn, pa, pc, lmz, mult
        integer :: n, nmax
        real(dp) :: t
        t = effective_tol(tol)
        nmax = effective_maxiter(maxiter)
        ua = a + real(m, dp)
        ub = 1.0_dp - c + a + real(m, dp)
        term = exp(complex_log_gamma(a + real(m, dp)) - complex_log_gamma(a) &
            + complex_log_gamma(1.0_dp - c + a + real(m, dp)) &
            - complex_log_gamma(1.0_dp - c + a) &
            - complex_log_gamma(cmplx(real(m + 1, dp), 0.0_dp, dp)))
        pm = complex_digamma(cmplx(real(m + 1, dp), 0.0_dp, dp))
        pn = complex_digamma(cmplx(1.0_dp, 0.0_dp, dp))
        pa = complex_digamma(a + real(m, dp))
        pc = complex_digamma(c - a - real(m, dp))
        lmz = log(-z)
        value = term * (lmz + pm + pn - pa - pc)
        do n = 1, nmax
            old = value
            term = term * ua * ub / (z * real(n * (n + m), dp))
            pm = pm + 1.0_dp / real(n + m, dp)
            pn = pn + 1.0_dp / real(n, dp)
            pa = pa + 1.0_dp / (a + real(m + n - 1, dp))
            pc = pc - 1.0_dp / (c - a - real(m + n, dp))
            value = old + term * (lmz + pm + pn - pa - pc)
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
        end do
        mult = (-z) ** (-a - real(m, dp)) * gamma_ratio3(c, a + real(m, dp), c - a)
        value = mult * value
    end function f15_3_14_bit1_impl

    function f15_3_14_bit2_impl(a, c, m, z, tol) result(value)
        complex(dp), intent(in) :: a, c, z
        integer, intent(in) :: m
        real(dp), intent(in), optional :: tol
        complex(dp) :: value, ua, term, old, mult
        integer :: n
        real(dp) :: t
        t = effective_tol(tol)
        mult = (-z) ** (-a) * gamma_ratio3(c, a + real(m, dp), cmplx(1.0_dp, 0.0_dp, dp))
        ua = a
        term = 1.0_dp
        value = gamma_ratio3(cmplx(real(m, dp), 0.0_dp, dp), c - a, &
            cmplx(1.0_dp, 0.0_dp, dp))
        do n = 1, m - 1
            old = value
            term = term * ua / (z * real(n, dp))
            value = old + term * gamma_ratio3(cmplx(real(m - n, dp), 0.0_dp, dp), &
                c - a - real(n, dp), cmplx(1.0_dp, 0.0_dp, dp))
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
        end do
        value = value * mult
    end function f15_3_14_bit2_impl

    function f15_3_10_11_12(a, b, m_in, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, b, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value
        integer :: m
        m = nint(real(m_in, dp))
        if (m == 0) then
            value = f15_3_10(a, b, z, tol=tol, maxiter=maxiter, method=method)
        else if (m > 0) then
            value = f15_3_11(a, b, cmplx(real(m, dp), 0.0_dp, dp), z, &
                tol=tol, maxiter=maxiter, method=method)
        else
            value = f15_3_12(a, b, cmplx(real(-m, dp), 0.0_dp, dp), z, &
                tol=tol, maxiter=maxiter, method=method)
        end if
    end function f15_3_10_11_12

    function f15_3_13_14(a, c, m_in, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, c, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value
        integer :: m
        m = nint(real(m_in, dp))
        if (m == 0) then
            value = f15_3_13(a, c, z, tol=tol, maxiter=maxiter, method=method)
        else if (m > 0) then
            value = f15_3_14(a, c, cmplx(real(m, dp), 0.0_dp, dp), z, &
                tol=tol, maxiter=maxiter, method=method)
        else
            value = f15_3_14(a + real(m, dp), c, cmplx(real(-m, dp), 0.0_dp, dp), z, &
                tol=tol, maxiter=maxiter, method=method)
        end if
    end function f15_3_13_14

    function w07_23_06_0029_01(a, n_in, m_in, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, n_in, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        integer :: n, m
        n = nint(real(n_in, dp))
        m = nint(real(m_in, dp))
        value = parity_sign(m) * complex_gamma(a - real(m, dp)) * factorial_int(m + n) &
            * (-z) ** (-a - real(n, dp)) / (complex_gamma(a) * factorial_int(n)) &
            * hypergeo_core(a + real(n, dp), cmplx(real(m + n + 1, dp), 0.0_dp, dp), &
                cmplx(real(n + 1, dp), 0.0_dp, dp), 1.0_dp / z, tol=tol, maxiter=maxiter)
    end function w07_23_06_0029_01

    function w07_23_06_0031_01(a, n_in, m_in, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, n_in, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        integer :: n, m
        n = nint(real(n_in, dp))
        m = nint(real(m_in, dp))
        value = w0031_bit1(a, n, m, z, tol) + w0031_bit2(a, n, m, z, tol, maxiter)
    end function w07_23_06_0031_01

    function w0031_bit1(a, n, m, z, tol) result(value)
        complex(dp), intent(in) :: a, z
        integer, intent(in) :: n, m
        real(dp), intent(in), optional :: tol
        complex(dp) :: value, ua, ub, lower, term, old, mult
        integer :: k
        real(dp) :: t
        t = effective_tol(tol)
        mult = gamma_ratio4(a + real(m, dp), cmplx(real(n, dp), 0.0_dp, dp), &
            cmplx(real(m, dp), 0.0_dp, dp), a + real(n, dp)) * (-z) ** (-a)
        ua = a
        ub = cmplx(1.0_dp - real(m, dp), 0.0_dp, dp)
        lower = cmplx(1.0_dp - real(n, dp), 0.0_dp, dp)
        term = 1.0_dp
        value = 1.0_dp
        do k = 1, m - 1
            old = value
            term = term * ua * ub / lower / (real(k, dp) * z)
            value = old + term
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
            lower = lower + 1.0_dp
        end do
        value = value * mult
    end function w0031_bit1

    function w0031_bit2(a, n, m, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, z
        integer, intent(in) :: n, m
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = parity_sign(m) * (-z) ** (-a - real(n, dp)) &
            * gamma_ratio4(a + real(m, dp), cmplx(real(n - m + 1, dp), 0.0_dp, dp), &
                a, cmplx(real(n + 1, dp), 0.0_dp, dp)) &
            * hypergeo_core(a + real(n, dp), cmplx(real(1 - m + n, dp), 0.0_dp, dp), &
                cmplx(real(n + 1, dp), 0.0_dp, dp), 1.0_dp / z, tol=tol, maxiter=maxiter)
    end function w0031_bit2

    function w07_23_06_0026_01(a, n_in, m_in, z, tol, maxiter, method) result(value)
        complex(dp), intent(in) :: a, n_in, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        complex(dp) :: value
        integer :: n, m
        n = nint(real(n_in, dp))
        m = nint(real(m_in, dp))
        if (m < n .or. m < 0 .or. n < 0) then
            value = nan_complex()
            return
        end if
        value = w0026_bit1(a, n, m, z, tol) + w0026_bit2(a, n, m, z, tol, maxiter) &
            + w0026_bit3(a, n, m, z)
        if (present(method)) continue
    end function w07_23_06_0026_01

    function w0026_bit1(a, n, m, z, tol) result(value)
        complex(dp), intent(in) :: a, z
        integer, intent(in) :: n, m
        real(dp), intent(in), optional :: tol
        complex(dp) :: value, ua, ub, lower, term, old, mult
        integer :: k
        real(dp) :: t
        if (n == 0) then
            value = (0.0_dp, 0.0_dp)
            return
        end if
        t = effective_tol(tol)
        mult = (-z) ** (-a) * gamma_ratio4(cmplx(real(n, dp), 0.0_dp, dp), a + real(m, dp), &
            cmplx(real(m, dp), 0.0_dp, dp), a + real(n, dp))
        ua = a
        ub = cmplx(1.0_dp - real(m, dp), 0.0_dp, dp)
        lower = cmplx(1.0_dp - real(n, dp), 0.0_dp, dp)
        term = 1.0_dp
        value = 1.0_dp
        do k = 1, n - 1
            old = value
            term = term * ua * ub / lower / (z * real(k, dp))
            value = old + term
            if (good_delta(value - old, value, t)) exit
            ua = ua + 1.0_dp
            ub = ub + 1.0_dp
            lower = lower + 1.0_dp
        end do
        value = value * mult
    end function w0026_bit1

    function w0026_bit2(a, n, m, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, z
        integer, intent(in) :: n, m
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, mult
        mult = parity_sign(n) * (-z) ** (-a - real(m, dp)) &
            * gamma_ratio3(a + real(m, dp), a, a + real(n, dp)) &
            * gamma_ratio3(a + real(m, dp), cmplx(real(m + 1, dp), 0.0_dp, dp), &
                cmplx(real(m - n + 1, dp), 0.0_dp, dp))
        value = mult * genhypergeo([cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            a + real(m, dp)], [cmplx(real(m + 1, dp), 0.0_dp, dp), &
            cmplx(real(m - n + 1, dp), 0.0_dp, dp)], 1.0_dp / z, tol=tol, maxiter=maxiter)
    end function w0026_bit2

    function w0026_bit3(a, n, m, z) result(value)
        complex(dp), intent(in) :: a, z
        integer, intent(in) :: n, m
        complex(dp) :: value, sumv, term, mult
        integer :: k
        if (m - n <= 0) then
            value = (0.0_dp, 0.0_dp)
            return
        end if
        mult = parity_sign(n) * gamma_ratio3(a + real(m, dp), a, &
            cmplx(real(m - n, dp), 0.0_dp, dp)) * (-z) ** (-a - real(n, dp))
        sumv = (0.0_dp, 0.0_dp)
        do k = 0, m - n - 1
            term = pochhammer_int(a + real(n, dp), k) &
                * pochhammer_int(cmplx(real(1 - m + n, dp), 0.0_dp, dp), k) &
                / (factorial_int(k) * factorial_int(k + n))
            term = term * (log(-z) &
                - complex_digamma(cmplx(real(m - n - k, dp), 0.0_dp, dp)) &
                + complex_digamma(cmplx(real(k + 1, dp), 0.0_dp, dp)) &
                + complex_digamma(cmplx(real(k + n + 1, dp), 0.0_dp, dp)) &
                - complex_digamma(a + real(k + n, dp))) * z ** (-real(k, dp))
            sumv = sumv + term
        end do
        value = mult * sumv
    end function w0026_bit3

    function hypergeo_taylor(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = genhypergeo_series([a, b], [c], z, tol=tol, maxiter=maxiter, &
            check_mod=.false., polynomial=.true.)
    end function hypergeo_taylor

    function hypergeo_a_nonpos_int(a, b, c, z, tol) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        complex(dp) :: value
        complex(dp), allocatable :: empty(:)
        integer :: degree
        if (.not. is_nonpos_integer(a)) then
            value = nan_complex()
            return
        end if
        if (is_near_integer(c) .and. abs(c - a) < 0.5_dp) then
            allocate(empty(0))
            value = genhypergeo_series([b], empty, z, tol=tol, check_mod=.false.)
        else
            degree = max(0, -nint(real(a, dp)))
            value = hypergeo_taylor(a, b, c, z, tol=tol, maxiter=degree)
        end if
    end function hypergeo_a_nonpos_int

    function hypergeo_aorb_nonpos_int(a, b, c, z, tol) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        complex(dp) :: value
        if (is_nonpos_integer(a) .and. is_nonpos_integer(b)) then
            if (real(a, dp) > real(b, dp)) then
                value = hypergeo_a_nonpos_int(a, b, c, z, tol)
            else
                value = hypergeo_a_nonpos_int(b, a, c, z, tol)
            end if
        else if (is_nonpos_integer(a)) then
            value = hypergeo_a_nonpos_int(a, b, c, z, tol)
        else if (is_nonpos_integer(b)) then
            value = hypergeo_a_nonpos_int(b, a, c, z, tol)
        else
            value = nan_complex()
        end if
    end function hypergeo_aorb_nonpos_int

    function hypergeo_general(a, b, c, z, tol, maxiter, choice_out) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        integer, intent(out), optional :: choice_out
        complex(dp) :: value, t(6)
        real(dp) :: mods(6)
        integer :: choice, i

        value = nan_complex()
        if (abs(z) <= tiny(1.0_dp)) then
            value = (1.0_dp, 0.0_dp)
            if (present(choice_out)) choice_out = 1
            return
        end if
        if (abs(1.0_dp - z) <= tiny(1.0_dp)) then
            ! Use the 1-z connection formula if finite.
            value = f15_3_6(a, b, c, z, tol, maxiter)
            if (present(choice_out)) choice_out = 3
            return
        end if
        t = thingfun(z)
        do i = 1, 6
            if (finite_complex(t(i))) then
                mods(i) = abs(t(i))
            else
                mods(i) = huge(1.0_dp)
            end if
        end do
        choice = minloc(mods, dim=1)
        select case (choice)
        case (1)
            value = f15_1_1(a, b, c, z, tol, maxiter)
        case (2)
            value = f15_3_4(a, b, c, z, tol, maxiter)
        case (3)
            value = f15_3_6(a, b, c, z, tol, maxiter)
        case (4)
            value = f15_3_7(a, b, c, z, tol, maxiter)
        case (5)
            value = f15_3_8(a, b, c, z, tol, maxiter)
        case (6)
            value = f15_3_9(a, b, c, z, tol, maxiter)
        end select
        if (present(choice_out)) choice_out = choice
    end function hypergeo_general

    function hypergeo_cover1(a, b, m_in, z, tol, maxiter, method, choice_out) result(value)
        complex(dp), intent(in) :: a, b, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        integer, intent(out), optional :: choice_out
        complex(dp) :: value, c, trans(6)
        real(dp) :: mods(6)
        logical :: bad(4)
        integer :: i, choice
        value = nan_complex()
        c = a + b + m_in
        trans = thingfun(z)
        do i = 1, 6
            mods(i) = merge(abs(trans(i)), huge(1.0_dp), finite_complex(trans(i)))
        end do
        bad = j15_3_7(a, b, c)
        if (any(bad)) mods(4) = huge(1.0_dp)
        bad = j15_3_8(a, b, c)
        if (any(bad)) mods(5) = huge(1.0_dp)
        bad = j15_3_9(a, b, c)
        if (any(bad)) mods(6) = huge(1.0_dp)
        choice = minloc(mods, dim=1)
        select case(choice)
        case(1); value = f15_1_1(a, b, c, z, tol, maxiter)
        case(2); value = f15_3_4(a, b, c, z, tol, maxiter)
        case(3); value = f15_3_10_11_12(a, b, m_in, z, tol, maxiter, method)
        case(4); value = f15_3_7(a, b, c, z, tol, maxiter)
        case(5); value = f15_3_8(a, b, c, z, tol, maxiter)
        case(6); value = f15_3_9(a, b, c, z, tol, maxiter)
        end select
        if (present(choice_out)) choice_out = choice
    end function hypergeo_cover1

    function hypergeo_cover2(a, c, m_in, z, tol, maxiter, method, choice_out) result(value)
        complex(dp), intent(in) :: a, c, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        integer, intent(out), optional :: choice_out
        complex(dp) :: value, b, trans(6)
        real(dp) :: mods(6)
        logical :: bad(4)
        integer :: i, choice
        value = nan_complex()
        b = a + m_in
        trans = thingfun(z)
        do i = 1, 6
            mods(i) = merge(abs(trans(i)), huge(1.0_dp), finite_complex(trans(i)))
        end do
        bad = j15_3_6(a, b, c)
        if (any(bad)) mods(3) = huge(1.0_dp)
        bad = j15_3_8(a, b, c)
        if (any(bad)) mods(5) = huge(1.0_dp)
        bad = j15_3_9(a, b, c)
        if (any(bad)) mods(6) = huge(1.0_dp)
        choice = minloc(mods, dim=1)
        select case(choice)
        case(1); value = f15_1_1(a, b, c, z, tol, maxiter)
        case(2); value = f15_3_4(a, b, c, z, tol, maxiter)
        case(3); value = f15_3_6(a, b, c, z, tol, maxiter)
        case(4); value = f15_3_13_14(a, c, m_in, z, tol, maxiter, method)
        case(5); value = f15_3_8(a, b, c, z, tol, maxiter)
        case(6); value = f15_3_9(a, b, c, z, tol, maxiter)
        end select
        if (present(choice_out)) choice_out = choice
    end function hypergeo_cover2

    function hypergeo_cover3(a, n_in, m_in, z, tol, maxiter, method, choice_out) result(value)
        complex(dp), intent(in) :: a, n_in, m_in, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        character(len=*), intent(in), optional :: method
        integer, intent(out), optional :: choice_out
        complex(dp) :: value
        integer :: n, m
        n = nint(real(n_in, dp))
        m = nint(real(m_in, dp))
        if (abs(z) <= 1.0_dp) then
            value = f15_1_1(a, a + real(n, dp), a + real(m, dp), z, tol, maxiter)
            if (present(choice_out)) choice_out = 1
        else if (m > n) then
            value = w07_23_06_0026_01(a, cmplx(real(n, dp), 0.0_dp, dp), &
                cmplx(real(m, dp), 0.0_dp, dp), z, tol, maxiter, method)
            if (present(choice_out)) choice_out = 2
        else
            value = w07_23_06_0031_01(a, cmplx(real(n, dp), 0.0_dp, dp), &
                cmplx(real(m, dp), 0.0_dp, dp), z, tol, maxiter)
            if (present(choice_out)) choice_out = 2
        end if
    end function hypergeo_cover3

    recursive function hypergeo_powerseries(a_in, b_in, c, z, tol, maxiter, depth) result(value)
        complex(dp), intent(in) :: a_in, b_in, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter, depth
        complex(dp) :: value, a, b, swap, m, n
        integer :: d

        d = 0
        if (present(depth)) d = depth
        if (d > 12) then
            value = nan_complex()
            return
        end if
        a = a_in
        b = b_in
        if (is_zero_parameter(a) .or. is_zero_parameter(b)) then
            if (is_zero_parameter(c)) then
                value = nan_complex()
            else
                value = (1.0_dp, 0.0_dp)
            end if
            return
        end if
        if (is_zero_parameter(c)) then
            value = nan_complex()
            return
        end if
        if (is_zero_parameter(a - c)) then
            value = (1.0_dp - z) ** (-b)
            return
        else if (is_zero_parameter(b - c)) then
            value = (1.0_dp - z) ** (-a)
            return
        end if
        if (is_nonpos_integer(a) .or. is_nonpos_integer(b)) then
            value = hypergeo_aorb_nonpos_int(a, b, c, z, tol)
            return
        end if
        if (is_nonpos_integer(c)) then
            value = nan_complex()
            return
        end if
        if (real(a, dp) > real(b, dp)) then
            swap = a; a = b; b = swap
        end if
        m = c - a
        n = b - a
        if (is_near_integer(m)) then
            if (real(m, dp) <= 0.0_dp) then
                value = (1.0_dp - z) ** (c - a - b) &
                    * hypergeo_powerseries(c - a, c - b, c, z, tol, maxiter, d + 1)
                return
            else if (is_near_integer(n)) then
                value = hypergeo_cover3(a, n, m, z, tol, maxiter)
                return
            end if
        end if
        m = c - a - b
        if (is_near_integer(m)) then
            value = hypergeo_cover1(a, b, m, z, tol, maxiter)
            return
        end if
        m = b - a
        if (is_near_integer(m)) then
            value = hypergeo_cover2(a, c, m, z, tol, maxiter)
            return
        end if
        value = hypergeo_general(a, b, c, z, tol, maxiter)
    end function hypergeo_powerseries

    recursive function hypergeo_core(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        complex(dp) :: cp(2)
        cp = crit_points()
        if (abs(z - cp(1)) < 0.1_dp .or. abs(z - cp(2)) < 0.1_dp) then
            value = hypergeo_gosper(a, b, c, z, tol, maxiter)
        else
            value = hypergeo_powerseries(a, b, c, z, tol, maxiter)
        end if
    end function hypergeo_core

    function hypergeo_gosper(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, d, e, f, dnew, enew, fnew, kk
        integer :: k, nmax
        real(dp) :: t
        d = 0.0_dp
        e = 1.0_dp
        f = 0.0_dp
        t = effective_tol(tol)
        nmax = effective_maxiter(maxiter)
        do k = 0, nmax
            kk = cmplx(real(k, dp), 0.0_dp, dp)
            dnew = (kk + a) * (kk + b) * z &
                * (e - (kk + c - b - a) * d * z / (1.0_dp - z)) &
                / (4.0_dp * (kk + 1.0_dp) * (kk + c / 2.0_dp) * (kk + (c + 1.0_dp) / 2.0_dp))
            enew = (kk + a) * (kk + b) * z &
                * (a * b * d * z / (1.0_dp - z) + (kk + c) * e) &
                / (4.0_dp * (kk + 1.0_dp) * (kk + c / 2.0_dp) * (kk + (c + 1.0_dp) / 2.0_dp))
            fnew = f - d * (kk * ((c - b - a) * z + kk * (z - 2.0_dp) - c) - a * b * z) &
                / (2.0_dp * (kk + c / 2.0_dp) * (1.0_dp - z)) + e
            if (good_delta(f - fnew, fnew, t)) then
                value = f
                return
            end if
            d = dnew
            e = enew
            f = fnew
        end do
        value = f
    end function hypergeo_gosper

    pure real(dp) function effective_tol(tol) result(t)
        real(dp), intent(in), optional :: tol
        t = 0.0_dp
        if (present(tol)) t = max(0.0_dp, tol)
    end function effective_tol

    pure integer function effective_maxiter(maxiter) result(n)
        integer, intent(in), optional :: maxiter
        n = 2000
        if (present(maxiter)) n = max(0, maxiter)
    end function effective_maxiter

    pure logical function good_delta(delta, value, tol) result(ok)
        complex(dp), intent(in) :: delta, value
        real(dp), intent(in) :: tol
        ok = abs(delta) <= tol
        if (tol <= 0.0_dp) ok = abs(delta) <= epsilon(1.0_dp) * max(1.0_dp, abs(value))
    end function good_delta

    pure real(dp) function parity_sign(n) result(s)
        integer, intent(in) :: n
        if (mod(abs(n), 2) == 0) then
            s = 1.0_dp
        else
            s = -1.0_dp
        end if
    end function parity_sign


    function f15_3_10_a(a,b,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,b,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_10(a,b,z,tol,maxiter,'a')
    end function f15_3_10_a
    function f15_3_10_b(a,b,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,b,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_10(a,b,z,tol,maxiter,'b')
    end function f15_3_10_b
    function f15_3_11_bit1(a,b,m,z,tol) result(v)
        complex(dp),intent(in)::a,b,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=f15_3_11_bit1_impl(a,b,nint(real(m,dp)),z,tol)
    end function f15_3_11_bit1
    function f15_3_11_bit2_a(a,b,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,b,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_11_bit2_impl(a,b,nint(real(m,dp)),z,tol,maxiter)
    end function f15_3_11_bit2_a
    function f15_3_11_bit2_b(a,b,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,b,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_11_bit2_impl(a,b,nint(real(m,dp)),z,tol,maxiter)
    end function f15_3_11_bit2_b
    function f15_3_12_bit1(a,b,m,z,tol) result(v)
        complex(dp),intent(in)::a,b,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=f15_3_12_bit1_impl(a,b,nint(real(m,dp)),z,tol)
    end function f15_3_12_bit1
    function f15_3_12_bit2_a(a,b,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,b,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_12_bit2_impl(a,b,nint(real(m,dp)),z,tol,maxiter)
    end function f15_3_12_bit2_a
    function f15_3_12_bit2_b(a,b,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,b,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_12_bit2_impl(a,b,nint(real(m,dp)),z,tol,maxiter)
    end function f15_3_12_bit2_b
    function f15_3_13_a(a,c,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,c,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_13(a,c,z,tol,maxiter,'a')
    end function f15_3_13_a
    function f15_3_13_b(a,c,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,c,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_13(a,c,z,tol,maxiter,'b')
    end function f15_3_13_b
    function f15_3_14_bit1_a(a,c,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,c,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_14_bit1_impl(a,c,nint(real(m,dp)),z,tol,maxiter)
    end function f15_3_14_bit1_a
    function f15_3_14_bit1_b(a,c,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,c,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=f15_3_14_bit1_impl(a,c,nint(real(m,dp)),z,tol,maxiter)
    end function f15_3_14_bit1_b
    function f15_3_14_bit2(a,c,m,z,tol) result(v)
        complex(dp),intent(in)::a,c,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=f15_3_14_bit2_impl(a,c,nint(real(m,dp)),z,tol)
    end function f15_3_14_bit2
    function w07_23_06_0031_01_bit1(a,n,m,z,tol) result(v)
        complex(dp),intent(in)::a,n,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=w0031_bit1(a,nint(real(n,dp)),nint(real(m,dp)),z,tol)
    end function w07_23_06_0031_01_bit1
    function w07_23_06_0031_01_bit2(a,n,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,n,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=w0031_bit2(a,nint(real(n,dp)),nint(real(m,dp)),z,tol,maxiter)
    end function w07_23_06_0031_01_bit2
    function w07_23_06_0026_01_bit1(a,n,m,z,tol) result(v)
        complex(dp),intent(in)::a,n,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=w0026_bit1(a,nint(real(n,dp)),nint(real(m,dp)),z,tol)
    end function w07_23_06_0026_01_bit1
    function w07_23_06_0026_01_bit2(a,n,m,z,tol,maxiter) result(v)
        complex(dp),intent(in)::a,n,m,z; real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter; complex(dp)::v
        v=w0026_bit2(a,nint(real(n,dp)),nint(real(m,dp)),z,tol,maxiter)
    end function w07_23_06_0026_01_bit2
    function w07_23_06_0026_01_bit3_a(a,n,m,z,tol) result(v)
        complex(dp),intent(in)::a,n,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=w0026_bit3(a,nint(real(n,dp)),nint(real(m,dp)),z); if(present(tol)) continue
    end function w07_23_06_0026_01_bit3_a
    function w07_23_06_0026_01_bit3_b(a,n,m,z,tol) result(v)
        complex(dp),intent(in)::a,n,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=w0026_bit3(a,nint(real(n,dp)),nint(real(m,dp)),z); if(present(tol)) continue
    end function w07_23_06_0026_01_bit3_b
    function w07_23_06_0026_01_bit3_c(a,n,m,z,tol) result(v)
        complex(dp),intent(in)::a,n,m,z; real(dp),intent(in),optional::tol; complex(dp)::v
        v=w0026_bit3(a,nint(real(n,dp)),nint(real(m,dp)),z); if(present(tol)) continue
    end function w07_23_06_0026_01_bit3_c

end module hypergeo_gauss_core
