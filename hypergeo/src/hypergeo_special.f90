! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_special
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use hypergeo_kinds, only : dp, pi, ci
    implicit none
    private

    real(dp), parameter :: near_int_tol = 1.0e-11_dp

    public :: complex_gamma, complex_log_gamma, complex_factorial
    public :: complex_digamma, gamma_ratio3, gamma_ratio4
    public :: is_near_integer, is_nonpos_integer, is_zero_parameter
    public :: factorial_int, pochhammer_int, log_pochhammer
    public :: nan_complex, finite_complex

contains

    pure function nan_complex() result(z)
        complex(dp) :: z
        real(dp) :: qnan
        qnan = ieee_value(0.0_dp, ieee_quiet_nan)
        z = cmplx(qnan, qnan, kind=dp)
    end function nan_complex

    pure logical function finite_complex(z) result(ok)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
        complex(dp), intent(in) :: z
        ok = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex

    pure logical function is_near_integer(z, tol) result(ans)
        complex(dp), intent(in) :: z
        real(dp), intent(in), optional :: tol
        real(dp) :: t
        t = near_int_tol
        if (present(tol)) t = max(tol, 0.0_dp)
        ans = abs(aimag(z)) <= t .and. abs(real(z, dp) - anint(real(z, dp))) <= t
    end function is_near_integer

    pure logical function is_nonpos_integer(z, tol) result(ans)
        complex(dp), intent(in) :: z
        real(dp), intent(in), optional :: tol
        if (present(tol)) then
            ans = is_near_integer(z, tol) .and. real(z, dp) < 0.5_dp
        else
            ans = is_near_integer(z) .and. real(z, dp) < 0.5_dp
        end if
    end function is_nonpos_integer

    pure logical function is_zero_parameter(z, tol) result(ans)
        complex(dp), intent(in) :: z
        real(dp), intent(in), optional :: tol
        if (present(tol)) then
            ans = is_near_integer(z, tol) .and. abs(z) < 0.5_dp
        else
            ans = is_near_integer(z) .and. abs(z) < 0.5_dp
        end if
    end function is_zero_parameter

    pure function lanczos_log(z) result(value)
        complex(dp), intent(in) :: z
        complex(dp) :: value
        complex(dp) :: zz, x, tee
        real(dp), parameter :: g = 7.0_dp
        real(dp), parameter :: p(9) = [ &
            0.99999999999980993227684700473478_dp, &
            676.520368121885098567009190444019_dp, &
            -1259.13921672240287047156078755283_dp, &
            771.3234287776530788486528258894_dp, &
            -176.61502916214059906584551354_dp, &
            12.507343278686904814458936853_dp, &
            -0.13857109526572011689554707_dp, &
            9.984369578019570859563e-6_dp, &
            1.50563273514931155834e-7_dp ]
        integer :: i

        zz = z - 1.0_dp
        x = cmplx(p(1), 0.0_dp, kind=dp)
        do i = 2, 9
            x = x + p(i) / (zz + real(i - 1, dp))
        end do
        tee = zz + g + 0.5_dp
        value = 0.5_dp * log(2.0_dp * pi) + (zz + 0.5_dp) * log(tee) - tee + log(x)
    end function lanczos_log

    recursive pure function complex_log_gamma(z) result(value)
        complex(dp), intent(in) :: z
        complex(dp) :: value
        if (real(z, dp) < 0.5_dp) then
            value = log(pi) - log(sin(pi * z)) - lanczos_log(1.0_dp - z)
        else
            value = lanczos_log(z)
        end if
    end function complex_log_gamma

    pure function complex_gamma(z) result(value)
        complex(dp), intent(in) :: z
        complex(dp) :: value
        value = exp(complex_log_gamma(z))
    end function complex_gamma

    pure function complex_factorial(z) result(value)
        complex(dp), intent(in) :: z
        complex(dp) :: value
        value = complex_gamma(z + 1.0_dp)
    end function complex_factorial

    recursive pure function complex_digamma(z) result(value)
        complex(dp), intent(in) :: z
        complex(dp) :: value, w, inv, inv2
        integer :: k

        if (real(z, dp) < 0.5_dp) then
            value = complex_digamma(1.0_dp - z) - pi / tan(pi * z)
            return
        end if

        w = z
        value = (0.0_dp, 0.0_dp)
        do k = 1, 100
            if (abs(w) >= 10.0_dp) exit
            value = value - 1.0_dp / w
            w = w + 1.0_dp
        end do
        inv = 1.0_dp / w
        inv2 = inv * inv
        value = value + log(w) - 0.5_dp * inv &
            - inv2 * (1.0_dp / 12.0_dp &
            - inv2 * (1.0_dp / 120.0_dp &
            - inv2 * (1.0_dp / 252.0_dp &
            - inv2 * (1.0_dp / 240.0_dp &
            - inv2 * (5.0_dp / 660.0_dp)))))
    end function complex_digamma

    pure function gamma_ratio3(a1, a2, a3) result(value)
        complex(dp), intent(in) :: a1, a2, a3
        complex(dp) :: value
        value = exp(complex_log_gamma(a1) - complex_log_gamma(a2) - complex_log_gamma(a3))
    end function gamma_ratio3

    pure function gamma_ratio4(a1, a2, a3, a4) result(value)
        complex(dp), intent(in) :: a1, a2, a3, a4
        complex(dp) :: value
        value = exp(complex_log_gamma(a1) + complex_log_gamma(a2) &
            - complex_log_gamma(a3) - complex_log_gamma(a4))
    end function gamma_ratio4

    pure real(dp) function factorial_int(n) result(value)
        integer, intent(in) :: n
        integer :: k
        if (n < 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        value = 1.0_dp
        do k = 2, n
            value = value * real(k, dp)
        end do
    end function factorial_int

    pure function pochhammer_int(z, n) result(value)
        complex(dp), intent(in) :: z
        integer, intent(in) :: n
        complex(dp) :: value
        integer :: k
        if (n < 0) then
            value = nan_complex()
            return
        end if
        value = (1.0_dp, 0.0_dp)
        do k = 0, n - 1
            value = value * (z + real(k, dp))
        end do
    end function pochhammer_int

    pure function log_pochhammer(z, n) result(value)
        complex(dp), intent(in) :: z
        integer, intent(in) :: n
        complex(dp) :: value
        if (n < 0) then
            value = nan_complex()
        else if (n == 0) then
            value = (0.0_dp, 0.0_dp)
        else
            value = complex_log_gamma(z + real(n, dp)) - complex_log_gamma(z)
        end if
    end function log_pochhammer

end module hypergeo_special
