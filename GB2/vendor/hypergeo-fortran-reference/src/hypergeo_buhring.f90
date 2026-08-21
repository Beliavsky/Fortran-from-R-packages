! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_buhring_mod
    use hypergeo_kinds, only : dp
    use hypergeo_special, only : complex_log_gamma, factorial_int, finite_complex, nan_complex
    use hypergeo_gauss_core, only : hypergeo_core
    implicit none
    private
    public :: lpham, buhring_eqn11, buhring_eqn12, buhring_eqn5_factors
    public :: buhring_eqn5_series, hypergeo_buhring
contains
    pure function lpham(x, n) result(v)
        complex(dp), intent(in) :: x
        integer, intent(in) :: n
        complex(dp) :: v
        v = complex_log_gamma(x + real(n, dp)) - complex_log_gamma(x)
    end function lpham

    function buhring_eqn11(n, s, a, b, c, z0) result(v)
        integer, intent(in) :: n
        complex(dp), intent(in) :: s, a, b, c
        complex(dp), intent(in), optional :: z0
        complex(dp) :: v, zz
        zz = cmplx(0.5_dp, 0.0_dp, dp); if (present(z0)) zz = z0
        v = exp(lpham(s, n) + lpham(1.0_dp + s - c, n) &
            - lpham(1.0_dp + 2.0_dp * s - a - b, n) - log(factorial_int(n))) &
            * hypergeo_core(cmplx(-real(n, dp), 0.0_dp, dp), a + b - 2.0_dp * s - real(n, dp), &
                c - s - real(n, dp), zz)
    end function buhring_eqn11

    function buhring_eqn12(n, s, a, b, c, z0) result(v)
        integer, intent(in) :: n
        complex(dp), intent(in) :: s, a, b, c
        complex(dp), intent(in), optional :: z0
        complex(dp) :: v, zz
        real(dp) :: sg
        zz = cmplx(0.5_dp, 0.0_dp, dp); if (present(z0)) zz = z0
        sg = merge(1.0_dp, -1.0_dp, mod(n, 2) == 0)
        v = sg * exp(lpham(s, n) + lpham(s + c - a - b, n) &
            - lpham(1.0_dp + 2.0_dp * s - a - b, n) - log(factorial_int(n))) &
            * hypergeo_core(cmplx(-real(n, dp), 0.0_dp, dp), a + b - 2.0_dp * s - real(n, dp), &
                1.0_dp + a + b - s - c - real(n, dp), 1.0_dp - zz)
    end function buhring_eqn12

    pure function buhring_eqn5_factors(a, b, c, z, z0) result(v)
        complex(dp), intent(in) :: a, b, c, z
        complex(dp), intent(in), optional :: z0
        complex(dp) :: v(2), zz
        zz = cmplx(0.5_dp, 0.0_dp, dp); if (present(z0)) zz = z0
        v(1) = exp(complex_log_gamma(c) + complex_log_gamma(b - a) - complex_log_gamma(b) &
            - complex_log_gamma(c - a) - a * log(zz - z))
        v(2) = exp(complex_log_gamma(c) + complex_log_gamma(a - b) - complex_log_gamma(a) &
            - complex_log_gamma(c - b) - b * log(zz - z))
    end function buhring_eqn5_factors

    function buhring_eqn5_series(s, a, b, c, z, z0, use11, tol, maxiter) result(v)
        complex(dp), intent(in) :: s, a, b, c, z
        complex(dp), intent(in), optional :: z0
        logical, intent(in), optional :: use11
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: v, old, zz, term
        logical :: u11
        real(dp) :: t
        integer :: n, nmax
        zz = cmplx(0.5_dp, 0.0_dp, dp); if (present(z0)) zz = z0
        u11 = .false.; if (present(use11)) u11 = use11
        t = 0.0_dp; if (present(tol)) t = max(0.0_dp, tol)
        nmax = 2000; if (present(maxiter)) nmax = max(1, maxiter)
        v = 1.0_dp
        do n = 1, nmax - 1
            old = v
            if (u11) then
                term = buhring_eqn11(n, s, a, b, c, zz)
            else
                term = buhring_eqn12(n, s, a, b, c, zz)
            end if
            v = old + term / (z - zz) ** n
            if (abs(v - old) <= max(t, epsilon(1.0_dp) * max(1.0_dp, abs(v)))) return
        end do
        if (.not. finite_complex(v)) v = nan_complex()
    end function buhring_eqn5_series

    function hypergeo_buhring(a, b, c, z, z0, tol, maxiter, use11) result(v)
        complex(dp), intent(in) :: a, b, c, z
        complex(dp), intent(in), optional :: z0
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        logical, intent(in), optional :: use11
        complex(dp) :: v, factors(2), zz
        zz = cmplx(0.5_dp, 0.0_dp, dp); if (present(z0)) zz = z0
        factors = buhring_eqn5_factors(a, b, c, z, zz)
        v = factors(1) * buhring_eqn5_series(a, a, b, c, z, zz, use11, tol, maxiter) &
            + factors(2) * buhring_eqn5_series(b, a, b, c, z, zz, use11, tol, maxiter)
    end function hypergeo_buhring
end module hypergeo_buhring_mod
