! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_special
    use goftest_kinds, only : dp, pi
    implicit none
    private

    public :: bessel_k_frac

contains

    function bessel_k_frac(nu, z) result(kv)
        real(dp), intent(in) :: nu, z
        real(dp) :: kv
        real(dp) :: tmax
        integer :: nseg, j

        if (z <= 0.0_dp) then
            kv = huge(1.0_dp)
            return
        end if

        if (z > 80.0_dp) then
            kv = bessel_k_asymptotic(nu, z)
            return
        end if

        tmax = max(8.0_dp, log(2.0_dp / z) + 8.0_dp)
        tmax = min(tmax, 40.0_dp)
        nseg = max(8, ceiling(tmax))
        kv = 0.0_dp
        do j = 0, nseg - 1
            kv = kv + gauss16_segment(nu, z, &
                tmax * real(j, dp) / real(nseg, dp), &
                tmax * real(j + 1, dp) / real(nseg, dp))
        end do
    end function bessel_k_frac

    pure function bessel_k_asymptotic(nu, z) result(kv)
        real(dp), intent(in) :: nu, z
        real(dp) :: kv
        real(dp) :: mu, term, series
        integer :: k

        mu = 4.0_dp * nu * nu
        series = 1.0_dp
        term = 1.0_dp
        do k = 1, 8
            term = term * (mu - real((2 * k - 1)**2, dp)) / &
                (real(k, dp) * 8.0_dp * z)
            series = series + term
        end do
        kv = sqrt(pi / (2.0_dp * z)) * exp(-z) * series
    end function bessel_k_asymptotic

    function gauss16_segment(nu, z, a, b) result(ans)
        real(dp), intent(in) :: nu, z, a, b
        real(dp) :: ans
        real(dp), parameter :: x(8) = [ &
            0.095012509837637440185_dp, 0.281603550779258913230_dp, &
            0.458016777657227386342_dp, 0.617876244402643748447_dp, &
            0.755404408355003033895_dp, 0.865631202387831743880_dp, &
            0.944575023073232576078_dp, 0.989400934991649932596_dp ]
        real(dp), parameter :: w(8) = [ &
            0.189450610455068496285_dp, 0.182603415044923588867_dp, &
            0.169156519395002538189_dp, 0.149595988816576732081_dp, &
            0.124628971255533872052_dp, 0.095158511682492784810_dp, &
            0.062253523938647892863_dp, 0.027152459411754094852_dp ]
        real(dp) :: mid, half, t1, t2
        integer :: i

        mid = 0.5_dp * (a + b)
        half = 0.5_dp * (b - a)
        ans = 0.0_dp
        do i = 1, 8
            t1 = mid - half * x(i)
            t2 = mid + half * x(i)
            ans = ans + w(i) * (k_integrand(nu, z, t1) + k_integrand(nu, z, t2))
        end do
        ans = half * ans
    end function gauss16_segment

    pure function k_integrand(nu, z, t) result(f)
        real(dp), intent(in) :: nu, z, t
        real(dp) :: f, logf, c

        c = cosh(t)
        logf = -z * c + log(cosh(nu * t))
        if (logf < log(tiny(1.0_dp))) then
            f = 0.0_dp
        else
            f = exp(logf)
        end if
    end function k_integrand

end module goftest_special
