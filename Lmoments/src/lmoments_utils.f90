! SPDX-License-Identifier: GPL-2.0-only
module lmoments_utils
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter, public :: pi = acos(-1.0_dp)

    public :: sort_real
    public :: choose_real
    public :: standard_normal_quantile

contains

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j, m, n
        real(dp) :: test

        n = size(x)
        if (n <= 1) return

        m = 1
        do while (3 * m + 1 < n)
            m = 3 * m + 1
        end do

        do while (m > 0)
            do i = m + 1, n
                test = x(i)
                j = i
                do while (j > m)
                    if (test >= x(j - m)) exit
                    x(j) = x(j - m)
                    j = j - m
                end do
                x(j) = test
            end do
            m = m / 3
        end do
    end subroutine sort_real


    pure real(dp) function choose_real(n, k) result(value)
        integer, intent(in) :: n, k
        integer :: i, kk

        if (n < 0 .or. k < 0 .or. k > n) then
            value = 0.0_dp
            return
        end if
        kk = min(k, n - k)
        value = 1.0_dp
        do i = 1, kk
            value = value * real(n - kk + i, dp) / real(i, dp)
        end do
    end function choose_real


    ! AS241, Wichura (1988), as used by J. R. M. Hosking's LMOMENTS 3.04.
    ! Hosking's original source is distributed under the IBM permission notice
    ! included in licenses/HOSKING-LMOMENTS.txt.
    pure real(dp) function standard_normal_quantile(p) result(qn)
        real(dp), intent(in) :: p
        real(dp), parameter :: split1 = 0.425_dp
        real(dp), parameter :: split2 = 5.0_dp
        real(dp), parameter :: const1 = 0.180625_dp
        real(dp), parameter :: const2 = 1.6_dp
        real(dp), parameter :: a(0:7) = [ &
            3.3871328727963666080_dp, 1.3314166789178437745e2_dp, &
            1.9715909503065514427e3_dp, 1.3731693765509461125e4_dp, &
            4.5921953931549871457e4_dp, 6.7265770927008700853e4_dp, &
            3.3430575583588128105e4_dp, 2.5090809287301226727e3_dp ]
        real(dp), parameter :: b(1:7) = [ &
            4.2313330701600911252e1_dp, 6.8718700749205790830e2_dp, &
            5.3941960214247511077e3_dp, 2.1213794301586595867e4_dp, &
            3.9307895800092710610e4_dp, 2.8729085735721942674e4_dp, &
            5.2264952788528545610e3_dp ]
        real(dp), parameter :: c(0:7) = [ &
            1.42343711074968357734_dp, 4.63033784615654529590_dp, &
            5.76949722146069140550_dp, 3.64784832476320460504_dp, &
            1.27045825245236838258_dp, 2.41780725177450611770e-1_dp, &
            2.27238449892691845833e-2_dp, 7.74545014278341407640e-4_dp ]
        real(dp), parameter :: d(1:7) = [ &
            2.05319162663775882187_dp, 1.67638483018380384940_dp, &
            6.89767334985100004550e-1_dp, 1.48103976427480074590e-1_dp, &
            1.51986665636164571966e-2_dp, 5.47593808499534494600e-4_dp, &
            1.05075007164441684324e-9_dp ]
        real(dp), parameter :: e(0:7) = [ &
            6.65790464350110377720_dp, 5.46378491116411436990_dp, &
            1.78482653991729133580_dp, 2.96560571828504891230e-1_dp, &
            2.65321895265761230930e-2_dp, 1.24266094738807843860e-3_dp, &
            2.71155556874348757815e-5_dp, 2.01033439929228813265e-7_dp ]
        real(dp), parameter :: f(1:7) = [ &
            5.99832206555887937690e-1_dp, 1.36929880922735805310e-1_dp, &
            1.48753612908506148525e-2_dp, 7.86869131145613259100e-4_dp, &
            1.84631831751005468180e-5_dp, 1.42151175831644588870e-7_dp, &
            2.04426310338993978564e-15_dp ]
        real(dp) :: q, r, num, den
        integer :: i

        if (p <= 0.0_dp) then
            qn = ieee_value(0.0_dp, ieee_negative_inf)
            return
        end if
        if (p >= 1.0_dp) then
            qn = ieee_value(0.0_dp, ieee_positive_inf)
            return
        end if

        q = p - 0.5_dp
        if (abs(q) <= split1) then
            r = const1 - q * q
            num = a(7)
            do i = 6, 0, -1
                num = num * r + a(i)
            end do
            den = b(7)
            do i = 6, 1, -1
                den = den * r + b(i)
            end do
            den = den * r + 1.0_dp
            qn = q * num / den
            return
        end if

        if (q < 0.0_dp) then
            r = p
        else
            r = 1.0_dp - p
        end if
        r = sqrt(-log(r))

        if (r <= split2) then
            r = r - const2
            num = c(7)
            do i = 6, 0, -1
                num = num * r + c(i)
            end do
            den = d(7)
            do i = 6, 1, -1
                den = den * r + d(i)
            end do
            den = den * r + 1.0_dp
        else
            r = r - split2
            num = e(7)
            do i = 6, 0, -1
                num = num * r + e(i)
            end do
            den = f(7)
            do i = 6, 1, -1
                den = den * r + f(i)
            end do
            den = den * r + 1.0_dp
        end if
        qn = num / den
        if (q < 0.0_dp) qn = -qn
    end function standard_normal_quantile

end module lmoments_utils
