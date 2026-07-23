! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

module longmemo_stats
    use longmemo_kinds, only : dp
    implicit none
    private

    public :: normal_cdf, normal_sf, student_t_cdf, sample_sd
    public :: random_normal_vector, random_exponential_vector, set_random_seed

contains

    elemental real(dp) function normal_cdf(x)
        real(dp), intent(in) :: x

        normal_cdf = 0.5_dp*erfc(-x/sqrt(2.0_dp))
    end function normal_cdf


    elemental real(dp) function normal_sf(x)
        real(dp), intent(in) :: x

        normal_sf = 0.5_dp*erfc(x/sqrt(2.0_dp))
    end function normal_sf


    pure real(dp) function student_t_cdf(t, nu)
        real(dp), intent(in) :: t, nu
        real(dp) :: x, ib

        if (nu <= 0.0_dp) then
            student_t_cdf = ieee_value_or_zero()
            return
        end if
        if (abs(t) <= tiny(1.0_dp)) then
            student_t_cdf = 0.5_dp
            return
        end if

        x = nu/(nu + t*t)
        ib = regularized_beta(x, 0.5_dp*nu, 0.5_dp)
        if (t > 0.0_dp) then
            student_t_cdf = 1.0_dp - 0.5_dp*ib
        else
            student_t_cdf = 0.5_dp*ib
        end if
    end function student_t_cdf


    pure real(dp) function sample_sd(x)
        real(dp), intent(in) :: x(:)
        real(dp) :: mean_x

        if (size(x) < 2) then
            sample_sd = 0.0_dp
            return
        end if
        mean_x = sum(x)/real(size(x), dp)
        sample_sd = sqrt(sum((x - mean_x)**2)/real(size(x) - 1, dp))
    end function sample_sd


    subroutine set_random_seed(seed)
        integer, intent(in) :: seed
        integer, allocatable :: put(:)
        integer :: n, i

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729*i, huge(1) - 1)
            if (put(i) == 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine set_random_seed


    subroutine random_normal_vector(x)
        real(dp), intent(out) :: x(:)
        real(dp) :: u1, u2, r, angle
        integer :: i

        i = 1
        do while (i <= size(x))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            r = sqrt(-2.0_dp*log(u1))
            angle = 2.0_dp*acos(-1.0_dp)*u2
            x(i) = r*cos(angle)
            if (i + 1 <= size(x)) x(i + 1) = r*sin(angle)
            i = i + 2
        end do
    end subroutine random_normal_vector


    subroutine random_exponential_vector(x)
        real(dp), intent(out) :: x(:)
        real(dp), allocatable :: u(:)

        allocate(u(size(x)))
        call random_number(u)
        x = -log(max(u, tiny(1.0_dp)))
    end subroutine random_exponential_vector


    pure real(dp) function regularized_beta(x, a, b)
        real(dp), intent(in) :: x, a, b
        real(dp) :: bt

        if (x <= 0.0_dp) then
            regularized_beta = 0.0_dp
            return
        else if (x >= 1.0_dp) then
            regularized_beta = 1.0_dp
            return
        end if

        bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
                 a*log(x) + b*log(1.0_dp - x))
        if (x < (a + 1.0_dp)/(a + b + 2.0_dp)) then
            regularized_beta = bt*beta_continued_fraction(a, b, x)/a
        else
            regularized_beta = 1.0_dp - bt*beta_continued_fraction(b, a, 1.0_dp - x)/b
        end if
    end function regularized_beta


    pure real(dp) function beta_continued_fraction(a, b, x)
        real(dp), intent(in) :: a, b, x
        integer, parameter :: max_iter = 300
        real(dp), parameter :: eps = 3.0e-14_dp
        real(dp), parameter :: fpmin = 1.0e-300_dp
        real(dp) :: aa, c, d, del, h, qab, qam, qap
        integer :: m, m2

        qab = a + b
        qap = a + 1.0_dp
        qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab*x/qap
        if (abs(d) < fpmin) d = fpmin
        d = 1.0_dp/d
        h = d

        do m = 1, max_iter
            m2 = 2*m
            aa = real(m, dp)*(b - real(m, dp))*x/ &
                 ((qam + real(m2, dp))*(a + real(m2, dp)))
            d = 1.0_dp + aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            h = h*d*c

            aa = -(a + real(m, dp))*(qab + real(m, dp))*x/ &
                 ((a + real(m2, dp))*(qap + real(m2, dp)))
            d = 1.0_dp + aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= eps) exit
        end do
        beta_continued_fraction = h
    end function beta_continued_fraction


    pure real(dp) function ieee_value_or_zero()
        ieee_value_or_zero = 0.0_dp/1.0_dp
    end function ieee_value_or_zero

end module longmemo_stats
