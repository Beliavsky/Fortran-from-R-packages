! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_random
    use lsmc_kinds, only : dp
    use lsmc_math, only : pi
    implicit none
    private

    public :: seed_random_number
    public :: fill_normal
    public :: fill_correlated_normals

contains

    subroutine seed_random_number(seed)
        integer, intent(in) :: seed
        integer, allocatable :: seed_values(:)
        integer :: i
        integer :: n
        integer(kind=8) :: state

        call random_seed(size=n)
        allocate(seed_values(n))
        state = abs(int(seed, kind=8)) + 104729_8
        do i = 1, n
            state = modulo(1664525_8 * state + 1013904223_8, 2147483647_8)
            seed_values(i) = int(max(1_8, state))
        end do
        call random_seed(put=seed_values)
    end subroutine seed_random_number

    subroutine fill_normal(x)
        real(dp), intent(out) :: x(:)
        real(dp) :: radius
        real(dp) :: theta
        real(dp) :: u1
        real(dp) :: u2
        integer :: i

        i = 1
        do while (i <= size(x))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            radius = sqrt(-2.0_dp * log(u1))
            theta = 2.0_dp * pi * u2
            x(i) = radius * cos(theta)
            if (i + 1 <= size(x)) x(i + 1) = radius * sin(theta)
            i = i + 2
        end do
    end subroutine fill_normal

    subroutine fill_correlated_normals(z1, z2, rho)
        real(dp), intent(out) :: z1(:)
        real(dp), intent(out) :: z2(:)
        real(dp), intent(in) :: rho
        real(dp), allocatable :: independent(:)

        if (size(z1) /= size(z2)) error stop "fill_correlated_normals: size mismatch"
        if (abs(rho) > 1.0_dp) error stop "fill_correlated_normals: abs(rho) must not exceed one"

        allocate(independent(size(z1)))
        call fill_normal(z1)
        call fill_normal(independent)
        z2 = rho * z1 + sqrt(max(0.0_dp, 1.0_dp - rho * rho)) * independent
    end subroutine fill_correlated_normals

end module lsmc_random
