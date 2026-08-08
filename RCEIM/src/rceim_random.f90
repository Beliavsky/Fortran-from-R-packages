! SPDX-License-Identifier: GPL-2.0-or-later
module rceim_random
    use rceim_kinds, only : dp
    implicit none
    private
    public :: rceim_set_seed, random_normal_scalar, fill_uniform, fill_normal
contains
    subroutine rceim_set_seed(seed_value)
        integer, intent(in) :: seed_value
        integer :: n, i
        integer, allocatable :: seed(:)
        integer(kind=8) :: z

        call random_seed(size=n)
        allocate(seed(n))
        z = int(seed_value, kind=8)
        if (z < 0_8) z = -z
        do i = 1, n
            z = modulo(1103515245_8*z + 12345_8 + 104729_8*int(i,8), 2147483647_8)
            seed(i) = int(max(1_8, z))
        end do
        call random_seed(put=seed)
    end subroutine rceim_set_seed

    real(dp) function random_normal_scalar() result(z)
        real(dp) :: u1, u2
        real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        z = sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
    end function random_normal_scalar

    subroutine fill_uniform(x, lo, hi)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in) :: lo, hi
        real(dp) :: u(size(x))
        call random_number(u)
        x = lo + (hi-lo)*u
    end subroutine fill_uniform

    subroutine fill_normal(x, mu, sigma)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in) :: mu, sigma
        integer :: i
        if (sigma <= 0.0_dp) then
            x = mu
            return
        end if
        do i = 1, size(x)
            x(i) = mu + sigma*random_normal_scalar()
        end do
    end subroutine fill_normal
end module rceim_random
