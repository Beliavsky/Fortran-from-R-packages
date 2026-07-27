! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_rng
    use iso_fortran_env, only: int64
    use derivmkts_kinds, only: dp, pi
    implicit none
    private
    public :: rng_state, seed_rng, uniform_rng, normal_rng, poisson_rng, binomial_rng

    type :: rng_state
        integer(int64) :: state = 88172645463393265_int64
        logical :: has_spare = .false.
        real(dp) :: spare = 0.0_dp
    end type rng_state

contains

    subroutine seed_rng(rng, seed)
        type(rng_state), intent(inout) :: rng
        integer, intent(in) :: seed
        rng%state = int(max(seed,1), int64)*6364136223846793005_int64 + 1442695040888963407_int64
        if (rng%state == 0_int64) rng%state = 88172645463393265_int64
        rng%has_spare = .false.
    end subroutine seed_rng

    real(dp) function uniform_rng(rng) result(u)
        type(rng_state), intent(inout) :: rng
        integer(int64) :: x
        x = rng%state
        x = ieor(x, shiftl(x, 13))
        x = ieor(x, shiftr(x, 7))
        x = ieor(x, shiftl(x, 17))
        rng%state = x
        u = real(iand(x, int(z'001FFFFFFFFFFFFF',int64)), dp)/real(int(z'0020000000000000',int64),dp)
        if (u <= 0.0_dp) u = epsilon(1.0_dp)
    end function uniform_rng

    real(dp) function normal_rng(rng) result(z)
        type(rng_state), intent(inout) :: rng
        real(dp) :: u1, u2, rad
        if (rng%has_spare) then
            z = rng%spare
            rng%has_spare = .false.
            return
        end if
        u1 = uniform_rng(rng)
        u2 = uniform_rng(rng)
        rad = sqrt(-2.0_dp*log(u1))
        z = rad*cos(2.0_dp*pi*u2)
        rng%spare = rad*sin(2.0_dp*pi*u2)
        rng%has_spare = .true.
    end function normal_rng

    integer function poisson_rng(rng, lambda) result(k)
        type(rng_state), intent(inout) :: rng
        real(dp), intent(in) :: lambda
        real(dp) :: l, p, z
        if (lambda <= 0.0_dp) then
            k = 0
        else if (lambda < 30.0_dp) then
            l = exp(-lambda)
            p = 1.0_dp
            k = 0
            do
                k = k + 1
                p = p*uniform_rng(rng)
                if (p <= l) exit
            end do
            k = k - 1
        else
            do
                z = lambda + sqrt(lambda)*normal_rng(rng)
                k = nint(z)
                if (k >= 0) exit
            end do
        end if
    end function poisson_rng

    integer function binomial_rng(rng, n, p) result(k)
        type(rng_state), intent(inout) :: rng
        integer, intent(in) :: n
        real(dp), intent(in) :: p
        integer :: i
        k = 0
        do i = 1, max(n,0)
            if (uniform_rng(rng) < p) k = k + 1
        end do
    end function binomial_rng

end module derivmkts_rng
