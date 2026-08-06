! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_random
    use iso_fortran_env, only : int64
    use ltsa_kinds, only : dp, two_pi
    implicit none
    private

    integer(int64), save :: state = 88172645463325252_int64
    logical, save :: have_spare = .false.
    real(dp), save :: spare = 0.0_dp

    public :: set_ltsa_seed, ltsa_uniform, ltsa_normal

contains

    subroutine set_ltsa_seed(seed)
        integer(int64), intent(in) :: seed
        if (seed == 0_int64) then
            state = 88172645463325252_int64
        else
            state = seed
        end if
        have_spare = .false.
    end subroutine set_ltsa_seed

    real(dp) function ltsa_uniform() result(u)
        integer(int64) :: x
        x = state
        x = ieor(x, shiftl(x, 13))
        x = ieor(x, shiftr(x, 7))
        x = ieor(x, shiftl(x, 17))
        state = x
        u = real(iand(x, int(z'001FFFFFFFFFFFFF', int64)), dp) / real(int(z'0020000000000000', int64), dp)
        if (u <= 0.0_dp) u = epsilon(1.0_dp)
    end function ltsa_uniform

    real(dp) function ltsa_normal() result(z)
        real(dp) :: u1, u2, radius
        if (have_spare) then
            z = spare
            have_spare = .false.
            return
        end if
        u1 = ltsa_uniform()
        u2 = ltsa_uniform()
        radius = sqrt(-2.0_dp*log(u1))
        z = radius*cos(two_pi*u2)
        spare = radius*sin(two_pi*u2)
        have_spare = .true.
    end function ltsa_normal

end module ltsa_random
