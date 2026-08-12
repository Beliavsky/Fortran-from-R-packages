! Portable deterministic RNG used instead of R's global RNG.
module genalg_rng
    use iso_fortran_env, only : int64
    use genalg_kinds, only : dp
    implicit none
    private

    integer(int64), parameter :: pm_mod = 2147483647_int64
    integer(int64), parameter :: pm_mul = 16807_int64

    type, public :: rng_state
        integer(int64) :: state = 123456789_int64
    contains
        procedure :: seed => rng_seed
        procedure :: uniform => rng_uniform
        procedure :: integer => rng_integer
    end type rng_state

contains

    subroutine rng_seed(self, seed)
        class(rng_state), intent(inout) :: self
        integer(int64), intent(in) :: seed

        self%state = modulo(abs(seed), pm_mod - 1_int64) + 1_int64
    end subroutine rng_seed

    function rng_uniform(self) result(u)
        class(rng_state), intent(inout) :: self
        real(dp) :: u

        self%state = modulo(pm_mul * self%state, pm_mod)
        if (self%state == 0_int64) self%state = 1_int64
        u = real(self%state, dp) / real(pm_mod, dp)
    end function rng_uniform

    function rng_integer(self, lo, hi) result(k)
        class(rng_state), intent(inout) :: self
        integer, intent(in) :: lo, hi
        integer :: k
        real(dp) :: u

        if (hi < lo) error stop "rng_integer: invalid bounds"
        u = self%uniform()
        k = lo + int(u * real(hi - lo + 1, dp))
        if (k > hi) k = hi
    end function rng_integer

end module genalg_rng
