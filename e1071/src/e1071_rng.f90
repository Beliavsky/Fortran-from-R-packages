module e1071_rng
    use e1071_kinds, only: dp, i8
    use e1071_constants, only: e1071_pi
    implicit none
    private

    type, public :: rng_state
        integer(i8) :: state = 104729_i8
    end type rng_state

    public :: rng_seed, rng_uniform, rng_normal, rng_integer

contains

    subroutine rng_seed(rng, seed)
        type(rng_state), intent(inout) :: rng !! Mutable pseudo-random-number generator state to initialize.
        integer, intent(in) :: seed !! Integer seed; zero and negative values are mapped into the valid Park-Miller state range.
        integer(i8), parameter :: modulus = 2147483647_i8
        integer(i8) :: s

        s = modulo(int(seed, i8), modulus - 1_i8)
        if (s < 0_i8) s = s + modulus - 1_i8
        rng%state = s + 1_i8
    end subroutine rng_seed

    function rng_uniform(rng) result(u)
        type(rng_state), intent(inout) :: rng !! Mutable generator state advanced by one uniform draw.
        real(dp) :: u
        integer(i8), parameter :: multiplier = 16807_i8
        integer(i8), parameter :: modulus = 2147483647_i8

        rng%state = modulo(multiplier * rng%state, modulus)
        if (rng%state <= 0_i8) rng%state = 1_i8
        u = real(rng%state, dp) / real(modulus, dp)
    end function rng_uniform

    function rng_normal(rng) result(z)
        type(rng_state), intent(inout) :: rng !! Mutable generator state advanced by two uniforms for one standard-normal draw.
        real(dp) :: z
        real(dp) :: u1
        real(dp) :: u2

        u1 = max(rng_uniform(rng), tiny(1.0_dp))
        u2 = rng_uniform(rng)
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * e1071_pi * u2)
    end function rng_normal

    function rng_integer(rng, n) result(k)
        type(rng_state), intent(inout) :: rng !! Mutable generator state advanced by one draw.
        integer, intent(in) :: n !! Number of equiprobable integer outcomes; must be positive.
        integer :: k

        if (n <= 0) error stop "rng_integer: n must be positive"
        k = min(n, 1 + int(rng_uniform(rng) * real(n, dp)))
    end function rng_integer

end module e1071_rng
