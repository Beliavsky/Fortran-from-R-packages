! Portable RNG support for the Fortran translation.
! DEoptim itself obtains random variates from R's RNG; this module supplies a
! deterministic standalone replacement because an FPM library cannot depend on R.
module deoptim_rng
    use deoptim_kinds, only : dp, i8
    implicit none
    private

    integer(i8), parameter :: pm_m = 2147483647_i8
    integer(i8), parameter :: pm_a = 16807_i8
    integer(i8), parameter :: pm_q = 127773_i8
    integer(i8), parameter :: pm_r = 2836_i8
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: de_rng
        integer(i8) :: state = 1_i8
        logical :: have_normal = .false.
        real(dp) :: normal_cache = 0.0_dp
    contains
        procedure :: seed => rng_seed
        procedure :: uniform => rng_uniform
        procedure :: normal => rng_normal
        procedure :: cauchy => rng_cauchy
        procedure :: randint => rng_randint
    end type de_rng

contains

    subroutine rng_seed(this, seed_value)
        class(de_rng), intent(inout) :: this
        integer(i8), intent(in) :: seed_value
        integer(i8) :: s
        integer :: count

        if (seed_value == 0_i8) then
            call system_clock(count=count)
            s = int(count, i8)
        else
            s = seed_value
        end if
        s = modulo(abs(s), pm_m - 1_i8) + 1_i8
        this%state = s
        this%have_normal = .false.
        this%normal_cache = 0.0_dp
    end subroutine rng_seed

    function rng_uniform(this) result(u)
        class(de_rng), intent(inout) :: this
        real(dp) :: u
        integer(i8) :: hi, lo, test

        hi = this%state / pm_q
        lo = modulo(this%state, pm_q)
        test = pm_a * lo - pm_r * hi
        if (test > 0_i8) then
            this%state = test
        else
            this%state = test + pm_m
        end if
        u = real(this%state, dp) / real(pm_m, dp)
    end function rng_uniform

    function rng_randint(this, n) result(k)
        class(de_rng), intent(inout) :: this
        integer, intent(in) :: n
        integer :: k

        if (n <= 1) then
            k = 1
        else
            k = 1 + int(this%uniform() * real(n, dp))
            if (k > n) k = n
        end if
    end function rng_randint

    function rng_normal(this, mean, sd) result(x)
        class(de_rng), intent(inout) :: this
        real(dp), intent(in) :: mean, sd
        real(dp) :: x
        real(dp) :: u1, u2, radius, z

        if (this%have_normal) then
            z = this%normal_cache
            this%have_normal = .false.
        else
            u1 = max(this%uniform(), tiny(1.0_dp))
            u2 = this%uniform()
            radius = sqrt(-2.0_dp * log(u1))
            z = radius * cos(2.0_dp * pi * u2)
            this%normal_cache = radius * sin(2.0_dp * pi * u2)
            this%have_normal = .true.
        end if
        x = mean + sd * z
    end function rng_normal

    function rng_cauchy(this, location, scale) result(x)
        class(de_rng), intent(inout) :: this
        real(dp), intent(in) :: location, scale
        real(dp) :: x
        real(dp) :: u

        u = this%uniform()
        x = location + scale * tan(pi * (u - 0.5_dp))
    end function rng_cauchy

end module deoptim_rng
