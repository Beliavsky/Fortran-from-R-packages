module ppso_rng
    use ppso_kinds, only : dp
    implicit none
    private

    type, public :: rng_state
        integer(kind=8) :: state = 88172645463393265_8
        logical :: have_spare = .false.
        real(dp) :: spare = 0.0_dp
    contains
        procedure :: seed => rng_seed
        procedure :: uniform => rng_uniform
        procedure :: normal => rng_normal
        procedure :: randint => rng_randint
    end type rng_state

    public :: latin_hypercube

contains

    subroutine rng_seed(self, seed_value)
        class(rng_state), intent(inout) :: self
        integer(kind=8), intent(in) :: seed_value
        integer(kind=8) :: s

        s = seed_value
        if (s == 0_8) s = 88172645463393265_8
        self%state = ieor(s, int(z'9E3779B97F4A7C15', kind=8))
        if (self%state == 0_8) self%state = 88172645463393265_8
        self%have_spare = .false.
        self%spare = 0.0_dp
    end subroutine rng_seed

    function next_u64(self) result(x)
        class(rng_state), intent(inout) :: self
        integer(kind=8) :: x

        x = self%state
        x = ieor(x, ishft(x, 13))
        x = ieor(x, ishft(x, -7))
        x = ieor(x, ishft(x, 17))
        self%state = x
    end function next_u64

    function rng_uniform(self) result(u)
        class(rng_state), intent(inout) :: self
        real(dp) :: u
        integer(kind=8) :: x
        integer(kind=8), parameter :: mask53 = int(z'001FFFFFFFFFFFFF', kind=8)

        x = iand(next_u64(self), mask53)
        u = (real(x, dp) + 0.5_dp) / real(mask53 + 1_8, dp)
    end function rng_uniform

    function rng_normal(self) result(z)
        class(rng_state), intent(inout) :: self
        real(dp) :: z
        real(dp) :: u1, u2, rad
        real(dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp

        if (self%have_spare) then
            z = self%spare
            self%have_spare = .false.
            return
        end if

        u1 = max(self%uniform(), tiny(1.0_dp))
        u2 = self%uniform()
        rad = sqrt(-2.0_dp * log(u1))
        z = rad * cos(twopi * u2)
        self%spare = rad * sin(twopi * u2)
        self%have_spare = .true.
    end function rng_normal

    function rng_randint(self, lo, hi) result(k)
        class(rng_state), intent(inout) :: self
        integer, intent(in) :: lo, hi
        integer :: k

        if (hi < lo) error stop "rng_randint: invalid interval"
        k = lo + int(self%uniform() * real(hi - lo + 1, dp))
        if (k > hi) k = hi
    end function rng_randint

    subroutine latin_hypercube(rng, x)
        type(rng_state), intent(inout) :: rng
        real(dp), intent(out) :: x(:,:)
        integer :: n, d, i, j, k, tmp
        integer, allocatable :: p(:)

        n = size(x, 2)
        d = size(x, 1)
        if (n <= 0 .or. d <= 0) return
        allocate(p(n))

        do j = 1, d
            p = [(i, i=1,n)]
            do i = n, 2, -1
                k = rng%randint(1, i)
                tmp = p(i)
                p(i) = p(k)
                p(k) = tmp
            end do
            do i = 1, n
                x(j,i) = (real(p(i)-1, dp) + rng%uniform()) / real(n, dp)
            end do
        end do
    end subroutine latin_hypercube

end module ppso_rng
