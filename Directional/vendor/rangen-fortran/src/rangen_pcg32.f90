module rangen_pcg32
    use rangen_kinds, only : dp, i8
    implicit none
    private

    integer(i8), parameter :: mask16 = int(z'FFFF', i8)
    integer(i8), parameter :: mask32 = int(z'FFFFFFFF', i8)
    integer(i8), parameter :: base16 = 65536_i8
    integer(i8), parameter :: mul_limb(0:3) = [32557_i8, 19605_i8, 62509_i8, 22609_i8]

    type, public :: pcg32_state
        integer(i8) :: state = 0_i8
        integer(i8) :: inc = 1442695040888963407_i8
    contains
        procedure :: seed => pcg_seed
        procedure :: next_uint32 => pcg_next_uint32
        procedure :: uniform_closed => pcg_uniform_closed
        procedure :: uniform_open => pcg_uniform_open
    end type pcg32_state

contains

    subroutine pcg_seed(self, seed_value, stream)
        class(pcg32_state), intent(inout) :: self
        integer(i8), intent(in) :: seed_value
        integer(i8), intent(in), optional :: stream

        self%state = seed_value
        if (present(stream)) then
            self%inc = ior(stream, 1_i8)
        else
            self%inc = 1442695040888963407_i8
        end if
    end subroutine pcg_seed

    function pcg_next_uint32(self) result(out)
        class(pcg32_state), intent(inout) :: self
        integer(i8) :: out
        integer(i8) :: oldstate, xorshifted
        integer :: rot

        oldstate = self%state
        self%state = pcg_muladd(oldstate, self%inc)
        xorshifted = iand(shiftr(ieor(shiftr(oldstate, 18), oldstate), 27), mask32)
        rot = int(iand(shiftr(oldstate, 59), 31_i8))
        if (rot == 0) then
            out = xorshifted
        else
            out = iand(ior(shiftr(xorshifted, rot), shiftl(xorshifted, 32 - rot)), mask32)
        end if
    end function pcg_next_uint32

    function pcg_uniform_closed(self) result(u)
        class(pcg32_state), intent(inout) :: self
        real(dp) :: u
        u = real(self%next_uint32(), dp) / 4294967295.0_dp
    end function pcg_uniform_closed

    function pcg_uniform_open(self) result(u)
        class(pcg32_state), intent(inout) :: self
        real(dp) :: u
        u = (real(self%next_uint32(), dp) + 0.5_dp) / 4294967296.0_dp
    end function pcg_uniform_open

    pure function pcg_muladd(x, inc) result(y)
        integer(i8), intent(in) :: x, inc
        integer(i8) :: y
        integer(i8) :: s(0:3), q(0:3), r(0:3)
        integer(i8) :: total, carry
        integer :: k, j

        do k = 0, 3
            s(k) = iand(shiftr(x, 16 * k), mask16)
            q(k) = iand(shiftr(inc, 16 * k), mask16)
        end do

        carry = 0_i8
        do k = 0, 3
            total = carry + q(k)
            do j = 0, k
                total = total + s(j) * mul_limb(k - j)
            end do
            r(k) = modulo(total, base16)
            carry = total / base16
        end do

        y = r(0)
        y = ior(y, shiftl(r(1), 16))
        y = ior(y, shiftl(r(2), 32))
        y = ior(y, shiftl(r(3), 48))
    end function pcg_muladd

end module rangen_pcg32
