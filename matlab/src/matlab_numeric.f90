! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module matlab_numeric
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use matlab_kinds, only : dp
    implicit none
    private

    public :: ceil
    public :: fix
    public :: matlab_mod
    public :: rem
    public :: nextpow2
    public :: pow2
    public :: pow2_scaled
    public :: linspace
    public :: logspace
    public :: std
    public :: std_cols
    public :: sum_cols

contains

    elemental function ceil(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = real(ceiling(x), dp)
    end function ceil

    elemental function fix(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = real(int(x), dp)
    end function fix

    elemental function matlab_mod(x, y) result(z)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp) :: z

        if (y == 0.0_dp) then
            z = x
        else
            z = modulo(x, y)
        end if
    end function matlab_mod

    elemental function rem(x, y) result(z)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y
        real(dp) :: z

        if (y == 0.0_dp) then
            z = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            z = mod(x, y)
        end if
    end function rem

    elemental function nextpow2(x) result(p)
        real(dp), intent(in) :: x
        integer :: p
        real(dp) :: ax

        ax = abs(x)
        if (ax == 0.0_dp) ax = 1.0_dp
        p = ceiling(log(ax) / log(2.0_dp))
    end function nextpow2

    elemental function pow2(e) result(y)
        real(dp), intent(in) :: e
        real(dp) :: y

        y = 2.0_dp ** e
    end function pow2

    elemental function pow2_scaled(f, e) result(y)
        real(dp), intent(in) :: f
        real(dp), intent(in) :: e
        real(dp) :: y

        y = f * 2.0_dp ** e
    end function pow2_scaled

    function linspace(a, b, n) result(x)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        integer, intent(in), optional :: n
        real(dp), allocatable :: x(:)
        integer :: i, m

        m = 100
        if (present(n)) m = n
        if (m < 2) then
            allocate(x(1))
            x(1) = b
            return
        end if

        allocate(x(m))
        do i = 1, m
            x(i) = a + real(i - 1, dp) * (b - a) / real(m - 1, dp)
        end do
        x(m) = b
    end function linspace

    function logspace(a, b, n) result(x)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        integer, intent(in), optional :: n
        real(dp), allocatable :: x(:)
        real(dp) :: bb, pi
        integer :: i, m

        pi = acos(-1.0_dp)
        bb = b
        if (b == pi) bb = log10(pi)
        m = 50
        if (present(n)) m = n
        if (m < 2) then
            allocate(x(1))
            x(1) = 10.0_dp ** bb
            return
        end if
        allocate(x(m))
        do i = 1, m
            x(i) = 10.0_dp ** (a + real(i - 1, dp) * (bb - a) / real(m - 1, dp))
        end do
        x(m) = 10.0_dp ** bb
    end function logspace

    function std(x, flag) result(s)
        real(dp), intent(in) :: x(:)
        integer, intent(in), optional :: flag
        real(dp) :: s
        real(dp) :: mean_x
        integer :: denominator, f

        f = 0
        if (present(flag)) f = flag
        if (size(x) == 0) then
            s = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        denominator = size(x) - 1
        if (f /= 0) denominator = size(x)
        if (denominator <= 0) then
            s = 0.0_dp
            return
        end if

        mean_x = sum(x) / real(size(x), dp)
        s = sqrt(sum((x - mean_x) ** 2) / real(denominator, dp))
    end function std

    function std_cols(a, flag) result(s)
        real(dp), intent(in) :: a(:, :)
        integer, intent(in), optional :: flag
        real(dp), allocatable :: s(:)
        integer :: j

        allocate(s(size(a, 2)))
        do j = 1, size(a, 2)
            if (present(flag)) then
                s(j) = std(a(:, j), flag)
            else
                s(j) = std(a(:, j))
            end if
        end do
    end function std_cols

    function sum_cols(a) result(s)
        real(dp), intent(in) :: a(:, :)
        real(dp), allocatable :: s(:)
        integer :: j

        allocate(s(size(a, 2)))
        do j = 1, size(a, 2)
            s(j) = sum(a(:, j))
        end do
    end function sum_cols
end module matlab_numeric
