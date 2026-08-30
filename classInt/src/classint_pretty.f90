! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_pretty
    use classint_kinds, only: dp
    implicit none
    private

    public :: pretty_breaks

contains

    pure function pretty_breaks(lo_in, hi_in, n) result(brks)
        !! Constructs aesthetically spaced breaks covering a finite numeric range.
        real(dp), intent(in) :: lo_in !! Lower endpoint covered by the generated breaks.
        real(dp), intent(in) :: hi_in !! Upper endpoint of the finite data range; must not be below lo_in.
        integer, intent(in) :: n !! Approximate requested number of intervals; values below one are promoted to one.
        real(dp), allocatable :: brks(:) !! Ascending break vector covering `lo_in:hi_in`.
        real(dp), parameter :: high_u_bias = 1.5_dp
        real(dp), parameter :: u5_bias = 2.75_dp
        real(dp), parameter :: shrink_small = 0.75_dp
        real(dp), parameter :: rounding_eps = 1.0e-10_dp
        real(dp), parameter :: f_min = 2.0_dp**(-20)
        real(dp) :: base, candidate, cell, dx, hi, lo, ns, nu, subsmall, unit, ulimit
        integer :: i, k, ndiv
        logical :: small

        if (hi_in < lo_in) error stop "pretty_breaks: decreasing range"
        lo = lo_in
        hi = hi_in
        ndiv = max(1, n)
        dx = hi - lo
        if (abs(dx) <= tiny(1.0_dp) .and. abs(hi) <= tiny(1.0_dp)) then
            cell = 1.0_dp
            small = .true.
        else
            cell = max(abs(lo), abs(hi))
            ulimit = 1.0_dp + 1.0_dp / (1.0_dp + high_u_bias)
            ulimit = ulimit * real(ndiv, dp) * epsilon(1.0_dp)
            small = dx < 3.0_dp * cell * ulimit
        end if
        if (small) then
            if (cell > 10.0_dp) cell = 9.0_dp + cell / 10.0_dp
            cell = cell * shrink_small
        else
            cell = dx / real(ndiv, dp)
        end if
        subsmall = f_min * tiny(1.0_dp)
        if (subsmall <= 0.0_dp) subsmall = tiny(1.0_dp)
        cell = max(cell, subsmall)
        cell = min(cell, huge(1.0_dp) / 1.25_dp)
        base = 10.0_dp**floor(log10(cell))
        unit = base
        candidate = 2.0_dp * base
        if (candidate - cell < high_u_bias * (cell - unit)) then
            unit = candidate
            candidate = 5.0_dp * base
            if (candidate - cell < u5_bias * (cell - unit)) then
                unit = candidate
                candidate = 10.0_dp * base
                if (candidate - cell < high_u_bias * (cell - unit)) unit = candidate
            end if
        end if
        ns = floor(lo / unit + rounding_eps)
        nu = ceiling(hi / unit - rounding_eps)
        do while (ns * unit > lo + rounding_eps * unit)
            ns = ns - 1.0_dp
        end do
        do while (nu * unit < hi - rounding_eps * unit)
            nu = nu + 1.0_dp
        end do
        k = max(1, nint(nu - ns))
        allocate (brks(k + 1))
        do i = 0, k
            brks(i + 1) = (ns + real(i, dp)) * unit
        end do
    end function pretty_breaks
end module classint_pretty
