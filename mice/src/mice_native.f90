! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Translation of mice src/legendre.cpp.
module mice_native
    use r_kinds, only : dp
    use mice_status, only : mice_ok, mice_invalid_argument
    implicit none
    private

    public :: legendre_basis

contains

    pure subroutine legendre_basis(x, degree, basis, info)
        real(dp), intent(in) :: x(:) !! Input values, conventionally scaled to `[0,1]` as in upstream `legendre()`.
        integer, intent(in), value :: degree !! Number of normalized Legendre columns to generate; must be at least one.
        real(dp), allocatable, intent(out) :: basis(:, :) !! Matrix with `size(x)` rows and `degree` normalized polynomial columns.
        integer, intent(out) :: info !! `mice_ok` on success or `mice_invalid_argument` for nonpositive degree.

        real(dp), allocatable :: z(:), raw(:, :)
        integer :: j

        if (degree < 1) then
            info = mice_invalid_argument
            return
        end if
        allocate(basis(size(x), degree), z(size(x)), raw(size(x), degree))
        z = 2.0_dp * x - 1.0_dp
        raw(:, 1) = z
        if (degree >= 2) raw(:, 2) = (3.0_dp * z * z - 1.0_dp) / 2.0_dp
        do j = 3, degree
            raw(:, j) = (real(2 * j - 1, dp) * z * raw(:, j - 1) - &
                         real(j - 1, dp) * raw(:, j - 2)) / real(j, dp)
        end do
        do j = 1, degree
            basis(:, j) = sqrt(real(2 * j + 1, dp)) * raw(:, j)
        end do
        info = mice_ok
    end subroutine legendre_basis

end module mice_native
