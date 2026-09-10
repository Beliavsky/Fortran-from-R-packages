! SPDX-License-Identifier: GPL-2.0-or-later
module test_bounds_support
    use dlm, only : dlm_indicator, dp
    implicit none
    private

    type, extends(dlm_indicator), public :: box_indicator
        real(dp) :: half_width(2) = [2.0_dp, 3.0_dp]
    contains
        procedure :: evaluate => box_evaluate
    end type box_indicator

contains

    pure logical function box_evaluate(self, x) result(inside)
        class(box_indicator), intent(in) :: self !! Box indicator carrying half-width limits for each coordinate.
        real(dp), intent(in) :: x(:) !! Candidate two-dimensional point tested for box membership.

        inside = size(x) == 2
        if (inside) inside = all(abs(x) <= self%half_width)
    end function box_evaluate

end module test_bounds_support

program test_bounds
    use dlm, only : convex_bounds, dlm_success, dp
    use test_bounds_support, only : box_indicator
    implicit none

    type(box_indicator) :: box
    real(dp) :: bounds(2), direction(2), x(2)
    integer :: info

    x = [0.5_dp, -0.5_dp]
    direction = [1.0_dp, 0.0_dp]
    call convex_bounds(x, direction, box, bounds, info, tolerance=1.0e-10_dp)
    if (info /= dlm_success) error stop "convex_bounds failed"
    if (abs(bounds(1) + 2.5_dp) > 2.0e-10_dp) error stop "convex_bounds lower failed"
    if (abs(bounds(2) - 1.5_dp) > 2.0e-10_dp) error stop "convex_bounds upper failed"

    direction = [1.0_dp, 1.0_dp]
    call convex_bounds(x, direction, box, bounds, info, tolerance=1.0e-10_dp)
    if (info /= dlm_success) error stop "convex_bounds diagonal failed"
    if (abs(bounds(1) + 2.5_dp) > 2.0e-10_dp) error stop "convex_bounds diagonal lower failed"
    if (abs(bounds(2) - 1.5_dp) > 2.0e-10_dp) error stop "convex_bounds diagonal upper failed"

    print *, "test_bounds: PASS"
end program test_bounds
