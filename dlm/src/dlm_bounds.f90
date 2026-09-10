! SPDX-License-Identifier: GPL-2.0-or-later
module dlm_bounds
    use dlm_types, only : dp, dlm_success, dlm_invalid_argument
    implicit none
    private

    type, abstract, public :: dlm_indicator
    contains
        procedure(dlm_indicator_evaluate), deferred, pass :: evaluate
    end type dlm_indicator

    abstract interface
        pure logical function dlm_indicator_evaluate(self, x) result(inside)
            import :: dlm_indicator, dp
            class(dlm_indicator), intent(in) :: self !! Indicator object holding any fixed support parameters.
            real(dp), intent(in) :: x(:) !! Point at which membership of the bounded convex support is tested.
        end function dlm_indicator_evaluate
    end interface

    public :: convex_bounds

contains

    pure subroutine convex_bounds(x, direction, indicator, bounds, info, tolerance, max_expansions)
        real(dp), intent(in) :: x(:) !! Point known to lie inside the bounded convex support.
        real(dp), intent(in) :: direction(:) !! Nonzero search direction with the same length as `x`.
        class(dlm_indicator), intent(in) :: indicator !! Convex-support membership object used during line searches.
        real(dp), intent(out) :: bounds(2) !! Signed lower and upper offsets `u` for points `x + u*direction`.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_argument` for invalid/unbounded searches.
        real(dp), intent(in), optional :: tolerance !! Absolute bisection tolerance; default `1e-7`.
        integer, intent(in), optional :: max_expansions !! Maximum doubling steps per side; default 100.

        real(dp) :: bracket_high, bracket_low, candidate, step, tol
        integer :: expansion_limit, i

        info = dlm_success
        bounds = 0.0_dp
        tol = 1.0e-7_dp
        if (present(tolerance)) tol = tolerance
        expansion_limit = 100
        if (present(max_expansions)) expansion_limit = max_expansions

        if (size(direction) /= size(x) .or. size(x) < 1) then
            info = dlm_invalid_argument
            return
        end if
        if (maxval(abs(direction)) <= tiny(1.0_dp) .or. tol <= 0.0_dp .or. expansion_limit < 1) then
            info = dlm_invalid_argument
            return
        end if
        if (.not. indicator%evaluate(x)) then
            info = dlm_invalid_argument
            return
        end if

        step = -2.0_dp
        do i = 1, expansion_limit
            if (.not. indicator%evaluate(x + step * direction)) exit
            step = 2.0_dp * step
        end do
        if (indicator%evaluate(x + step * direction)) then
            info = dlm_invalid_argument
            return
        end if
        bracket_low = step
        bracket_high = 0.0_dp
        do while (bracket_high - bracket_low >= tol)
            candidate = 0.5_dp * (bracket_low + bracket_high)
            if (indicator%evaluate(x + candidate * direction)) then
                bracket_high = candidate
            else
                bracket_low = candidate
            end if
        end do
        bounds(1) = bracket_high

        step = 2.0_dp
        do i = 1, expansion_limit
            if (.not. indicator%evaluate(x + step * direction)) exit
            step = 2.0_dp * step
        end do
        if (indicator%evaluate(x + step * direction)) then
            info = dlm_invalid_argument
            return
        end if
        bracket_low = 0.0_dp
        bracket_high = step
        do while (bracket_high - bracket_low >= tol)
            candidate = 0.5_dp * (bracket_low + bracket_high)
            if (indicator%evaluate(x + candidate * direction)) then
                bracket_low = candidate
            else
                bracket_high = candidate
            end if
        end do
        bounds(2) = bracket_low
    end subroutine convex_bounds

end module dlm_bounds
