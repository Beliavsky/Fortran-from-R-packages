! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_spline
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use r_kinds, only : dp
    use r_linalg, only : symmetric_eigen
    implicit none
    private

    public :: spline_lrtp

contains

    subroutine spline_lrtp(x, basis, info, k, knots)
        real(dp), intent(in) :: x(:) !! Scalar predictor values at which the low-rank thin-plate spline basis is evaluated.
        real(dp), allocatable, intent(out) :: basis(:, :) !! Observation-by-knot low-rank thin-plate spline basis matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid knot settings, missing data, or eigensolver failure.
        integer, intent(in), optional :: k !! Number of default quantile knots; default is 10 when knots is absent.
        real(dp), intent(in), optional :: knots(:) !! Optional explicit knot positions; when present, k is ignored.
        real(dp), allocatable :: b_matrix(:, :)
        real(dp), allocatable :: eigenvalues(:)
        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: knot_values(:)
        real(dp), allocatable :: transform(:, :)
        real(dp), allocatable :: z_matrix(:, :)
        real(dp) :: coefficient
        real(dp) :: tolerance
        integer :: i
        integer :: j
        integer :: knot_count
        integer :: status

        info = 0
        if (size(x) < 1) then
            allocate(basis(0, 0))
            info = 1
            return
        end if
        if (present(knots)) then
            knot_count = size(knots)
            if (knot_count < 1 .or. any(ieee_is_nan(knots))) then
                allocate(basis(0, 0))
                info = 2
                return
            end if
            knot_values = knots
        else
            knot_count = 10
            if (present(k)) knot_count = k
            if (knot_count < 1) then
                allocate(basis(0, 0))
                info = 2
                return
            end if
            call default_quantile_knots(x, knot_count, knot_values, status)
            if (status /= 0) then
                allocate(basis(0, 0))
                info = status
                return
            end if
        end if
        allocate(z_matrix(size(x), knot_count), b_matrix(knot_count, knot_count))
        do j = 1, knot_count
            do i = 1, size(x)
                if (ieee_is_nan(x(i))) then
                    z_matrix(i, j) = 0.0_dp
                else
                    z_matrix(i, j) = abs(x(i) - knot_values(j)) ** 3
                end if
            end do
            do i = 1, knot_count
                b_matrix(i, j) = abs(knot_values(i) - knot_values(j)) ** 3
            end do
        end do
        call symmetric_eigen(b_matrix, eigenvalues, eigenvectors, status, descending=.true.)
        if (status /= 0) then
            allocate(basis(0, 0))
            info = 3
            return
        end if
        tolerance = sqrt(epsilon(1.0_dp))
        allocate(transform(knot_count, knot_count))
        transform = 0.0_dp
        do i = 1, knot_count
            if (abs(eigenvalues(i)) <= tolerance) cycle
            coefficient = sign(1.0_dp, eigenvalues(i)) / sqrt(abs(eigenvalues(i)))
            do j = 1, knot_count
                transform(:, j) = transform(:, j) + coefficient * eigenvectors(:, i) * eigenvectors(j, i)
            end do
        end do
        basis = matmul(z_matrix, transform)
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) basis(i, :) = 0.0_dp
        end do
    end subroutine spline_lrtp

    pure subroutine default_quantile_knots(x, knot_count, knots, info)
        real(dp), intent(in) :: x(:) !! Predictor values; NaNs are omitted to match R quantile(..., na.rm=TRUE) behavior.
        integer, intent(in) :: knot_count !! Positive number of equally spaced interior probability knots.
        real(dp), allocatable, intent(out) :: knots(:) !! R type-7 quantiles at probabilities i/(k+1), i=1,...,k.
        integer, intent(out) :: info !! Zero on success; nonzero when no non-NaN predictor values are available.
        real(dp), allocatable :: finite_values(:)
        real(dp) :: probability
        integer :: count
        integer :: i

        count = count_non_nan(x)
        if (count < 1 .or. knot_count < 1) then
            allocate(knots(0))
            info = 1
            return
        end if
        allocate(finite_values(count))
        count = 0
        do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            count = count + 1
            finite_values(count) = x(i)
        end do
        call insertion_sort(finite_values)
        allocate(knots(knot_count))
        do i = 1, knot_count
            probability = real(i, dp) / real(knot_count + 1, dp)
            knots(i) = quantile_type7_sorted(finite_values, probability)
        end do
        info = 0
    end subroutine default_quantile_knots

    pure integer function count_non_nan(x) result(count)
        real(dp), intent(in) :: x(:) !! Numeric vector whose non-NaN entries are counted.
        integer :: i

        count = 0
        do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) count = count + 1
        end do
    end function count_non_nan

    pure subroutine insertion_sort(values)
        real(dp), intent(inout) :: values(:) !! Numeric vector sorted increasingly in place.
        real(dp) :: key
        integer :: i
        integer :: j

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine insertion_sort

    pure real(dp) function quantile_type7_sorted(values, probability) result(value)
        real(dp), intent(in) :: values(:) !! Increasing finite data values.
        real(dp), intent(in) :: probability !! Quantile probability in the closed interval [0,1].
        real(dp) :: fraction
        real(dp) :: h
        integer :: lower
        integer :: n

        n = size(values)
        if (probability <= 0.0_dp .or. n == 1) then
            value = values(1)
            return
        end if
        if (probability >= 1.0_dp) then
            value = values(n)
            return
        end if
        h = 1.0_dp + real(n - 1, dp) * probability
        lower = floor(h)
        fraction = h - real(lower, dp)
        if (lower >= n) then
            value = values(n)
        else
            value = values(lower) + fraction * (values(lower + 1) - values(lower))
        end if
    end function quantile_type7_sorted

end module mcmcglmm_spline
