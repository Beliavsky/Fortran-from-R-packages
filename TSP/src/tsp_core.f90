! Modern Fortran translation of computational code from TSP 1.2.7.
! Original Copyright (C) Michael Hahsler and Kurt Hornik.
! SPDX-License-Identifier: GPL-3.0-only
! See LICENSE, COPYING, and UPSTREAM.md for provenance and licensing.

module tsp_core
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, &
        ieee_positive_inf, ieee_negative_inf, ieee_quiet_nan
    use tsp_kinds, only : dp
    implicit none
    private

    public :: is_square_matrix, is_symmetric_matrix, is_valid_tour
    public :: tour_length, etsp_tour_length, euclidean_distance_matrix
    public :: insertion_cost, replace_infinite
    public :: positive_infinity, negative_infinity, quiet_nan

contains

    pure logical function is_square_matrix(x) result(ok)
        real(dp), intent(in) :: x(:,:)
        ok = size(x, 1) == size(x, 2)
    end function is_square_matrix

    pure logical function is_symmetric_matrix(x, tol) result(ok)
        real(dp), intent(in) :: x(:,:)
        real(dp), intent(in), optional :: tol
        real(dp) :: eps
        integer :: i, j, n

        ok = .false.
        if (.not. is_square_matrix(x)) return
        eps = 1.0e-12_dp
        if (present(tol)) eps = tol
        n = size(x, 1)
        do j = 1, n
            do i = j + 1, n
                if (ieee_is_finite(x(i,j)) .neqv. ieee_is_finite(x(j,i))) return
                if (ieee_is_finite(x(i,j))) then
                    if (abs(x(i,j) - x(j,i)) > eps * max(1.0_dp, abs(x(i,j)), abs(x(j,i)))) return
                else
                    if (ieee_is_nan(x(i,j)) .or. ieee_is_nan(x(j,i))) return
                    if ((x(i,j) > 0.0_dp) .neqv. (x(j,i) > 0.0_dp)) return
                end if
            end do
        end do
        ok = .true.
    end function is_symmetric_matrix

    pure logical function is_valid_tour(order, n) result(ok)
        integer, intent(in) :: order(:)
        integer, intent(in) :: n
        logical, allocatable :: seen(:)
        integer :: i

        ok = .false.
        if (size(order) /= n) return
        allocate(seen(n), source=.false.)
        do i = 1, n
            if (order(i) < 1 .or. order(i) > n) return
            if (seen(order(i))) return
            seen(order(i)) = .true.
        end do
        ok = .true.
    end function is_valid_tour

    pure real(dp) function positive_infinity() result(x)
        x = ieee_value(0.0_dp, ieee_positive_inf)
    end function positive_infinity

    pure real(dp) function negative_infinity() result(x)
        x = ieee_value(0.0_dp, ieee_negative_inf)
    end function negative_infinity

    pure real(dp) function quiet_nan() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function quiet_nan

    pure real(dp) function tour_length(cost, order) result(total)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: order(:)
        integer :: i, n
        logical :: posinf, neginf
        real(dp) :: edge

        n = size(cost, 1)
        if (.not. is_square_matrix(cost) .or. .not. is_valid_tour(order, n)) then
            total = quiet_nan()
            return
        end if
        if (n <= 1) then
            total = 0.0_dp
            return
        end if

        total = 0.0_dp
        posinf = .false.
        neginf = .false.
        do i = 1, n - 1
            edge = cost(order(i), order(i+1))
            call accumulate_edge(edge, total, posinf, neginf)
        end do
        edge = cost(order(n), order(1))
        call accumulate_edge(edge, total, posinf, neginf)

        if (posinf .and. neginf) then
            total = quiet_nan()
        else if (posinf) then
            total = positive_infinity()
        else if (neginf) then
            total = negative_infinity()
        end if
    end function tour_length

    pure subroutine accumulate_edge(edge, total, posinf, neginf)
        real(dp), intent(in) :: edge
        real(dp), intent(inout) :: total
        logical, intent(inout) :: posinf, neginf

        if (ieee_is_nan(edge)) then
            total = quiet_nan()
        else if (.not. ieee_is_finite(edge)) then
            if (edge > 0.0_dp) then
                posinf = .true.
            else
                neginf = .true.
            end if
        else
            total = total + edge
        end if
    end subroutine accumulate_edge

    pure real(dp) function etsp_tour_length(coords, order) result(total)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: order(:)
        integer :: i, n

        n = size(coords, 1)
        if (.not. is_valid_tour(order, n)) then
            total = quiet_nan()
            return
        end if
        total = 0.0_dp
        if (n <= 1) return
        do i = 1, n - 1
            total = total + sqrt(sum((coords(order(i),:) - coords(order(i+1),:))**2))
        end do
        total = total + sqrt(sum((coords(order(n),:) - coords(order(1),:))**2))
    end function etsp_tour_length

    pure subroutine euclidean_distance_matrix(coords, cost)
        real(dp), intent(in) :: coords(:,:)
        real(dp), allocatable, intent(out) :: cost(:,:)
        integer :: i, j, n

        n = size(coords, 1)
        allocate(cost(n,n), source=0.0_dp)
        do j = 1, n
            do i = j + 1, n
                cost(i,j) = sqrt(sum((coords(i,:) - coords(j,:))**2))
                cost(j,i) = cost(i,j)
            end do
        end do
    end subroutine euclidean_distance_matrix

    pure subroutine insertion_cost(cost, order, k, delta)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: order(:)
        integer, intent(in) :: k
        real(dp), allocatable, intent(out) :: delta(:)
        integer :: i, m
        real(dp) :: a, b, c

        m = size(order)
        allocate(delta(m))
        if (m == 0) return
        if (m == 1) then
            delta(1) = cost(order(1), k)
            return
        end if

        do i = 1, m - 1
            a = cost(order(i), k)
            b = cost(k, order(i+1))
            c = cost(order(i), order(i+1))
            delta(i) = safe_add_remove(a, b, c, .false.)
        end do
        a = cost(order(m), k)
        b = cost(k, order(1))
        c = cost(order(m), order(1))
        delta(m) = safe_add_remove(a, b, c, .true.)
    end subroutine insertion_cost

    pure real(dp) function safe_add_remove(a, b, c, closing) result(v)
        real(dp), intent(in) :: a, b, c
        logical, intent(in) :: closing

        if (closing) then
            if ((.not. ieee_is_finite(a) .and. a > 0.0_dp) .or. &
                (.not. ieee_is_finite(b) .and. b > 0.0_dp)) then
                v = positive_infinity()
            else if (.not. ieee_is_finite(c) .and. c > 0.0_dp) then
                v = negative_infinity()
            else
                v = a + b - c
            end if
        else
            if ((.not. ieee_is_finite(a) .and. a < 0.0_dp) .or. &
                (.not. ieee_is_finite(b) .and. b < 0.0_dp) .or. &
                (.not. ieee_is_finite(c) .and. c > 0.0_dp)) then
                v = negative_infinity()
            else if ((.not. ieee_is_finite(a) .and. a > 0.0_dp) .or. &
                (.not. ieee_is_finite(b) .and. b > 0.0_dp) .or. &
                (.not. ieee_is_finite(c) .and. c < 0.0_dp)) then
                v = positive_infinity()
            else
                v = a + b - c
            end if
        end if
    end function safe_add_remove

    subroutine replace_infinite(x, y, positive_value, negative_value, ierr)
        real(dp), intent(in) :: x(:,:)
        real(dp), allocatable, intent(out) :: y(:,:)
        real(dp), intent(in), optional :: positive_value, negative_value
        integer, intent(out), optional :: ierr
        real(dp) :: xmin, xmax, diff_range, pvalue, nvalue
        integer :: i, j
        logical :: found

        if (present(ierr)) ierr = 0
        allocate(y, source=x)
        found = .false.
        do j = 1, size(x,2)
            do i = 1, size(x,1)
                if (.not. ieee_is_finite(x(i,j)) .and. .not. ieee_is_nan(x(i,j))) found = .true.
            end do
        end do
        if (.not. found) return

        found = .false.
        xmin = huge(1.0_dp)
        xmax = -huge(1.0_dp)
        do j = 1, size(x,2)
            do i = 1, size(x,1)
                if (ieee_is_finite(x(i,j))) then
                    found = .true.
                    xmin = min(xmin, x(i,j))
                    xmax = max(xmax, x(i,j))
                end if
            end do
        end do
        if (.not. found) then
            if (present(ierr)) ierr = 1
            return
        end if

        diff_range = xmax - xmin
        if (abs(diff_range) <= tiny(1.0_dp)) then
            if (abs(xmin) <= tiny(1.0_dp)) then
                diff_range = 1.0_dp
            else
                diff_range = abs(2.0_dp * xmin)
            end if
        end if
        pvalue = xmax + 2.0_dp * diff_range
        nvalue = xmin - 2.0_dp * diff_range
        if (present(positive_value)) pvalue = positive_value
        if (present(negative_value)) nvalue = negative_value

        do j = 1, size(y,2)
            do i = 1, size(y,1)
                if (.not. ieee_is_finite(y(i,j)) .and. .not. ieee_is_nan(y(i,j))) then
                    if (y(i,j) > 0.0_dp) then
                        y(i,j) = pvalue
                    else
                        y(i,j) = nvalue
                    end if
                end if
            end do
        end do
    end subroutine replace_infinite

end module tsp_core
