! Modern Fortran translation of computational code from TSP 1.2.7.
! Original Copyright (C) Michael Hahsler and Kurt Hornik.
! SPDX-License-Identifier: GPL-3.0-only
! See LICENSE, COPYING, and UPSTREAM.md for provenance and licensing.

module tsp_transform
    use tsp_kinds, only : dp
    use tsp_types, only : tsp_path_collection
    use tsp_core, only : tour_length, positive_infinity, negative_infinity
    implicit none
    private

    public :: insert_dummy, reformulate_atsp_as_tsp, filter_atsp_tour
    public :: cut_tour_single, cut_tour_multiple

contains

    subroutine insert_dummy(cost, n_dummy, output, const, dummy_distance)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: n_dummy
        real(dp), allocatable, intent(out) :: output(:,:)
        real(dp), intent(in), optional :: const, dummy_distance
        real(dp) :: c, inf
        integer :: n, i, j

        n = size(cost,1)
        c = 0.0_dp
        inf = positive_infinity()
        if (present(const)) c = const
        if (present(dummy_distance)) inf = dummy_distance
        if (n_dummy <= 0) then
            allocate(output, source=cost)
            return
        end if

        allocate(output(n+n_dummy,n+n_dummy), source=c)
        output(:n,:n) = cost
        do i = n + 1, n + n_dummy
            output(i,i) = 0.0_dp
        end do
        if (n_dummy > 1) then
            do j = n + 1, n + n_dummy
                do i = n + 1, n + n_dummy
                    if (i /= j) output(i,j) = inf
                end do
            end do
        end if
    end subroutine insert_dummy

    subroutine reformulate_atsp_as_tsp(atsp, tsp, infeasible, cheap, ierr)
        real(dp), intent(in) :: atsp(:,:)
        real(dp), allocatable, intent(out) :: tsp(:,:)
        real(dp), intent(in), optional :: infeasible, cheap
        integer, intent(out), optional :: ierr
        real(dp), allocatable :: m(:,:)
        real(dp) :: bad, low, min_off, max_value
        integer :: i, j, n

        if (present(ierr)) ierr = 0
        n = size(atsp,1)
        if (size(atsp,2) /= n) then
            allocate(tsp(0,0))
            if (present(ierr)) ierr = 1
            return
        end if
        bad = positive_infinity()
        low = negative_infinity()
        if (present(infeasible)) bad = infeasible
        if (present(cheap)) low = cheap

        min_off = huge(1.0_dp)
        max_value = -huge(1.0_dp)
        do j = 1, n
            do i = 1, n
                max_value = max(max_value, atsp(i,j))
                if (i /= j) min_off = min(min_off, atsp(i,j))
            end do
        end do
        if (low >= min_off) then
            allocate(tsp(0,0))
            if (present(ierr)) ierr = 2
            return
        end if
        if (bad < max_value) then
            allocate(tsp(0,0))
            if (present(ierr)) ierr = 3
            return
        end if

        allocate(m, source=atsp)
        do i = 1, n
            m(i,i) = low
        end do
        allocate(tsp(2*n,2*n), source=bad)
        tsp(1:n,n+1:2*n) = transpose(m)
        tsp(n+1:2*n,1:n) = m
    end subroutine reformulate_atsp_as_tsp

    subroutine filter_atsp_tour(transformed_order, atsp, order)
        integer, intent(in) :: transformed_order(:)
        real(dp), intent(in) :: atsp(:,:)
        integer, allocatable, intent(out) :: order(:)
        integer, allocatable :: forward(:), backward(:)
        integer :: i, k, n

        n = size(atsp,1)
        allocate(forward(n))
        k = 0
        do i = 1, size(transformed_order)
            if (transformed_order(i) <= n) then
                k = k + 1
                if (k <= n) forward(k) = transformed_order(i)
            end if
        end do
        if (k /= n) then
            allocate(order(0))
            return
        end if
        allocate(backward(n))
        backward = forward(n:1:-1)
        if (tour_length(atsp, backward) < tour_length(atsp, forward)) then
            allocate(order, source=backward)
        else
            allocate(order, source=forward)
        end if
    end subroutine filter_atsp_tour

    subroutine cut_tour_single(tour, cut_city, path, exclude_cut, ierr)
        integer, intent(in) :: tour(:)
        integer, intent(in) :: cut_city
        integer, allocatable, intent(out) :: path(:)
        logical, intent(in), optional :: exclude_cut
        integer, intent(out), optional :: ierr
        integer :: pos, i, n, start, m
        logical :: exclude

        if (present(ierr)) ierr = 0
        exclude = .true.
        if (present(exclude_cut)) exclude = exclude_cut
        n = size(tour)
        pos = 0
        do i = 1, n
            if (tour(i) == cut_city) then
                pos = i
                exit
            end if
        end do
        if (pos == 0) then
            allocate(path(0))
            if (present(ierr)) ierr = 1
            return
        end if
        if (exclude) then
            m = max(0, n - 1)
            start = pos + 1
        else
            m = n
            start = pos
        end if
        allocate(path(m))
        do i = 1, m
            path(i) = tour(1 + modulo(start + i - 2, n))
        end do
    end subroutine cut_tour_single

    subroutine cut_tour_multiple(tour, cut_cities, paths, exclude_cut, ierr)
        integer, intent(in) :: tour(:)
        integer, intent(in) :: cut_cities(:)
        type(tsp_path_collection), intent(out) :: paths
        logical, intent(in), optional :: exclude_cut
        integer, intent(out), optional :: ierr
        integer, allocatable :: cut_pos(:), sorted(:)
        integer :: n, nc, i, j, p, q, len, tmp, out_index
        logical :: exclude

        if (present(ierr)) ierr = 0
        exclude = .true.
        if (present(exclude_cut)) exclude = exclude_cut
        n = size(tour)
        nc = size(cut_cities)
        allocate(paths%path(nc), cut_pos(nc), sorted(nc))
        do i = 1, nc
            cut_pos(i) = 0
            do j = 1, n
                if (tour(j) == cut_cities(i)) then
                    cut_pos(i) = j
                    exit
                end if
            end do
            if (cut_pos(i) == 0) then
                if (present(ierr)) ierr = 1
                return
            end if
        end do
        sorted = cut_pos
        do i = 2, nc
            tmp = sorted(i)
            j = i - 1
            do while (j >= 1)
                if (sorted(j) <= tmp) exit
                sorted(j+1) = sorted(j)
                j = j - 1
            end do
            sorted(j+1) = tmp
        end do

        ! R's cut_tour returns the wrapped segment first: last cut -> first cut,
        ! followed by first -> second, second -> third, and so on.
        do i = 0, nc - 1
            if (i == 0) then
                p = sorted(nc)
                q = sorted(1) + n
                out_index = 1
            else
                p = sorted(i)
                q = sorted(i+1)
                out_index = i + 1
            end if
            if (exclude) then
                len = max(0, q - p - 1)
                allocate(paths%path(out_index)%city(len))
                do j = 1, len
                    paths%path(out_index)%city(j) = tour(1 + modulo(p + j - 1, n))
                end do
            else
                len = q - p
                allocate(paths%path(out_index)%city(len))
                do j = 1, len
                    paths%path(out_index)%city(j) = tour(1 + modulo(p + j - 2, n))
                end do
            end if
        end do
    end subroutine cut_tour_multiple

end module tsp_transform
