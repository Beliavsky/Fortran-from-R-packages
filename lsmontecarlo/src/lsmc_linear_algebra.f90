! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_linear_algebra
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use lsmc_kinds, only : dp
    implicit none
    private

    public :: least_squares

contains

    subroutine least_squares(x, y, beta, ok)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: beta(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: column_norms(:)
        real(dp), allocatable :: q(:, :)
        real(dp), allocatable :: r(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: scales(:)
        real(dp), allocatable :: solution(:)
        real(dp), allocatable :: v(:)
        real(dp) :: norm_v
        real(dp) :: temp
        real(dp) :: tolerance
        integer, allocatable :: permutation(:)
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: p
        integer :: pivot
        integer :: rank
        integer :: temp_index

        n = size(x, 1)
        p = size(x, 2)
        if (size(y) /= n .or. size(beta) /= p) error stop "least_squares: inconsistent dimensions"
        beta = 0.0_dp
        if (n == 0 .or. p == 0) then
            ok = .false.
            return
        end if

        allocate(a(n, p), q(n, p), r(p, p), rhs(p), scales(p), solution(p), v(n))
        allocate(column_norms(p), permutation(p))
        q = 0.0_dp
        r = 0.0_dp
        rhs = 0.0_dp
        solution = 0.0_dp
        do j = 1, p
            scales(j) = sqrt(sum(x(:, j)**2))
            if (scales(j) > tiny(1.0_dp)) then
                a(:, j) = x(:, j) / scales(j)
            else
                a(:, j) = 0.0_dp
                scales(j) = 1.0_dp
            end if
            column_norms(j) = sum(a(:, j)**2)
            permutation(j) = j
        end do

        tolerance = 100.0_dp * epsilon(1.0_dp) * real(max(n, p), dp)
        rank = 0
        do k = 1, p
            pivot = k - 1 + maxloc(column_norms(k:p), dim=1)
            if (pivot /= k) then
                call swap_columns(a, k, pivot)
                temp = column_norms(k)
                column_norms(k) = column_norms(pivot)
                column_norms(pivot) = temp
                temp_index = permutation(k)
                permutation(k) = permutation(pivot)
                permutation(pivot) = temp_index
                do i = 1, k - 1
                    temp = r(i, k)
                    r(i, k) = r(i, pivot)
                    r(i, pivot) = temp
                end do
            end if

            v = a(:, k)
            norm_v = sqrt(sum(v**2))
            if (norm_v <= tolerance) exit
            rank = k
            q(:, k) = v / norm_v
            r(k, k) = norm_v
            do j = k + 1, p
                r(k, j) = dot_product(q(:, k), a(:, j))
                a(:, j) = a(:, j) - r(k, j) * q(:, k)
                column_norms(j) = sum(a(:, j)**2)
            end do
        end do

        if (rank == 0) then
            ok = .false.
            return
        end if

        rhs(1:rank) = matmul(transpose(q(:, 1:rank)), y)
        solution(rank) = rhs(rank) / r(rank, rank)
        do i = rank - 1, 1, -1
            solution(i) = (rhs(i) - dot_product(r(i, i + 1:rank), solution(i + 1:rank))) / r(i, i)
        end do
        do i = 1, rank
            beta(permutation(i)) = solution(i) / scales(permutation(i))
        end do
        ok = all(ieee_is_finite(beta))
        if (.not. ok) beta = 0.0_dp
    end subroutine least_squares

    subroutine swap_columns(a, first, second)
        real(dp), intent(inout) :: a(:, :)
        integer, intent(in) :: first
        integer, intent(in) :: second
        real(dp), allocatable :: temporary(:)

        allocate(temporary(size(a, 1)))
        temporary = a(:, first)
        a(:, first) = a(:, second)
        a(:, second) = temporary
    end subroutine swap_columns

end module lsmc_linear_algebra
