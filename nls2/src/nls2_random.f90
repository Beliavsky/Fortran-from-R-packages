! SPDX-License-Identifier: GPL-2.0-only
module nls2_random
    use nls2_kinds, only : dp
    implicit none
    private
    public :: seed_rng, random_uniform, latin_hypercube

contains

    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)
        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729 * (i - 1), huge(1) - 1)
            if (put(i) <= 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine seed_rng

    subroutine random_uniform(lower, upper, points)
        real(dp), intent(in) :: lower(:), upper(:)
        real(dp), intent(out) :: points(:,:)
        real(dp), allocatable :: u(:,:)
        integer :: j
        allocate(u(size(points,1), size(points,2)))
        call random_number(u)
        do j = 1, size(points,2)
            points(:,j) = lower(j) + u(:,j) * (upper(j) - lower(j))
        end do
    end subroutine random_uniform

    subroutine latin_hypercube(lower, upper, points)
        real(dp), intent(in) :: lower(:), upper(:)
        real(dp), intent(out) :: points(:,:)
        integer :: n, p, i, j, k
        integer, allocatable :: perm(:)
        real(dp), allocatable :: u(:)
        real(dp) :: r

        n = size(points,1)
        p = size(points,2)
        allocate(perm(n), u(n))
        do j = 1, p
            do i = 1, n
                perm(i) = i
            end do
            do i = n, 2, -1
                call random_number(r)
                k = 1 + int(r * real(i, dp))
                if (k > i) k = i
                call swap_int(perm(i), perm(k))
            end do
            call random_number(u)
            do i = 1, n
                points(i,j) = lower(j) + &
                    ((real(perm(i)-1,dp) + u(i)) / real(n,dp)) * (upper(j)-lower(j))
            end do
        end do
    end subroutine latin_hypercube

    subroutine swap_int(a, b)
        integer, intent(inout) :: a, b
        integer :: t
        t = a
        a = b
        b = t
    end subroutine swap_int

end module nls2_random
