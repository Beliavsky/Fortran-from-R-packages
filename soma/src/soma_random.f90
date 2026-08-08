! SPDX-License-Identifier: GPL-2.0-only
module soma_random
    use iso_fortran_env, only : int64
    use soma_kinds, only : dp
    implicit none
    private
    public :: soma_set_seed, random_index, sample_without_replacement

contains

    subroutine soma_set_seed(seed)
        integer, intent(in) :: seed
        integer :: i, n
        integer, allocatable :: put(:)
        integer(int64) :: state

        call random_seed(size=n)
        allocate(put(n))
        state = modulo(int(seed, int64), 2147483646_int64) + 1_int64
        do i = 1, n
            state = modulo(1664525_int64 * state + 1013904223_int64, 2147483647_int64)
            if (state == 0_int64) state = int(i, int64)
            put(i) = int(state, kind(put))
        end do
        call random_seed(put=put)
    end subroutine soma_set_seed

    integer function random_index(n) result(index_value)
        integer, intent(in) :: n
        real(dp) :: u

        if (n <= 1) then
            index_value = 1
            return
        end if
        call random_number(u)
        index_value = 1 + int(u * real(n, dp))
        if (index_value > n) index_value = n
    end function random_index

    subroutine sample_without_replacement(n, k, sample)
        integer, intent(in) :: n, k
        integer, intent(out) :: sample(k)
        integer, allocatable :: work(:)
        integer :: i, j, tmp

        allocate(work(n))
        work = [(i, i=1,n)]
        do i = 1, k
            j = i - 1 + random_index(n - i + 1)
            tmp = work(i)
            work(i) = work(j)
            work(j) = tmp
            sample(i) = work(i)
        end do
    end subroutine sample_without_replacement

end module soma_random
