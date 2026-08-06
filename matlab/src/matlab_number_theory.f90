! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module matlab_number_theory
    use, intrinsic :: iso_fortran_env, only : int64
    implicit none
    private

    public :: primes
    public :: isprime
    public :: factors

contains

    function primes(n) result(p)
        integer, intent(in) :: n
        integer, allocatable :: p(:)
        logical, allocatable :: candidate(:)
        integer :: i, j, count_p

        if (n < 2) then
            allocate(p(0))
            return
        end if
        allocate(candidate(n))
        candidate = .true.
        candidate(1) = .false.
        do i = 2, int(sqrt(real(n)))
            if (candidate(i)) then
                do j = i * i, n, i
                    candidate(j) = .false.
                end do
            end if
        end do
        count_p = count(candidate)
        allocate(p(count_p))
        j = 0
        do i = 2, n
            if (candidate(i)) then
                j = j + 1
                p(j) = i
            end if
        end do
    end function primes

    elemental function isprime(n) result(answer)
        integer(int64), intent(in) :: n
        logical :: answer
        integer(int64) :: d

        if (n < 2_int64) then
            answer = .false.
            return
        end if
        if (n == 2_int64) then
            answer = .true.
            return
        end if
        if (modulo(n, 2_int64) == 0_int64) then
            answer = .false.
            return
        end if
        d = 3_int64
        do while (d <= n / d)
            if (modulo(n, d) == 0_int64) then
                answer = .false.
                return
            end if
            d = d + 2_int64
        end do
        answer = .true.
    end function isprime

    function factors(n) result(f)
        integer(int64), intent(in) :: n
        integer(int64), allocatable :: f(:)
        integer(int64), allocatable :: work(:)
        integer(int64) :: remaining, d
        integer :: count_f

        if (n < 4_int64) then
            allocate(f(1))
            f(1) = n
            return
        end if
        allocate(work(64))
        count_f = 0
        remaining = n
        d = 2_int64
        do while (d <= remaining / d)
            do while (modulo(remaining, d) == 0_int64)
                count_f = count_f + 1
                work(count_f) = d
                remaining = remaining / d
            end do
            if (d == 2_int64) then
                d = 3_int64
            else
                d = d + 2_int64
            end if
        end do
        if (remaining > 1_int64) then
            count_f = count_f + 1
            work(count_f) = remaining
        end if
        allocate(f(count_f))
        f = work(:count_f)
    end function factors
end module matlab_number_theory
