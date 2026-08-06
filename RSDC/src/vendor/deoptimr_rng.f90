! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2 of the License, or any later version.
module deoptimr_rng
    use deoptimr_kinds, only: dp
    implicit none
    private

    public :: seed_rng, random_uniform, random_integer, sample_without_replacement

contains

    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729*i + 8191*i*i, huge(1) - 1)
            if (put(i) == 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine seed_rng

    function random_uniform(lower, upper) result(value)
        real(dp), intent(in), optional :: lower, upper
        real(dp) :: value, lo, hi

        lo = 0.0_dp
        hi = 1.0_dp
        if (present(lower)) lo = lower
        if (present(upper)) hi = upper
        call random_number(value)
        value = lo + (hi - lo)*value
    end function random_uniform

    function random_integer(lower, upper) result(value)
        integer, intent(in) :: lower, upper
        integer :: value
        real(dp) :: u

        if (upper < lower) error stop "random_integer: invalid interval"
        call random_number(u)
        value = lower + int(u*real(upper - lower + 1, dp))
        if (value > upper) value = upper
    end function random_integer

    subroutine sample_without_replacement(candidates, k, sample)
        integer, intent(in) :: candidates(:)
        integer, intent(in) :: k
        integer, intent(out) :: sample(k)
        integer, allocatable :: work(:)
        integer :: i, j, tmp, n

        n = size(candidates)
        if (k < 0 .or. k > n) error stop "sample_without_replacement: invalid sample size"
        allocate(work(n))
        work = candidates
        do i = 1, k
            j = random_integer(i, n)
            tmp = work(i)
            work(i) = work(j)
            work(j) = tmp
            sample(i) = work(i)
        end do
    end subroutine sample_without_replacement
end module deoptimr_rng
