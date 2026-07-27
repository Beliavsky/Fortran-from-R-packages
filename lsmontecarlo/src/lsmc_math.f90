! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_math
    use lsmc_kinds, only : dp
    implicit none
    private

    real(dp), parameter, public :: pi = acos(-1.0_dp)

    public :: normal_cdf
    public :: mean_value
    public :: sample_standard_error

contains

    pure elemental function normal_cdf(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value

        value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
    end function normal_cdf

    pure function mean_value(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value

        if (size(x) == 0) then
            value = 0.0_dp
        else
            value = sum(x) / real(size(x), dp)
        end if
    end function mean_value

    pure function sample_standard_error(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        real(dp) :: avg
        integer :: n

        n = size(x)
        if (n <= 1) then
            value = 0.0_dp
            return
        end if

        avg = sum(x) / real(n, dp)
        value = sqrt(sum((x - avg)**2) / real(n - 1, dp) / real(n, dp))
    end function sample_standard_error

end module lsmc_math
