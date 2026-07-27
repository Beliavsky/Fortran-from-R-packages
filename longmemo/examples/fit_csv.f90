! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

program fit_csv
    use longmemo_kinds, only : dp
    use longmemo_io, only : read_index_value_csv
    use longmemo, only : whittle_result, fexp_result, whittle_estimate, fexp_estimate
    implicit none

    character(len=512) :: filename
    real(dp), allocatable :: x(:)
    type(whittle_result) :: whittle
    type(fexp_result) :: fexp

    if (command_argument_count() /= 1) then
        print '(a)', "usage: fit_csv FILE.csv"
        stop 1
    end if
    call get_command_argument(1, filename)
    call read_index_value_csv(trim(filename), x)

    call whittle_estimate(x, "fGn", whittle, start_eta=[0.6_dp], covariance_m=4096)
    call fexp_estimate(x, 3, 0.5_dp, fexp)

    print '(a,i0)', "n = ", size(x)
    print '(a,f12.6)', "Whittle H = ", whittle%eta(1)
    print '(a,f12.6)', "Whittle SE = ", whittle%std_error(1)
    print '(a,f12.6)', "FEXP H = ", fexp%hurst
    print '(a,i0)', "FEXP order = ", fexp%order_poly
end program fit_csv
