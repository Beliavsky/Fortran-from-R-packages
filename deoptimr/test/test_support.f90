! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Portions derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This program is free software: you may redistribute it and/or modify it
! under the terms of GNU GPL version 2, or (at your option) any later version.
module test_support
    use deoptimr_kinds, only: dp
    implicit none
    private
    public :: assert_true, assert_close, assert_all_close
contains
    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            write(*, '(a)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message
        if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) then
            write(*, '(a,2(1x,es24.16))') 'FAIL: '//trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_all_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual(:), expected(:), tolerance
        character(len=*), intent(in) :: message
        if (size(actual) /= size(expected)) then
            write(*, '(a)') 'FAIL: '//trim(message)//' size mismatch'
            error stop 1
        end if
        if (any(abs(actual - expected) > tolerance*max(1.0_dp, abs(expected)))) then
            write(*, '(a)') 'FAIL: '//trim(message)
            write(*, '(a,*(1x,es16.8))') 'actual:', actual
            write(*, '(a,*(1x,es16.8))') 'expected:', expected
            error stop 1
        end if
    end subroutine assert_all_close
end module test_support
