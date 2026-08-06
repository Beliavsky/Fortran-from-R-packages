! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2 of the License, or any later version.
module deoptimr_interfaces
    use deoptimr_kinds, only: dp
    implicit none
    private

    public :: objective_function, constraint_function

    abstract interface
        function objective_function(x) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function objective_function

        subroutine constraint_function(x, values)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: values(:)
        end subroutine constraint_function
    end interface
end module deoptimr_interfaces
