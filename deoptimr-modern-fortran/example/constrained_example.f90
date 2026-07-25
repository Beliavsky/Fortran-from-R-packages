! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Portions derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This program is free software: you may redistribute it and/or modify it
! under the terms of GNU GPL version 2, or (at your option) any later version.
module constrained_functions
    use deoptimr, only: dp
    implicit none
contains
    function objective(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = sum(x*x)
    end function objective

    subroutine equality(x, values)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        values(1) = x(1) + x(2) - 1.0_dp
    end subroutine equality
end module constrained_functions

program constrained_example
    use deoptimr, only: dp, seed_rng, jde_control, de_result, jde_optimize
    use constrained_functions, only: objective, equality
    implicit none

    type(jde_control) :: control
    type(de_result) :: result
    real(dp) :: equality_tolerance(1)

    call seed_rng(7411)
    control%population_size = 70
    control%max_iterations = 1800
    control%tolerance = 1.0e-8_dp
    equality_tolerance = 1.0e-5_dp

    call jde_optimize([-2.0_dp, -2.0_dp], [2.0_dp, 2.0_dp], objective, result, &
        control, equality, 1, 1, equality_tolerance)

    write(*, '(a,2(1x,es16.8))') 'parameters:', result%parameters
    write(*, '(a,1x,es16.8)') 'objective:', result%value
    write(*, '(a,1x,es16.8)') 'equality residual:', sum(result%parameters) - 1.0_dp
end program constrained_example
