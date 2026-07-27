! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Portions derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This program is free software: you may redistribute it and/or modify it
! under the terms of GNU GPL version 2, or (at your option) any later version.
module demo_objectives
    use deoptimr, only: dp
    implicit none
contains
    function sphere(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = sum(x*x)
    end function sphere

    function aluffi(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = 0.25_dp*x(1)**4 - 0.5_dp*x(1)**2 + 0.1_dp*x(1) + 0.5_dp*x(2)**2
    end function aluffi

    function becker_lago(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = sum((abs(x) - 5.0_dp)**2)
    end function becker_lago
end module demo_objectives

program demo_deoptimr
    use deoptimr, only: dp, seed_rng, jde_control, ncde_control, de_result, ncde_result
    use deoptimr, only: jde_optimize, spjde_optimize, ncde_optimize
    use demo_objectives, only: sphere, aluffi, becker_lago
    implicit none

    type(jde_control) :: jde_cfg
    type(ncde_control) :: ncde_cfg
    type(de_result) :: jde_result, spjde_result
    type(ncde_result) :: ncde_result_value

    call seed_rng(2345)
    jde_cfg%population_size = 50
    jde_cfg%max_iterations = 900
    jde_cfg%tolerance = 1.0e-9_dp
    call jde_optimize([-100.0_dp, -100.0_dp], [100.0_dp, 100.0_dp], &
        sphere, jde_result, jde_cfg)

    call seed_rng(42)
    jde_cfg%population_size = 40
    jde_cfg%max_iterations = 1600
    jde_cfg%cr_upper = 1.1_dp
    call spjde_optimize([-10.0_dp, -10.0_dp], [10.0_dp, 10.0_dp], &
        aluffi, spjde_result, jde_cfg)

    call seed_rng(2345)
    ncde_cfg%population_size = 120
    ncde_cfg%niche_radius = 4.0_dp
    ncde_cfg%archive_size = 20
    ncde_cfg%max_iterations = 350
    ncde_cfg%neighbor_lower = 6
    ncde_cfg%neighbor_upper = 24
    call ncde_optimize([-10.0_dp, -10.0_dp], [10.0_dp, 10.0_dp], &
        becker_lago, ncde_result_value, ncde_cfg)

    write(*, '(a,2(1x,es14.6),1x,a,1x,es14.6)') &
        'JDE sphere:', jde_result%parameters, 'f=', jde_result%value
    write(*, '(a,2(1x,es14.6),1x,a,1x,es14.6)') &
        'SPJDE Aluffi:', spjde_result%parameters, 'f=', spjde_result%value
    write(*, '(a,1x,i0)') 'NCDE archived solutions:', size(ncde_result_value%objective_archive)
end program demo_deoptimr
