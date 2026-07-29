! SPDX-License-Identifier: GPL-2.0-or-later
program garch_forecast
    use, intrinsic :: iso_fortran_env, only : int64
    use cvar, only : dp, cvar_ok, garch11_model, garch11_simulation, &
                     garch11_forecast, make_garch11, simulate_garch11, &
                     forecast_garch11
    implicit none

    type(garch11_model) :: model
    type(garch11_simulation) :: simulation
    type(garch11_forecast) :: forecast
    integer :: i

    model = make_garch11(omega=0.4_dp, alpha=0.3_dp, beta=0.5_dp)
    call simulate_garch11(model, 1000, simulation, burnin=100, seed=20220305_int64)
    if (simulation%status /= cvar_ok) error stop "simulation failed"

    call forecast_garch11(model, simulation%eps, simulation%h, 5, forecast, &
                          nsim=10000, seed=1234_int64)
    if (forecast%status /= cvar_ok) error stop "forecast failed"

    print '(a)', "GARCH(1,1) volatility and 95% predictive intervals"
    print '(a)', " horizon      E[h]       plugin lower  plugin upper  simulated lower  simulated upper"
    do i = 1, 5
        print '(i5,5(2x,f13.6))', i, forecast%h(i), &
            forecast%plugin_interval(i, 1), forecast%plugin_interval(i, 2), &
            forecast%simulation_interval(i, 1), forecast%simulation_interval(i, 2)
    end do
end program garch_forecast
