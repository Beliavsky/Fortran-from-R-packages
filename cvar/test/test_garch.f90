! SPDX-License-Identifier: GPL-2.0-or-later
program test_garch
    use, intrinsic :: iso_fortran_env, only : int64
    use cvar, only : dp, cvar_ok, innovation_std_t, innovation_ged, &
                     garch11_model, garch11_simulation, garch11_forecast, &
                     make_garch11, simulate_garch11, forecast_garch11, &
                     innovation_quantile
    implicit none

    type(garch11_model) :: model, t_model, ged_model
    type(garch11_simulation) :: sim1, sim2, sim_t, sim_ged
    type(garch11_forecast) :: forecast
    real(dp) :: eta_mean, eta_variance
    integer :: i

    model = make_garch11(omega=0.4_dp, alpha=0.3_dp, beta=0.5_dp)
    if (.not. model%valid()) error stop "valid normal GARCH rejected"
    call assert_close(model%unconditional_variance(), 2.0_dp, 1.0e-14_dp, &
                      "unconditional variance")

    call simulate_garch11(model, 2000, sim1, burnin=100, seed=1234_int64)
    call assert_status(sim1%status, cvar_ok, "normal simulation status")
    call assert_close(sim1%eps(1), sqrt(sim1%h(1)) * sim1%eta(1), &
                      1.0e-14_dp, "simulation observation identity")
    do i = 2, size(sim1%eps)
        call assert_close(sim1%h(i), model%omega + model%alpha * sim1%eps(i - 1)**2 + &
                          model%beta * sim1%h(i - 1), 2.0e-13_dp, "GARCH recurrence")
    end do

    call simulate_garch11(model, 2000, sim2, burnin=100, seed=1234_int64)
    call assert_status(sim2%status, cvar_ok, "repeat simulation status")
    if (maxval(abs(sim1%eps - sim2%eps)) > 0.0_dp) &
        error stop "same seed was not reproducible"
    if (maxval(abs(sim1%h - sim2%h)) > 0.0_dp) &
        error stop "same seed variance was not reproducible"

    eta_mean = sum(sim1%eta) / real(size(sim1%eta), dp)
    eta_variance = sum((sim1%eta - eta_mean)**2) / real(size(sim1%eta) - 1, dp)
    if (abs(eta_mean) > 0.08_dp) error stop "normal innovation mean is implausible"
    if (abs(eta_variance - 1.0_dp) > 0.12_dp) error stop "normal innovation variance is implausible"

    call forecast_garch11(model, sim1%eps, sim1%h, 5, forecast, nsim=4000, &
                          seed=9876_int64)
    call assert_status(forecast%status, cvar_ok, "forecast status")
    call assert_close(forecast%h(1), model%omega + model%alpha * sim1%eps(2000)**2 + &
                      model%beta * sim1%h(2000), 2.0e-14_dp, "one-step variance forecast")
    do i = 2, 5
        call assert_close(forecast%h(i), model%omega + (model%alpha + model%beta) * &
                          forecast%h(i - 1), 2.0e-14_dp, "multi-step variance forecast")
    end do
    call assert_close(forecast%plugin_interval(1, 1), &
                      innovation_quantile(model, 0.025_dp) * sqrt(forecast%h(1)), &
                      2.0e-14_dp, "plugin lower interval")
    if (any(forecast%simulation_interval(:, 1) >= forecast%simulation_interval(:, 2))) &
        error stop "simulation interval ordering failed"

    t_model = make_garch11(0.1_dp, 0.08_dp, 0.88_dp, &
                           innovation=innovation_std_t, shape=5.0_dp)
    call simulate_garch11(t_model, 30000, sim_t, burnin=500, seed=314159_int64)
    call assert_status(sim_t%status, cvar_ok, "standardized-t simulation status")
    eta_mean = sum(sim_t%eta) / real(size(sim_t%eta), dp)
    eta_variance = sum((sim_t%eta - eta_mean)**2) / real(size(sim_t%eta) - 1, dp)
    if (abs(eta_mean) > 0.04_dp) error stop "standardized-t mean is implausible"
    if (abs(eta_variance - 1.0_dp) > 0.08_dp) error stop "standardized-t variance is implausible"

    ged_model = make_garch11(0.1_dp, 0.08_dp, 0.88_dp, &
                             innovation=innovation_ged, shape=1.5_dp)
    call simulate_garch11(ged_model, 30000, sim_ged, burnin=500, seed=271828_int64)
    call assert_status(sim_ged%status, cvar_ok, "GED simulation status")
    eta_mean = sum(sim_ged%eta) / real(size(sim_ged%eta), dp)
    eta_variance = sum((sim_ged%eta - eta_mean)**2) / real(size(sim_ged%eta) - 1, dp)
    if (abs(eta_mean) > 0.04_dp) error stop "GED mean is implausible"
    if (abs(eta_variance - 1.0_dp) > 0.08_dp) error stop "GED variance is implausible"

    model = make_garch11(-0.1_dp, 0.1_dp, 0.8_dp)
    if (model%valid()) error stop "negative omega accepted"
    model = make_garch11(0.1_dp, 0.3_dp, 0.8_dp)
    if (model%valid()) error stop "nonstationary model accepted"

    print '(a)', "test_garch: all tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        if (abs(actual - expected) > tolerance) then
            print '(a)', trim(label)
            print '(a,es24.16)', "actual:   ", actual
            print '(a,es24.16)', "expected: ", expected
            error stop "assert_close failed"
        end if
    end subroutine assert_close

    subroutine assert_status(actual, expected, label)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        if (actual /= expected) then
            print '(a,2(1x,i0))', trim(label), actual, expected
            error stop "assert_status failed"
        end if
    end subroutine assert_status

end program test_garch
