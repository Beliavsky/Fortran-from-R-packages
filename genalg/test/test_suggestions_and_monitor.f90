program test_suggestions_and_monitor
    use genalg, only : dp, rbga_control, rbga_result, rbga
    implicit none
    type(rbga_control) :: ctl
    type(rbga_result) :: res
    real(dp) :: lo(3), hi(3), suggestions(2,3)
    integer :: monitor_calls

    lo = -5.0_dp
    hi = 5.0_dp
    suggestions(1,:) = [1.0_dp, -2.0_dp, 0.5_dp]
    suggestions(2,:) = [3.0_dp, 3.0_dp, 3.0_dp]
    ctl%pop_size = 30
    ctl%iters = 12
    ctl%elitism = 3
    ctl%mutation_chance = 0.1_dp
    ctl%seed = 4444
    monitor_calls = 0

    call rbga(lo, hi, objective, res, ctl, suggestions, monitor)

    if (abs(res%best_value) > 1.0e-14_dp) error stop "elite suggestion was not retained"
    if (monitor_calls /= ctl%iters) error stop "monitor callback count mismatch"
    if (abs(res%best(1)) > 1.0e-14_dp) error stop "suggestion not evaluated in first generation"
    print *, "test_suggestions_and_monitor: PASS"

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = sum((x - [1.0_dp, -2.0_dp, 0.5_dp])**2)
    end function objective

    subroutine monitor(iter, population, evaluations, best, mean)
        integer, intent(in) :: iter
        real(dp), intent(in) :: population(:,:), evaluations(:), best(:), mean(:)
        if (size(population,1) /= ctl%pop_size) error stop "monitor population size"
        if (size(evaluations) /= ctl%pop_size) error stop "monitor evaluation size"
        if (size(best) /= iter .or. size(mean) /= iter) error stop "monitor history size"
        monitor_calls = monitor_calls + 1
    end subroutine monitor

end program test_suggestions_and_monitor
