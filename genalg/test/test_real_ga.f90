program test_real_ga
    use genalg, only : dp, rbga_control, rbga_result, rbga
    implicit none
    type(rbga_control) :: ctl
    type(rbga_result) :: res
    real(dp) :: lo(2), hi(2)

    lo = [1.0_dp, 1.0_dp]
    hi = [5.0_dp, 10.0_dp]
    ctl%pop_size = 160
    ctl%iters = 120
    ctl%mutation_chance = 0.02_dp
    ctl%seed = 11731

    call rbga(lo, hi, objective, res, ctl)

    if (res%best_value > 5.0e-2_dp) error stop "real GA did not converge"
    if (any(res%population(:,1) < lo(1)) .or. any(res%population(:,1) > hi(1))) then
        error stop "real GA violated first bound"
    end if
    if (any(res%population(:,2) < lo(2)) .or. any(res%population(:,2) > hi(2))) then
        error stop "real GA violated second bound"
    end if
    if (res%nfe < ctl%pop_size) error stop "invalid evaluation count"
    print *, "test_real_ga: PASS", res%best_value

contains

    function objective(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = abs(x(1) - acos(-1.0_dp)) + abs(x(2) - sqrt(50.0_dp))
    end function objective

end program test_real_ga
