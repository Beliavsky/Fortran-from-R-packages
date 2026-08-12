program test_binary_ga
    use genalg, only : dp, rbga_bin_control, rbga_bin_result, rbga_bin
    implicit none
    type(rbga_bin_control) :: ctl
    type(rbga_bin_result) :: res
    integer :: i

    ctl%pop_size = 140
    ctl%iters = 100
    ctl%mutation_chance = 0.03_dp
    ctl%zero_to_one_ratio = 1.0_dp
    ctl%seed = 99173

    call rbga_bin(30, objective, res, ctl)

    if (sum(res%best_chromosome) < 29) error stop "binary GA did not approach OneMax"
    if (any(res%population < 0) .or. any(res%population > 1)) error stop "non-binary gene"
    do i = 1, size(res%population,1)
        if (sum(res%population(i,:)) < 0) error stop "impossible population sum"
    end do
    print *, "test_binary_ga: PASS", sum(res%best_chromosome), res%best_value

contains

    function objective(x) result(f)
        integer, intent(in) :: x(:)
        real(dp) :: f
        f = real(size(x) - sum(x), dp)
    end function objective

end program test_binary_ga
