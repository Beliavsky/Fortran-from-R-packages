program test_single_gene
    use genalg, only : dp, rbga_control, rbga_bin_control, rbga_result, rbga_bin_result, rbga, rbga_bin
    implicit none
    type(rbga_control) :: rctl
    type(rbga_bin_control) :: bctl
    type(rbga_result) :: rres
    type(rbga_bin_result) :: bres

    rctl%pop_size = 40
    rctl%iters = 30
    rctl%elitism = 4
    rctl%mutation_chance = 0.15_dp
    rctl%seed = 1234
    call rbga([-2.0_dp], [3.0_dp], real_objective, rres, rctl)
    if (rres%best_value > 5.0e-2_dp) error stop "single-gene real GA failed"

    bctl%pop_size = 30
    bctl%iters = 10
    bctl%elitism = 3
    bctl%mutation_chance = 0.0_dp
    bctl%seed = 9876
    bctl%legacy_binary_eval_cache = .false.
    call rbga_bin(1, binary_objective, bres, bctl)
    if (abs(bres%best_value) > 1.0e-15_dp) error stop "single-gene binary GA failed"
    print *, "test_single_gene: PASS"

contains

    function real_objective(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = abs(x(1) - 0.75_dp)
    end function real_objective

    function binary_objective(x) result(f)
        integer, intent(in) :: x(:)
        real(dp) :: f
        f = real(1 - x(1), dp)
    end function binary_objective

end program test_single_gene
