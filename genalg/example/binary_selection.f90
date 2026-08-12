program binary_selection
    use genalg, only : dp, rbga_bin_control, rbga_bin_result, rbga_bin
    implicit none
    type(rbga_bin_control) :: control
    type(rbga_bin_result) :: result

    control%pop_size = 120
    control%iters = 80
    control%mutation_chance = 0.02_dp
    control%zero_to_one_ratio = 2.0_dp
    control%seed = 54321
    ! Use corrected fitness caching for new applications.
    control%legacy_binary_eval_cache = .false.

    call rbga_bin(20, evaluate, result, control)
    print '(a,20i2)', "best chromosome: ", result%best_chromosome
    print '(a,es14.6)', "objective:       ", result%best_value

contains
    function evaluate(x) result(f)
        integer, intent(in) :: x(:)
        real(dp) :: f
        integer, parameter :: target(20) = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
        f = real(sum(abs(x-target)), dp)
    end function evaluate
end program binary_selection
