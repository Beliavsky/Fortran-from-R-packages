program pi_sqrt50
    use genalg, only : dp, rbga_control, rbga_result, rbga
    implicit none
    type(rbga_control) :: control
    type(rbga_result) :: result

    control%pop_size = 200
    control%iters = 100
    control%mutation_chance = 0.01_dp
    control%seed = 12345
    call rbga([1.0_dp, 1.0_dp], [5.0_dp, 10.0_dp], evaluate, result, control)

    print '(a,2f14.8)', "best chromosome: ", result%best_chromosome
    print '(a,es14.6)', "objective:       ", result%best_value

contains
    function evaluate(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = abs(x(1)-acos(-1.0_dp)) + abs(x(2)-sqrt(50.0_dp))
    end function evaluate
end program pi_sqrt50
