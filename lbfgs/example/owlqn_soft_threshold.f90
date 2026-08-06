program owlqn_soft_threshold_example
    use lbfgs
    implicit none

    real(dp), parameter :: target(6) = [3.0_dp, -1.0_dp, 0.2_dp, &
        -4.0_dp, 0.0_dp, 1.5_dp]
    real(dp) :: x(6)
    type(lbfgs_parameter_t) :: parameters
    type(lbfgs_result_t) :: result

    x = 0.0_dp
    parameters = lbfgs_parameter_t()
    parameters%orthantwise_c = 0.5_dp
    parameters%linesearch = lbfgs_linesearch_backtracking
    parameters%epsilon = 1.0e-10_dp

    call lbfgs_minimize(squared_error, x, result, parameters)

    print '(a,6(1x,f9.4))', 'target:  ', target
    print '(a,6(1x,f9.4))', 'solution:', x
    print '(a)', 'The exact solution is sign(target) * max(abs(target) - 0.5, 0).'
    print '(a,1x,a)', 'status:', result%message

contains

    subroutine squared_error(x, f, g, step, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(:)
        real(dp), intent(in) :: step
        class(*), intent(inout), optional :: user_data

        g = x - target
        f = 0.5_dp * dot_product(g, g)
    end subroutine squared_error

end program owlqn_soft_threshold_example
